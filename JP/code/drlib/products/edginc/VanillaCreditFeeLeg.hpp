//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : VanillaCreditFeeLeg.hpp
//
//   Description : A Fee Leg for a Vanilla Credit Instrument
//
//   Author      : Doris Morris
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#ifndef VANILLA_CREDIT_FEE_LEG_HPP
#define VANILLA_CREDIT_FEE_LEG_HPP


#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
#include "edginc/ICreditFeeLeg.hpp"
#include "edginc/IFixedRateCreditFeeLeg.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Stub.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(BadDayConvention);
FORWARD_DECLARE(DayCountConvention);
FORWARD_DECLARE(IForwardRatePricer);
FORWARD_DECLARE_WRAPPER(Holiday);

/* Base Fee Leg class for a Vanilla Credit Instrument*/
class PRODUCTS_DLL VanillaCreditFeeLeg : public CObject,
                                         public virtual ICreditFeeLeg,
                                         public virtual ICreditFeeLegGenerator,
                                         public virtual IFixedRateCreditFeeLeg,
                                         public virtual Theta::Shift
{
public:

    //----------------------------
    // VanillaCreditFeeLeg methods
    //----------------------------

    /* TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /* Constructor */
    VanillaCreditFeeLeg(double                      baseNotional,  // initial notional
                        const  DateTime&            valueDate,     // when to value the leg as of
                        const  DateTime&            effectiveDate, // protection start date
                        const  DateTime&            maturityDate,  // protection end date
                        double                      feeRate,       // fixed fee spread
                        const  MaturityPeriodSP     paymentFreq,   // payment frequency
                        const  DateTime&            fcd,           // first full coupon
                        const  DateTime&            lcd,           // last full coupon
                        const  DayCountConventionSP pmtDcc,        // payment day count convention
                        const  DayCountConventionSP accrualDcc,    // accrual day count convention
                        const  BadDayConventionSP   pmtBdc,        // payment bad day convention
                        const  BadDayConventionSP   accrualBdc,    // accrual bad day convention
                        const  BadDayConventionSP   valuationBdc,  // valuation bad day convention
                        const  StubSP               stubType,      // Stub payment types: none, simple, bond
                        const  HolidayWrapper       pmtHol,        // holiday calendar for payments
                        const  HolidayWrapper       accrualHol,    // accrual holidays
                        const  YieldCurveWrapper    discountCrv,   // curve to discount fee payments
                        bool                        payAccruedFee,
                        double                      delay);

    virtual ~VanillaCreditFeeLeg();

    /**
    * Accrued interest (cash amount, not percentage) for settlement on settlementDate.
    * Assumes fee payments are in chronological order.
    */
    double getAccruedInterest(const DateTime& settlementDate) const;

    /**Set the base notional of an instrument. This is here because some derivatives
       need to be able to control the notional of their underlyings to guarantee
       that the valuations make sense. In general the instrument will be loaded with
       the correct notional, so this method will not be commonly used.*/
    void setNotional( double newNotional);

    //----------------
    // CObject methods
    //----------------

    /** validate object inputs */
    virtual void validatePop2Object();

    /** Called once before the initial pricing */
    virtual void Validate();

    //----------------------
    // ICreditFeeLeg methods
    //----------------------

    virtual double price(
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
        IForwardRatePricerSP        model) const;

    /** Return all cash flow dates */
    virtual DateTimeArraySP getCashFlowDates() const;

    /** Return risky cash flow dates */
    virtual DateTimeArraySP getRiskyCashFlowDates() const;

    /** Return riskfree cash flow dates */
    virtual DateTimeArraySP getRisklessCashFlowDates() const;

    /** Returns the accrual periods (start and end dates) in this fee leg. */
    virtual AccrualPeriodArrayConstSP getAccrualPeriods() const;

    /** Returns risky accrual periods (start and end dates) in this fee leg. */
    virtual AccrualPeriodArrayConstSP getRiskyAccrualPeriods() const;

    /** Returns risky accrual periods (start and end dates) in this fee leg. */
    virtual AccrualPeriodArrayConstSP getRisklessAccrualPeriods() const;

    /** Returns all cash flows */
    virtual AbstractCashFlowArrayConstSP getCashFlows(IForwardRatePricerSP model) const;

    /** Returns risky cash flows only */
    virtual CashFlowArraySP getRiskyCashFlows(IForwardRatePricerSP model) const;

    /** Returns risk free cash flows only */
    virtual CashFlowArraySP getRisklessCashFlows(IForwardRatePricerSP model) const;

    /** Returns risky notional dates */
    virtual DateTimeArraySP getRiskyNotionalDates(IForwardRatePricerSP model) const;

	/** Returns risky coupon notional types */
	virtual CouponNotionalTypesArraySP getRiskyCouponNotionalTypes() const;

	/** Returns risky observation dates */
	virtual DateTimeArraySP getRiskyObservationDates() const;
   
	/** Return known cash flows corresponding to a CDO tranche */
    virtual CashFlowArraySP generateKnownCashFlows(
         const DateTime       today,
         const double         initialTrancheSize,
         const DateTimeArray  pastTrancheLossDates,
         const DoubleArray    pastTrancheLosses,
         const double         pastTrancheLoss,
         IForwardRatePricerSP model);

    /** Return known cash flows corresponding to a non-defaulted CDS */
    virtual CashFlowArraySP generateKnownCashFlows(
         IForwardRatePricerSP model) const;

    /** Estimates the known cash flows (so they are not really "known")
     * corresponding to a CDO tranche - takes into account estimated losses
     * in the future */
    virtual CashFlowArraySP estimateKnownCashFlows(
         const double         initialTrancheSize,
         const DateTimeArray  pastTrancheLossDates,
         const DoubleArray    pastTrancheLosses,
         IForwardRatePricerSP model);

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
    virtual CashFlowArraySP generateKnownCashFlowsGivenDefault(
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
        const DateTime&            lastTriggerDate) const;

    /** Returns the pay date which is last in terms of time */
    virtual DateTime getLastPayDate() const;

    /** Returns the observation date which is last in terms of time */
    virtual DateTime getLastObservationDate() const;

    /** When to stop tweaking for Yield Curve type tweaks */
    virtual DateTime lastYCSensDate(const DateTime& currentLastDate) const;

    /** Feedback method for getting information about a fee cashflow */
    virtual void getActiveFee(const DateTime&      withRespectTo,       // (I) get the fee whose accrual period contains this date
                              const DateTime&      earliestAccrualDate, // (I) for fee legs that dont specify accrue start dates, and we are interested in the first fee
                              IForwardRatePricerSP model,               // (I) for calculating the amount
                              DateTime&            accrueStartDate,     // (O) when the fee starts accruing
                              DateTime&            accrueEndDate,       // (O) when the fee finishes accruing
                              DateTime&            paymentDate,         // (O) when the fee is paid
                              double&              amount) const;       // (O) the cashflow amount

    /** Price this fee leg assuming it corresponds to a defaulted CDS */
    virtual double pv(const DateTime&         valueDate,
                      const DateTime&         defaultDeterminationDate,
                      const DateTime&         accrualPaymentDate,
                      const YieldCurveConstSP discount,
                      const bool              computeAccrual,
                      IForwardRatePricerSP    model) const;

    /** Retrieve market data for this fee leg */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Compute PV of fee leg at valuation date. */
    virtual double getFeeLegPV(const DateTime&              valuationDate, 
                               const DateTime&              earliestRiskyDate,
                               const DateTime&              latestRiskyDate,
                               const IDiscountCurve&        discount,
                               const IDiscountCurveRisky&   crv,
                               const IDecretionCurveConstSP prepay,
                               const bool                   includeAccrued,
                               const DayCountConventionSP   dcc,
                               IForwardRatePricerSP         model) const;

    /**
    * Compute PV of fee leg at valuationDate, but for unconditional payment at
    * paymentDate, conditional on no new defaults before valuationDate.
    */
    virtual double getFeeLegPV(const DateTime&              valuationDate, 
                               const DateTime&              paymentDate, 
                               const DateTime&              earliestRiskyDate,
                               const DateTime&              latestRiskyDate,
                               const IDiscountCurve&        discount,
                               const IDiscountCurveRisky&   crv,
                               const IDecretionCurveConstSP prepay,
                               const bool                   includeAccrued,
                               const DayCountConventionSP   dcc,
                               IForwardRatePricerSP         model) const;

    virtual double getFeeLegPV(const DateTime&              valuationDate, 
                               const DateTime&              paymentDate, 
                               const DateTime&              earliestRiskyDate,
                               const DateTime&              latestRiskyDate,
                               const IDiscountCurve&        discount,
                               const IDiscountCurveRisky&   crv,
                               const IDecretionCurveConstSP prepay,
                               const bool                   includeAccrued,
                               const DayCountConventionSP   dcc,
                               bool                         defaultValueOnly,
                               IForwardRatePricerSP         model) const;

    /** Returns the accrued interest */
    virtual double getFeeLegAI(const DateTime&              valuationDate, 
                               const DateTime&              paymentDate, 
                               const DateTime&              earliestAccrualStart,
                               const DateTime&              latestAccrualEnd,
                               const DayCountConventionSP   dcc, //allows an override to be specified
                               const IDiscountCurveRisky&   crv,
                               const IDiscountCurveConstSP  discount,
                               const IDecretionCurveConstSP prepay,
                               IForwardRatePricerSP         model) const;

    /** Prices the leg sufferring a default, under the (optional) credit event */
    virtual double getFeeLegDefaultedPV(
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
        DateTime                     lastTriggerDate) const;

    /* Computes the interest accrued up to valuationDate */
    virtual double getAccruedInterest(const DateTime& valuationDate,
                                      IForwardRatePricerSP model) const;

    /** Returns the leg notional */
    virtual double getFeeLegNotional() const;

    /** Returns the leg notional */
    virtual void setFeeLegNotional(double newNotional);    

    //-------------------------------
    // ICreditFeeLegGenerator methods
    //-------------------------------
    virtual ICreditFeeLegSP generateCreditFeeLeg(const DateTime& startDate,
                                                 const DateTime& endDate,
                                                 double feeRate) const;

    //-------------------------------
    // IFixedRateCreditFeeLeg methods
    //-------------------------------

    /** Get the fee rate of the fee leg.*/
    virtual double getRate() const;

    /** Change the fee rate of the fee leg.*/
    virtual void setRate(double newRate);

    //---------------------
    // Theta::Shift methods
    //---------------------
    virtual bool sensShift(Theta* shift);

	/** creates the SV Generator
		Note this method throws exception - ASarma
	*/
	virtual ICreditFeeLegSVGenSP createSVGen(
		const DateTime& valueDate,
		const DateTimeArray& modifiedFeeLegObservationDates,
		const DateTimeArray& productTimeline, //product timeline
		const IntArray& dateToDiscFactorIndex, //map from a date to an index on the DiscountFactor SV
		double lossConfigNotional
		) const;


    /****************** ICreditVanillaInstrument ******************/
    /**Accrued interest (cash amount, not percentage) for settlement on settlementDate.*/
//    virtual double getAccruedInterest(const DateTime& settlementDate) const;

    /**Calculates PV of the instrument, given a valuation date.
       Settlement is whatever the instrument determines it to be (e.g. T+1 CDS; T+3 bond, etc.)
       relative to the valuation date. Value should be conditional on no new defaults
       before valuationDate.*/
//    virtual double getPV(const DateTime&                 valuationDate,
 //                        const IDiscountCurveRisky&      crv) const;

    /**Calculates the PV of the instrument on a given valuationDate for unconditional settlement
       and payment on paymentDate. Used to calculate forward PVs and differentiate between
       conditional and unconditional settlement. Value should be conditional on no new defaults
       before valuationDate if valuationDate is in the future.*/
//    virtual double getPV(const DateTime&                 valuationDate,
  //                       const DateTime&                 paymentDate,
    //                     const IDiscountCurveRisky&      crv) const;


    /** Is perReturns true if the instrument has a finite maturity. This will be false if
        it is perpetual.*/
//    virtual bool hasFiniteMaturity() const;

    /**Returns the maturity of the leg. */
//    virtual DateTime getMaturity() const;

    /**Returns the earliest possible "interesting" date: this will be the earliest of the
       start of the first accrual period, the start of the contingent leg, the first
       cash-flow, etc.*/
//    virtual DateTime getStartDate() const;

    /** Return the base notional. */
//    virtual double getNotional() const;


    /* Get cashflows */
//    virtual CashFlowArraySP getCashFlows() const;
    /*************************************************************/

    /************************IHasRiskyFeeLeg **********************/
    /**Compute PV of fee leg at valuation date.*/
//    virtual double getFeeLegPV(const DateTime&                 valuationDate,
  //                             const IDiscountCurveRisky&      crv) const;

    /**Compute PV of fee leg at valuationDate, but for unconditional payment at
       paymentDate, conditional on no new defaults before valuationDate.*/
//    virtual double getFeeLegPV(const DateTime&                 valuationDate,
  //                             const DateTime&                 paymentDate,
    //                           const IDiscountCurveRisky&      crv) const;

//    virtual double getFeeLegPV(const DateTime&                 valuationDate,
  //                             const DateTime&                 paymentDate,
    //                           const IDiscountCurveRisky&      crv,
      //                         bool                            defaultValueOnly) const;

    /** override a control shift (eg for delta on trees) - may return
        null to use original. */
//    virtual CSensControl* AlterControl( const IModel*          modelParams,
  //                                      const CSensControl*    sensControl) const;

    /** Returns the value date (aka today) the instrument is currently
        pricing for */
//    virtual DateTime getValueDate() const;

    /** price a dead instrument until settlement - exercised, expired,
        knocked out etc.  returns true if it is dead (and priced), false
        if it is not dead */
  //  virtual bool priceDeadInstrument(CControl* control,
    //                                 CResults* results) const;

    /** Returns the name of the instrument's discount currency. */

 //   virtual string discountYieldCurveName() const;

    /*************************************************************/


private:

    VanillaCreditFeeLeg();

    // To create a VanillaCreditFeeLeg from a small set of input data
    VanillaCreditFeeLeg(const DateTime& startDate,
                        const DateTime& endDate,
                        double feeRate);

    // Class load mechanism
    static void load(CClassSP& clazz);

    // To create a default constructor with the EMPTY_SHELL_METHOD
    static IObject* defaultConstructor();

   /* generate the fee cash flow(s) */
    virtual CashFlowArraySP generateFeePayments();
    virtual DateTime getFirstCouponDate()const;

    //-------
    // Fields
    //-------
    mutable double          baseNotional;          // initial notional
    DateTime                valueDate;             // Valuation date
    DateTime                effectiveDate;         // aka protection start date
	                                               // aka for pyramid - initial accrual date
    DateTime                maturityDate;          // Maturity date of the fee
	                                               // aka for pyramid - final accrual date
    mutable double          feeRate;               // either a fixed fee/spread or the current floating fee rate

    /* data field required for pyramid to generate coupon/pmt schedule */
    MaturityPeriodSP        paymentFreq;           // payment frequency
    mutable DateTime        fcd;                   // First full coupon date aka FirstRoll in pyramid
    DateTime                lcd;                   // Last coupon date aka LastRoll in pyramid

    DayCountConventionSP    pmtDcc;                // Payment Day Count Convention
    DayCountConventionSP    accrualDcc;            // Accrual Day Count Convention

    BadDayConventionSP      pmtBdc;                // Payment Bad Day Convention
    BadDayConventionSP      accrualBdc;            // Accrual Bad Day Convention
    BadDayConventionSP      valuationBdc;          // Valuation Bad Day Convention

    HolidayWrapper          pmtHol;                // Payment holiday
    HolidayWrapper          accrualHol;            // Accrual holiday

    YieldCurveWrapper       discount;              // discount

    bool                    payAccruedFee;        // Whether to pay accrued fee on default or not
	                                              // Accrual dcc and bdc is delegated to the fee leg
    double                  delay;                // payment offset
    StubSP                  stubPaymentType;      // stub payment type to determine accrued interest
                                                  // and first cpn amount
    /* generated results */
    mutable CashFlowArraySP feePayments;           // cashflow array of fee payments

};

typedef smartPtr<VanillaCreditFeeLeg> VanillaCreditFeeLegSP;
typedef smartConstPtr<VanillaCreditFeeLeg> VanillaCreditFeeLegConstSP;


DRLIB_END_NAMESPACE

#endif // VANILLA_CREDIT_FEE_LEG_HPP
