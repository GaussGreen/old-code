//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CreditFeeLegWithPV.hpp
//
//   Description : ICreditFeeLeg able to pv itself
//
//   Author      : Antoine Gregoire
//
//   Date        : 3 Nov 2004
//
//----------------------------------------------------------------------------

#ifndef CREDIT_FEE_LEG_WITH_PV_HPP
#define CREDIT_FEE_LEG_WITH_PV_HPP

#include "edginc/ICreditFeeLeg.hpp"

DRLIB_BEGIN_NAMESPACE


/** base class for credit fee legs that implement pv */
class PRODUCTS_DLL CreditFeeLegWithPV : public CObject, public virtual ICreditFeeLeg{
public :

	/** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Return the value of this leg */
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

	/** values all  cashflows after cfCutOffDate and before or on endDate*/
	/** value of risky cashflows only */
	virtual double pvRisky(const DateTimeArray & timeline,
		const DoubleArray & effectiveCurve,
		const DateTime& cfCutOffDate, // ignore all cashflows on an before this date
		const DateTime & endDate,
		YieldCurveConstSP discCurve,
        IForwardRatePricerSP model);

	// overloaded version of pvRisky without end date
	double pvRisky(const DateTimeArray& timeline,
                   const DoubleArray&   effCurve,
                   const DateTime&      cfCutOffDate,
                   YieldCurveConstSP    discCurve,
                   IForwardRatePricerSP model);

    /** Returns the risky or riskless abstract cash flows, depending on the
     * value of the "wantRisky" parameter */
    AbstractCashFlowArrayConstSP getAbstractCashFlows(bool wantRisky,
                                                      IForwardRatePricerSP model) const;

    /** Returns risky cash flows only */
    virtual CashFlowArraySP getRiskyCashFlows(IForwardRatePricerSP model) const;

    /** Returns risk free cash flows only */
    virtual CashFlowArraySP getRisklessCashFlows(IForwardRatePricerSP model) const;

    /** Returns risky notional dates */
    virtual DateTimeArraySP getRiskyNotionalDates(IForwardRatePricerSP model) const;

    /** Return known cash flows corresponding to a non-defaulted CDS */
    virtual CashFlowArraySP generateKnownCashFlows(IForwardRatePricerSP model) const;

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

	static CreditFeeLegWithPV * buildCreditFeeLeg(
        double                     coupon,
        const DateTimeArray&       observationDates, // observation dates
        const DateTimeArray&       accrualStartDates, // accrual dates : can be non contiguous
        const DateTimeArray&       accrualEndDates,
        const DateTimeArray&       payDates,
        const CouponNotionalTypes& couponNotionalType,
        const string&              payDCC,
        const IModel*              model,
        const MarketData*          market);

    /** Estimates the known cash flows (so they are not really "known")
     * corresponding to a CDO tranche - takes into account estimated losses
     * in the future */
    virtual CashFlowArraySP estimateKnownCashFlows(
         const double         initialTrancheSize,
         const DateTimeArray  pastTrancheLossDates,
         const DoubleArray    pastTrancheLosses,
         IForwardRatePricerSP model);

    /** Returns the accrued interest of a non-defaulted CDS*/
    virtual double getFeeLegAI(const DateTime&              valuationDate, 
                               const DateTime&              paymentDate, 
                               const DateTime&              earliestAccrualStart,
                               const DateTime&              latestAccrualEnd,
                               const DayCountConventionSP   dcc, //allows an override to be specified
                               const IDiscountCurveRisky&   crv,
                               const IDiscountCurveConstSP  discount,
                               const IDecretionCurveConstSP prepay,
                               IForwardRatePricerSP         model) const;

    /* Gets the value of the fee leg corresponging to a defaulted CDS
       under the (optional) credit event */
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

    virtual double getFeeLegPV(const DateTime&              valuationDate, 
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

protected:
	CreditFeeLegWithPV(CClassConstSP clazz);

private:
	static void load(CClassSP& clazz);

    /** Prices the rebate cashflows by computing the pay date of each of them,
     * according to the "pay as you go" configuration */
    double priceRebatePayments(
        const IDiscountCurveRiskySP effectiveCurve,
        const DateTime&             valDateCF,
        CashFlowArraySP             rebatePayments,      //
        BoolArrayConstSP            payAsYouGoArray,     // These params are
        IntArrayConstSP             numDelayDaysArray,   // only required if
        DateTimeArrayConstSP        startDates,          // there are rebate
        DateTimeArrayConstSP        endDates,            // payments. Otherwise
        DateTimeArrayConstSP        paymentDates,        // NULL is fine
        IBadDayAdjusterConstSP      bda) const;

    /** Returns the determination and accrual pay date corresponding
        to a defaulted CDS, using the event override if required */
    void getDeterminationAndPayDate(
        const DateTime&            protectionStartDate,
        const DateTime&            protectionEndDate,
        const DateTime&            valuationDate,
        const DateTime&            defaultDate,
        ICreditEventOverrideNameSP creditEventOverride,
        const DateTime&            lastTriggerDate,
        CIntSP                     triggerDelay,
        CIntSP                     defaultToSettlementDelay,
        IBadDayAdjusterConstSP     badDayAdjuster,
        DateTime&                  determinationDate,     // (O)
        DateTime&                  accrualPayDate) const; // (O)
};

typedef smartPtr<CreditFeeLegWithPV> CreditFeeLegWithPVSP;
typedef smartConstPtr<CreditFeeLegWithPV> CreditFeeLegWithPVConstSP;

DRLIB_END_NAMESPACE

#endif // CREDIT_FEE_LEG_WITH_PV_HPP
