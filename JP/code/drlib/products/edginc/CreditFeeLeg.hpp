//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CreditFeeLeg.hpp
//
//   Description : Simple implementation of ICreditFeeLeg
//
//   Author      : Antoine Gregoire
//
//   Date        : 3 Nov 2004
//
//----------------------------------------------------------------------------

#ifndef CREDIT_FEE_LEG_HPP
#define CREDIT_FEE_LEG_HPP

#include "edginc/CreditFeeLegWithPV.hpp"
#include "edginc/Theta.hpp"
#include "edginc/CreditFeeLegCoupon.hpp"
#include "edginc/IFixedRateCreditFeeLeg.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implementation of ICreditFeeLeg without using ICreditCashFlow */
class PRODUCTS_DLL CreditFeeLeg : public CreditFeeLegWithPV,
                                virtual public IRestorableWithRespectTo<CreditFeeLegCoupon>,
                                virtual public Theta::Shift,
                                virtual public IFixedRateCreditFeeLeg //temporary
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    CreditFeeLeg();

    /** public constructor for fixed coupon and unit notional */
    CreditFeeLeg(double  coupon,
                 const DateTimeArray       &observationDates, // observation dates
                 const DateTimeArray       &accrualStartDates, // accrual dates : can be non contiguous
                 const DateTimeArray       &accrualEndDates,
                 const DateTimeArray       &payDates,
                 const CouponNotionalTypes &couponNotionalType,
                 const string              &payDCC,
                 const IModel* model, const MarketData* market);

    virtual ~CreditFeeLeg();

    /** Overrides default */
    virtual void validatePop2Object();

    /** Theta shift implementation*/
    bool sensShift(Theta* shift) ;

    /** GetMarket implementation*/
    void getMarket(const IModel* model, const MarketData* market);

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

	/** Returns risky coupon notional types */
	virtual CouponNotionalTypesArraySP getRiskyCouponNotionalTypes() const;

	/** Returns risky observation dates */
	virtual DateTimeArraySP getRiskyObservationDates() const;

    /** Return known cash flows */
    virtual CashFlowArraySP generateKnownCashFlows(
          const DateTime       today,
          const double         initialTrancheSize,
          const DateTimeArray  pastTrancheLossDates,
          const DoubleArray    pastTrancheLosses,
          const double         pastTrancheLoss,
          IForwardRatePricerSP model);

    /** Returns the pay date which is last in terms of time */
    virtual DateTime getLastPayDate() const;

    /** Returns the observation date which is last in terms of time */
    virtual DateTime getLastObservationDate() const;

    /** When to stop tweaking for Yield Curve type tweaks. Returns
        max(currentLastDate, this max date) where this max date =
        max(rate maturity, last pay date) */
    virtual DateTime lastYCSensDate(const DateTime& currentLastDate) const;

    /** Feedback method for getting information about a fee cashflow */
    virtual void getActiveFee(const DateTime&      withRespectTo,       // (I) get the fee whose accrual period contains this date
                              const DateTime&      earliestAccrualDate, // (I) for fee legs that dont specify accrue start dates, and we are interested in the first fee
                              IForwardRatePricerSP model,               // (I) for calculating the amount
                              DateTime&            accrueStartDate,     // (O) when the fee starts accruing
                              DateTime&            accrueEndDate,       // (O) when the fee finishes accruing
                              DateTime&            paymentDate,         // (O) when the fee is paid
                              double&              amount) const;       // (O) the cashflow amount

    /** Returns the leg notional */
    virtual double getFeeLegNotional() const;    

    /** Returns the leg notional */
    virtual void setFeeLegNotional(double newNotional);    

    //// Implementing tweakable with respect to CreditFeelegCoupon
    //// - a shift in coupon
    TweakOutcome sensShift(const PropertyTweak<CreditFeeLegCoupon>& shift);
    string sensName(const CreditFeeLegCoupon* shift) const;
    void sensRestore(const PropertyTweak<CreditFeeLegCoupon>& shift);


	/** creates the SV Generator  */
	virtual ICreditFeeLegSVGenSP createSVGen(
		const DateTime& valueDate,
		const DateTimeArray& modifiedFeeLegObservationDates,
		const DateTimeArray& productTimeline, //product timeline
		const IntArray& dateToDiscFactorIndex, //mapping from date to discFactorIndex
		double lossConfigNotional
		) const;

    //-------------------------------
    // IFixedRateCreditFeeLeg methods
    //-------------------------------
    /**Get the fee rate of the fee leg.*/
    virtual double getRate() const;
    
    /**Change the fee rate of the fee leg. Note that, if the rate is currently zero,
       this will throw an exception.*/
    virtual void setRate(double newRate);

private:
	// Methods
    static IObject* defaultConstructor();

    double calculateRate(DateTime             refixDate) const;

    void setFixingforThetaShift(const DateTime&      valueDate,
                                const YieldCurve*    discount,
                                const DateTime&      rollDate);

    bool isCouponNotionalTypeValid() const;

    ///// Fields

    bool                isFixed;  /** TRUE means that the whole leg is fixed */
    CDoubleArray        notionals;
    DateTimeArray       observationDates; // observation dates

    DateTimeArray       refixDates;
    DateTimeArray       accrualStartDates;  // accrual dates : can be non contiguous
    DateTimeArray       accrualEndDates;
    DateTimeArray       payDates;
    CDoubleArray        spreads; /* rate used is supplied fixings if isFixed is true
                                    otherwise it's weight * (historic fixing or rate
                                    from curve) + spread */
    CDoubleArray        weights;
    CDoubleArray        fixings;     // historic fixings or absolute levels if isFixed
    string              refixInterval;      //e.g. 3M or 6M....
    CouponNotionalTypes couponNotionalType;
    string              payDCC;
    string              rateDCC;
    string              badDayConvention;

    CashFlowArraySP     upfrontFees;  // Optional upfront payments


    // to get market data
    YieldCurveWrapper   couponCurve;
    DateTime            valueDate;

    // transient
    MaturityPeriodSP        refixIntervalSP;
    DayCountConventionSP    payDCCSP;
    DayCountConventionSP    rateDCCSP;
    BadDayConventionConstSP badDayConventionSP;
};

DRLIB_END_NAMESPACE

#endif // CREDIT_FEE_LEG_HPP
