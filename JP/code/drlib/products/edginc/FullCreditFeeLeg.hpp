//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : FullCreditFeeLeg.hpp
//
//   Description : Implementation of ICreditFeeLeg as a list of cashflows
//
//   Author      : Antoine Gregoire
//
//   Date        : 3 Nov 2004
//
//----------------------------------------------------------------------------

#ifndef FULL_CREDIT_FEE_LEG_HPP
#define FULL_CREDIT_FEE_LEG_HPP

#include "edginc/CreditFeeLegWithPV.hpp"
#include "edginc/IFixedRateCreditFeeLeg.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implementation of ICreditFeeLeg using ICreditCashFlow */
class PRODUCTS_DLL FullCreditFeeLeg :
    public CreditFeeLegWithPV,
    virtual public IFixedRateCreditFeeLeg //temporary
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    /** Explicit constructor */
    FullCreditFeeLeg(AbstractCashFlowArraySP creditCashFlows);

    /** Destructor */
    virtual ~FullCreditFeeLeg();

    /** GetMarket implementation*/
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

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

    /** creates the SV Generator  */
	virtual ICreditFeeLegSVGenSP createSVGen(
		const DateTime& valDate,
		const DateTimeArray& modifiedFeeLegObservationDates,
		const DateTimeArray& productTimeline, //product timeline
		const IntArray& dateToDiscFactorIndex, //map from a date to an index on the DiscountFactor SV
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
    /** Constructor (for reflection) */
    FullCreditFeeLeg();

    static IObject* defaultConstructor();

    ///// Fields
    AbstractCashFlowArraySP creditCashFlows;
};

DRLIB_END_NAMESPACE

#endif // FULL_CREDIT_FEE_LEG_HPP
