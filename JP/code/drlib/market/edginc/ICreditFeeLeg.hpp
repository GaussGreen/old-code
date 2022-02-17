//----------------------------------------------------------------------------
//
//   Filename    : ICreditFeeLeg.hpp
//
//   Description : Interface for credit fee leg
//
//----------------------------------------------------------------------------

#ifndef ICREDIT_FEE_LEG_HPP
#define ICREDIT_FEE_LEG_HPP

#include "edginc/Atomic.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/AbstractCashFlow.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
#include "edginc/ICreditLossConfig.hpp"
#include "edginc/StateVariable.hpp"
#include "edginc/ElemStateVariableGen.hpp"
#include "edginc/CreditCashFlow.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(IBadDayAdjuster);
FORWARD_DECLARE(IForwardRatePricer);
FORWARD_DECLARE(IDecretionCurve);
FORWARD_DECLARE(ICreditContingentLeg);
FORWARD_DECLARE(ICreditEventOverrideName);

/** Interface for credit fee leg */
class MARKET_DLL ICreditFeeLeg : public virtual IGetMarket {
public:
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
        IForwardRatePricerSP        model) const = 0;

    /** Return all cash flow dates */
    virtual DateTimeArraySP getCashFlowDates() const = 0;

    /** Return risky cash flow dates */
    virtual DateTimeArraySP getRiskyCashFlowDates() const = 0;

    /** Return riskfree cash flow dates */
    virtual DateTimeArraySP getRisklessCashFlowDates() const = 0;

    /** Returns the accrual periods (start and end dates) in this fee leg. */
    virtual AccrualPeriodArrayConstSP getAccrualPeriods() const = 0;

    /** Returns risky accrual periods (start and end dates) in this fee leg. */
    virtual AccrualPeriodArrayConstSP getRiskyAccrualPeriods() const = 0;

    /** Returns risky accrual periods (start and end dates) in this fee leg. */
    virtual AccrualPeriodArrayConstSP getRisklessAccrualPeriods() const = 0;

    /** Returns all cash flows */
    virtual AbstractCashFlowArrayConstSP getCashFlows(IForwardRatePricerSP model) const = 0;

    /** Returns risky cash flows only */
    virtual CashFlowArraySP getRiskyCashFlows(IForwardRatePricerSP model) const = 0;

    /** Returns risk free cash flows only */
    virtual CashFlowArraySP getRisklessCashFlows(IForwardRatePricerSP model) const = 0;

    /** Returns risky notional dates */
    virtual DateTimeArraySP getRiskyNotionalDates(IForwardRatePricerSP model) const = 0;

	/** Returns risky coupon notional types */
	virtual CouponNotionalTypesArraySP getRiskyCouponNotionalTypes() const = 0;

	/** Returns risky observation dates */
	virtual DateTimeArraySP getRiskyObservationDates() const = 0;

    /** Return known cash flows corresponding to a CDO tranche */
    virtual CashFlowArraySP generateKnownCashFlows(
         const DateTime       today,
         const double         initialTrancheSize,
         const DateTimeArray  pastTrancheLossDates,
         const DoubleArray    pastTrancheLosses,
         const double         pastTrancheLoss,
         IForwardRatePricerSP model) = 0;

    /** Return known cash flows corresponding to a non-defaulted CDS */
    virtual CashFlowArraySP generateKnownCashFlows(IForwardRatePricerSP model) const = 0;

    /** Estimates the known cash flows (so they are not really "known")
     * corresponding to a CDO tranche - takes into account estimated losses
     * in the future */
    virtual CashFlowArraySP estimateKnownCashFlows(
         const double         initialTrancheSize,
         const DateTimeArray  pastTrancheLossDates,
         const DoubleArray    pastTrancheLosses,
         IForwardRatePricerSP model) = 0;

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
        const DateTime&            lastTriggerDate) const = 0;

    /** Returns the pay date which is last in terms of time */
    virtual DateTime getLastPayDate() const = 0;

    /** Returns the observation date which is last in terms of time */
    virtual DateTime getLastObservationDate() const = 0;

    /** When to stop tweaking for Yield Curve type tweaks */
    virtual DateTime lastYCSensDate(const DateTime& currentLastDate) const = 0;

    /** Feedback method for getting information about a fee cashflow */
    virtual void getActiveFee(const DateTime&      withRespectTo,       // (I) get the fee whose accrual period contains this date
                              const DateTime&      earliestAccrualDate, // (I) for fee legs that dont specify accrue start dates, and we are interested in the first fee
                              IForwardRatePricerSP model,               // (I) for calculating the amount
                              DateTime&            accrueStartDate,     // (O) when the fee starts accruing
                              DateTime&            accrueEndDate,       // (O) when the fee finishes accruing
                              DateTime&            paymentDate,         // (O) when the fee is paid
                              double&              amount) const = 0;   // (O) the cashflow amount

    /** Compute PV of fee leg at valuation date.
     *  Optionally, include accrued interest in the result */
    virtual double getFeeLegPV(const DateTime&              valuationDate, 
                               const DateTime&              earliestRiskyDate,
                               const DateTime&              latestRiskyDate,
                               const IDiscountCurve&        discount,
                               const IDiscountCurveRisky&   crv,
                               const IDecretionCurveConstSP prepay,
                               const bool                   includeAccrued,
                               const DayCountConventionSP   dcc,
                               IForwardRatePricerSP         model) const = 0;

    /**
    * Compute PV of fee leg at valuationDate, but for unconditional payment at
    * paymentDate, conditional on no new defaults before valuationDate.
    * Optionally, include accrued interest in the result
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
                               IForwardRatePricerSP         model) const = 0;

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
                               IForwardRatePricerSP         model) const = 0;

    /** Returns the accrued interest */
    virtual double getFeeLegAI(const DateTime&              valuationDate, 
                               const DateTime&              paymentDate, 
                               const DateTime&              earliestAccrualStart,
                               const DateTime&              latestAccrualEnd,
                               const DayCountConventionSP   dcc, //allows an override to be specified
                               const IDiscountCurveRisky&   crv,
                               const IDiscountCurveConstSP  discount,
                               const IDecretionCurveConstSP prepay,
                               IForwardRatePricerSP         model) const = 0;

    /** Prices the leg sufferring a default, under the (optional) credit event */
    virtual double getFeeLegDefaultedPV(const DateTime&              valuationDate,
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
                                        DateTime                     lastTriggerDate) const = 0;

    /* Computes the interest accrued up to valuationDate */
    virtual double getAccruedInterest(const DateTime& valuationDate,
                                      IForwardRatePricerSP model) const = 0;

    /** Returns the leg notional */
    virtual double getFeeLegNotional() const = 0;    

    /** Returns the leg notional */
    virtual void setFeeLegNotional(double newNotional) = 0;    

    /** Inner class: declaration of State Variable Gen */
	FORWARD_DECLARE(ISVGen);

	/** creates the SV Generator  */
	virtual smartPtr<ICreditFeeLeg::ISVGen> createSVGen(
		const DateTime& valueDate,
		const DateTimeArray& modifiedFeeLegObservationDates,
		const DateTimeArray& productTimeline, //product timeline
		const IntArray& dateToDiscFactorIndex, //map from a date to an index on the DiscountFactor SV
		double lossConfigNotional
		) const = 0;

/** Static utility method to be used for converting productTimeline to
	a map that facilitates quick conversion from any Date to a productTimeline index.
	if Date does not lie on the productTimeline, then it returns -1
*/
static void
	buildDateToProductTimelineIndexMap(
		IntArray& dateToProductTimelineIndexMap,
		const DateTimeArray& productTimeline); //input

/** Static utility method for creating scaling factors for the risky cashflows
	of the ICreditFeeLeg.
	scalingFactor = Coupon X YearFraction /lossConfigNotional * notional

*/

static void
buildScalingFactors(
	DoubleArray& scalingFactor, //output
	DoubleArray& historicalScalingFactor, //output
	const DateTime& valueDate, //input
	const AccrualPeriodArray& accrualPeriods, //input
	const DateTimeArray& payDates, //input
	const CouponNotionalTypesArray& couponNotionalTypes,	//input
	const DayCountConventionArray& dayCountConventions, //input
	const DoubleArray& coupons, //input
	const DoubleArray& notionals, //input
	double lossConfigNotional); //input
};

typedef smartPtr<ICreditFeeLeg> ICreditFeeLegSP;
typedef smartConstPtr<ICreditFeeLeg> ICreditFeeLegConstSP;

/** Describes something (like a CDS or bond) which has a fee leg which is
  * publicly accessible.
  */
class MARKET_DLL IHasCreditFeeLeg : public virtual IObject {
public:

    /** Return the fee leg */
    virtual ICreditFeeLegSP getFeeLeg() const = 0;

    static CClassConstSP const TYPE;
    virtual ~IHasCreditFeeLeg();
private:
    static void load(CClassSP& clazz);
};

// ############################################################################

/** A fee leg generator is an object which is capable of generating a fee
  * leg, if you give it new start and end dates and a fee rate. */
class MARKET_DLL ICreditFeeLegGenerator : virtual public IObject {
public:
    /**Generate a new fee leg with specified start and end dates and fee rate whose other
       characteristics are defined by this.*/
    virtual ICreditFeeLegSP generateCreditFeeLeg(
        const DateTime& startDate,
        const DateTime& endDate,
        double feeRate
        ) const = 0;

    static CClassConstSP const TYPE;
    virtual ~ICreditFeeLegGenerator();
private:
    static void load(CClassSP& clazz);
};

/** Represents the SV generator for the ICreditFeeLeg  */
class MARKET_DLL ICreditFeeLeg::ISVGen :	public virtual IElemStateVariableGen,
											public virtual VirtualDestructorBase

{
public:
	/** forward declaration of the inner class representing the actual state variable   */
	FORWARD_DECLARE(ISV);

	/** constructor */
	ISVGen() {}

	/** virtual destructor */
	virtual ~ISVGen() {}

	/** Same as create but avoids dynamic cast */
	virtual smartPtr<ICreditFeeLeg::ISVGen::ISV> createNewSV(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const = 0;

// ##  Methods of the parent interface, IStateVariable
	/** Fetches the state variable from the stateGenerator */
	virtual IStateVariableSP create(
		IStateVariableSP oldStateVar,
		IStateVariableGen::IStateGen* stateGen) const = 0; //the second argument above is PathGenerator

// ##  Methods of the parent interface, IElemStateVariableGen - revisit
	virtual void attachSVGen(IElemStateVariableGenVisitor*) const {}

private:
    ISVGen(const ISVGen&);
    ISVGen& operator=(const ISVGen&);

}; // ICreditFeeLeg::SVGen

typedef smartPtr<ICreditFeeLeg::ISVGen> ICreditFeeLegSVGenSP;
typedef smartConstPtr<ICreditFeeLeg::ISVGen> ICreditFeeLegSVGenConstSP;

// ############################################################################

/** Interface for state variables representing the ICreditFeeLeg
	the credit event time into a bucket
*/
FORWARD_DECLARE(SVDiscFactor)

class MARKET_DLL ICreditFeeLeg::ISVGen::ISV : public virtual IStateVariable
{
public:
	/** Constructor  */
	ISV() {};

	/** Virtual Destructor  */
	virtual ~ISV() {};

	/** computes the pv - given the result of the ICreditLossConfig::IIndexedSVGen::ISV
		later, we may need another method that computes the pv for ICreditLossConfig::ISVGen::ISV
		DiscountFactor SV would be a member of ISVGen and perhaps passed through the constructor
	*/
	virtual double pv(
		double& riskFreeCashFlow,
		double& unitCouponPV,
		bool computeUnitCouponPV,
		const DateTime& valueDate,
		const CreditLossConfigIndexedSVResultPath& svResult,
		const SVDiscFactor& discountCurve
		) const = 0;

	/** method of the parent interface, IStateVariable */
	virtual bool doingPast() const = 0;


private:
    ISV(const ISV&);
    ISV& operator=(const ISV&);
}; // ICreditFeeLeg::ISVGen::ISV

typedef smartPtr<ICreditFeeLeg::ISVGen::ISV> ICreditFeeLegSVSP;
typedef smartConstPtr<ICreditFeeLeg::ISVGen::ISV> ICreditFeeLegVConstSP;

// ############################################################################

DRLIB_END_NAMESPACE

#endif // ICREDIT_FEE_LEG_HPP
