//----------------------------------------
//
//   Group       : Credit Hybrids QR&D
//
//   Filename    : CreditCashFlow.hpp
//
//   Description : Credit style cashflows
//
//   Author      : Gordon Stephens
//
//   Date        : 8 June 2005
//----------------------------------------

#ifndef CREDITCASHFLOW_HPP
#define CREDITCASHFLOW_HPP

#include "edginc/AbstractCashFlow.hpp"
#include "edginc/RiskFreeCashFlow.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
#include "edginc/CCMPriceUtil.hpp"

DRLIB_BEGIN_NAMESPACE

enum CouponNotionalTypes
{
    OBSERVATION_DATE = 0, // fixed observation date (observation date must be separately specified)
    AVERAGE = 1           // average between accrued start and accrued end date
};

/** This goes in the header file. It is only needed if you want an array of
    your enums. In particular, this follows the usual paradigm of defining
    arrays in that this template specialisation is only need if you don't want
    an array of smart pointers. */
template <> class MARKET_DLL arrayObjectCast<CouponNotionalTypes>
{
public:
    /** Casts array element to an IObject */
    static IObjectSP toIObject(CouponNotionalTypes value);

    /** Sets the value of the indexed component of the specified
        array object to the specified new value. */
    static CouponNotionalTypes fromIObject(IObjectSP value);
};

//required for arrays (implicitly of SP's)
//typedef smartPtr<BoxedEnum<CouponNotionalTypes> > CouponNotionalTypesSP;
//typedef array<CouponNotionalTypesSP> CouponNotionalTypesArray;
typedef array<CouponNotionalTypes> CouponNotionalTypesArray;
typedef smartPtr<CouponNotionalTypesArray> CouponNotionalTypesArraySP;

class MARKET_DLL CreditCashFlow : public AbstractCashFlow
{
public:
    static CClassConstSP const TYPE;

    //CreditCashFlow methods
    DateTime getObservationDate() const;

    //// When to stop tweaking for Yield Curve type tweaks
    virtual DateTime lastYCSensDate() const;

    //// Return known cash flow corresponding to a CDO tranche
    CashFlow getKnownCashFlow(const double         initialTrancheSize,
                              const DateTimeArray  pastTrancheLossDates,
                              const DoubleArray    pastTrancheLosses,
                              IForwardRatePricerSP model);

    //// Returns known cash flow corresponding to a defaulted CDS, computing
    //// accrued as required.
    //// An empty defaultDeterminationDate means the name has not been 
    //// triggered, so the whole cashflow will be accrued.
    virtual CashFlow getKnownCashFlowGivenDefault(
        const DateTime&      defaultDeterminationDate,
        const DateTime&      accrualPaymentDate,
        IForwardRatePricerSP model) const;

    /** Returns the fraction of this cashflow accrued in case of default.
        Adhoc cashflows do not accrue (since they have no accrue period) and
        an empty defaultDeterminationDate is interpreted as meaning that the 
        name did not default, so the complete cashflow is returned. */
    virtual CashFlow getAccruedCashFlow(
        const DateTime&      valuationDate,
        IForwardRatePricerSP model) const;

    //// Returns the outstanding notional on a tranche
    //// based off cashflow observation date(s)
    //// against which this cashflow may be paid
    double expectedTrancheOutstanding(
        const double                            initialNotional,
        const double                            outstandingNotional,
        const DateTime&                         today,
        const CashFlowArray&                    pastTrancheLosses,
        const IDiscountCurveRiskySP             effectiveCurve,
        const IDiscountCurveRisky::RecoveryType recTyp,
        const double                            recoveryDelay) const;

    //AbstractCashFlow methods
    virtual DateTime getPayDate()  const;
    virtual double   getAmount(IForwardRatePricerSP model) const;
    virtual double   getNotional() const;
    virtual void     setNotional(double newNotional);
    virtual double   getRate() const; //implicitly a fixed rate
    virtual void     setRate(double newRate); //implicitly a fixed rate
    
    /** Returns cash flow amount using given rate override (if relevant) */
    virtual double getAmount(double rateOverride) const;

    //CObject methods
    virtual void validatePop2Object();

    /** GetMarket implementation*/
    void getMarket(const IModel* model, const MarketData* market);
    
    /** Says if this cash flow is actually an adhoc cash flow */
    virtual bool isAdhoc() const;

    /** Says if this cash flow is risk free */
    virtual bool isRiskFree() const;

    /** Returns the accrual period (start and end dates) corresponding to 
     *  this cash flow. */
    virtual AccrualPeriodSP getAccrualPeriod() const;

	/** Access method for the coupon notional type
     */
	virtual const CouponNotionalTypes& getCouponNotionalType() const; 

	/** Access method for the riskFree cash flow
     */
    virtual RiskFreeCashFlowConstSP getRiskFreeCashFlow() const; 
	
	CreditCashFlow(RiskFreeCashFlowSP  rfCashFlow,
                   CouponNotionalTypes couponNotionalType,
                   DateTimeSP          observationDate);

    CreditCashFlow();

    ~CreditCashFlow();

private:
    //infrastructure support
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    /** Checks if couponNotionalType is valid */
    bool isCouponNotionalTypeValid();

    /** Computes the average notional in this accrue period, taking
        the given losses into account */
    double computeAverageNotional(const DateTimeArray& pastLossDates,
                                  const DoubleArray&   pastLosses,
                                  const double         initialNotional) const;

    /** Returns the (fraction of) this cashflow accrued in case of default.
        Adhoc cashflows do not accrue (since they have no accrue period) and
        a defaultDeterminationDate is interpreted as meaning that the name did 
        not default, so the complete cashflow is returned.
        CAUTION: Risky fee amounts are obtained by calling getAmount();
        If the cashflow is not an adhoc cashflow, the amounts are internally 
        computed by multiplying the coupon by the notional and therefore 
        both will typically have the same sign (negative coupons are allowed
        though). For CDSs this is wrong (notional > 0 means long protection, 
        so fees should typically be negative). 
        Therefore non-adhoc cashflow amounts will be multiplied 
        by -1 here, even though this is not consistent with other methods 
        in this class! */
    CashFlow getAccruedCashFlow(
        const DateTime&      defaultDeterminationDate,
        const DateTime&      accrualPaymentDate,
        const DateTime&      excludeCashFlowsBeforeThisDate,
        IForwardRatePricerSP model) const;


    // FIELDS
    RiskFreeCashFlowSP  rfCashFlow;         //the underlying cashflow description
    CouponNotionalTypes couponNotionalType; //Coupon notional type, ie when do we observe
                                            //a credit risky quantity (eg: a tranche notional)
    DateTimeSP          observationDate;    //required for OBSERVATION_DATE flows
};

typedef smartPtr<CreditCashFlow> CreditCashFlowSP;
typedef smartConstPtr<CreditCashFlow> CreditCashFlowConstSP;
typedef array<CreditCashFlowSP, CreditCashFlow> CreditCashFlowArray;
typedef smartPtr<CreditCashFlowArray> CreditCashFlowArraySP;

DRLIB_END_NAMESPACE

#endif
