//----------------------------------------------------------------------------
//
//   Group       : New York Credit QR
//
//   Author      : Charles Morcom
//
//   Filename    : IDiscountCurveRisky.hpp
//
//   Description : Interface for a credit curve. This is something that has
//                 reasonable notions of default probabilties, zero-coupon
//                 bond prices, credit protection, and annuity pricing
//                 with recovery of accrued-interest in case of default.
//
//   Date        : 10 January 2006
//
//----------------------------------------------------------------------------

#ifndef QR_IDISCOUNTCURVERISKY_HPP
#define QR_IDISCOUNTCURVERISKY_HPP

#include "edginc/IDiscountCurve.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/DefaultRates.hpp"

//#ifdef DEBUG
#define EXPOSE_IDISCOUNT_CURVE_RISKY_ADDINS
//#endif

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(IDiscountCurveRisky)


/**An abstract interface for a credit clean spread curve plus interest-rate discounting. 
 * It includes risk-free discounting,
 * default risk, and recovery (and hence loss) risk, and is sufficient for calculation of
 * forwards and simple bond/CDS-type instrument prices using simple deterministic discounting. 
 *
 * Methods named "pv" are inherited from IDiscountCurve. In this context, these refer to
 * PVs of risky cash-flows with no default recovery. This is natural for certain kinds of
 * "riskifications" of normally risk-free instruments; any instrument which has any kind
 * of default-contingent behaviour will probably need to use this interface instead.
 * */
class MARKET_DLL IDiscountCurveRisky : public virtual IDiscountCurve {
public:
    static CClassConstSP const TYPE;
    virtual ~IDiscountCurveRisky();

    /**Returns the probability that there are no default events from d1 to d2, 
       conditional that there are no default events up to d1*/
    virtual double survivalProb(
        const DateTime& d1, 
        const DateTime& d2) 
        const = 0;
    /**Probability of no default in [this->baseDate(),dt], given no default at base date.*/
    virtual double survivalProb(
        const DateTime& dt
        ) const = 0;

    /**Returns the market recovery rate on defaulted assets given default at 
       time defaultDate. This allows the possibility of a term-structure of recovery rates. */
    virtual double getRecovery(
        const DateTime& defaultDate
        ) const = 0;    

    /**Recovery-rate for immediate default.*/
    virtual double getRecovery() const = 0;
    
    /**Defines different types of recovery in case of default. In all cases 'R' means
     * the market recovery rate defined by this curve. If you want a specific recovery
     * rate, you should engineer it yourself using RECOVER_1. */
    enum RecoveryType {
        /**Recover the market recovery rate - like a bond's principal & accrued.*/
        RECOVER_R,
        /**Recover 100% of face at default - like a CDS's accrued interest.*/
        RECOVER_1,
        /**Recover (1-R) times face value - like a CDS's protection payment.*/
        RECOVER_1_MINUS_R,
        /**Recover nothing in case of default*/
        RECOVER_0
    };

    /**Returns the value at paymentDate (and conditional on no default before then) 
     * of a contingent claim that pays recTyp (1, R, 1-R or, trivially, 0)
     * in case of default between startDate and endDate, and zero otherwise.
     * The default payment is made with a delay of recoveryDelay calendar days.
     * Whether integration is done continuously or discretely and what, 
     * if any approximations,
     * are used is up to the curve implementation. */
    virtual double protectionPV(
        const DateTime&     paymentDate, 
        const DateTime&     startDt, 
        const DateTime&     endDt,
        RecoveryType        recTyp,
        double              recoveryDelay=0
        ) const = 0;

    /**Returns the value at paymentDate (and conditional on no default before then) 
     * of a contingent claim that pays recTyp (1, R, 1-R or, trivially, 0)
     * in case of default between startDate and endDate, and zero otherwise.
     * The default payment is made on recoveryDate, no matter when the default happens.
     * Whether integration is done continuously or discretely and what, 
     * if any approximations,
     * are used is up to the curve implementation. */
    virtual double protectionPV(
        const DateTime&     paymentDate, 
        const DateTime&     startDt, 
        const DateTime&     endDt,
        RecoveryType        recTyp,
        const DateTime&     recoveryDate
        ) const = 0;

    /**Returns the value at paymentDate (and conditional on no default before then) 
     * of a sequence of payments, with simple linear accrued-interest
     * recovery (so the claim on a coupon, C, payable at the end of accrual period (S,T) given
     * default at t is C(t-S)/(T_S)) in case of default defined by accruedRecType as for 
     * protectionPV. This allows the curve to compute default-accrual PVs itself
     * in a way which reflects the type of curve (e.g. flat forwards, etc.)
     * The accrual periods are the intervals between consecutive payment dates. To have a first
     * accrual period, you should set the first payment to zero, and then the first date is
     * the beginning of the first accrual period.
     * Payments before paymentDate are ignored, except to the extent that they affect accrued
     * interest due.
     * Any default payments are made with a delay of recoveryDelay calendar days after default.
     * Unless recType==RECOVER_0, the cash-flow dates MUST BE IN INCREASING ORDER, or the
     * default accrual calculations will not be very meaningful! if recType==RECOVER_0,
     * this method should return the same value as pv(payments, paymentDate), inherited
     * from IDiscountCurve. */
    virtual double annuityPV(
        const CashFlowArray&    payments,
        const DateTime&         paymentDate,
        RecoveryType            accruedRecTyp,
        double                  recoveryDelay=0,
	DateTime                accrueStartDate = DateTime()
        ) const = 0;

    /**Returns the value at paymentDate (and conditional on no default before then) 
     * of a sequence of payments, with simple linear accrued-interest
     * recovery (so the claim on a coupon, C, payable at the end of accrual period (S,T) given
     * default at t is C(t-S)/(T_S)) in case of default defined by accruedRecType as for 
     * protectionPV. This allows the curve to compute default-accrual PVs itself
     * in a way which reflects the type of curve (e.g. flat forwards, etc.)
     * The accrual periods are the intervals between consecutive payment dates. To have a first
     * accrual period, you should set the first payment to zero, and then the first date is
     * the beginning of the first accrual period.
     * Payments before paymentDate are ignored, except to the extent that they affect accrued
     * interest due.
     * Any default payments are made on recoveryDate, no matter when the default happens.
     * Unless recType==RECOVER_0, the cash-flow dates MUST BE IN INCREASING ORDER, or the
     * default accrual calculations will not be very meaningful! if recType==RECOVER_0,
     * this method should return the same value as pv(payments, paymentDate), inherited
     * from IDiscountCurve. */
    virtual double annuityPV(
        const CashFlowArray&    payments,
        const DateTime&         paymentDate,
        RecoveryType            accruedRecTyp,
        const DateTime&         recoveryDate,
	DateTime                accrueStartDate = DateTime()
	) const = 0;         

    //------------------------------------------------------
    // Additional Riskless PV methods for flexibility
    //------------------------------------------------------

    /** Compute discount factor between two dates. This is the price
     * at time date1 of a zero-coupon (discount) bond that pays 1
     * at time date2 if it still exists then. 
     * Note that, if this bond is risky, this method
     * in IDiscountCurveRisky means the PV of such a bond which knocks
     * out on default with no recovery at all, and the price is 
     * contingent on no default before date1.
     * @param date1 payment settlement date (conditional on existence at date1)
     * @param date2 payment/zero-coupon bond maturity date
     * @return Discount factor between date1 & date2
     */
    virtual double risklessPV(const DateTime& date1, 
                              const DateTime& date2) const = 0;
    
    /** Compute price for settlement today of a zero-coupon bond
     * maturing at date. Note settlement is to TODAY, not to
     * the curve spot date. This is because some curves
     * may have ambiguous spot-dates - for example should a combined
     * credit and rates curve have spot date T+1 or T+2?
     * @param date To get discount factor/PV for
     * @return Discount factor between today & given date
     */
    virtual double risklessPV(const DateTime& date) const = 0;
    
    /** Calculates present value to baseDate of supplied cash flows conditional
        on continued existence (i.e. no default for a risky curve)
        Cash-flows on or before baseDate are ignored. No
        ordering of the cashflows is assumed. */
    virtual double risklessPV(const CashFlowArray& cashFlows,
                              const DateTime&      baseDate) const = 0;

    /** accessor methods for logOfDiscFactorKey */
    virtual IDiscountCurve::IKey* getDiscountKey() const = 0;
    virtual DefaultRates::IKey* getRiskyKey() const = 0;

    /** Returns a DefaultRates object which gives access to 
        useful functionality including "default rates", aka clean default
        spreads. The information needed to bootstrap the clean spreads is
        obtained from this object.
        The DefaultRate returned is immutable and therefore it will not be
        changed - which means that there is no need to clone it */
    virtual DefaultRatesSP getDefaultRates() const = 0;

private:
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE
#endif
