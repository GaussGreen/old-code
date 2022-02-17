//----------------------------------------------------------------------------
//
// Group       : Credit Hybrids QR
//
// Description : Risky discount curve to price FtD contingent legs
//
// Date        : November 2006
//
//----------------------------------------------------------------------------

#ifndef QR_FTDEFFECTIVECURVE_HPP
#define QR_FTDEFFECTIVECURVE_HPP

#include "edginc/DefaultRates.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/CCMPriceUtil.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/ICondLossDistributionsGen.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"
#include "edginc/NToDefaultLossConfig.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ConvolutionProduct);
//using SP's within the class methods
FORWARD_DECLARE(FtDEffectiveCurve);


/** Class to compute the pv of the contingent leg of an FtD.
    It requires "ICondLossDistributionsGenKeys" because it is not possible 
    to pre-compute the data (survival probs, etc) required to be used here.
    So, instead, the data required to compute those is passed in, and they
    are computed on the fly.
    The issue is that the recoveryDelay (which is only know at the point
    of pricing the contingent leg) is required in the integration accross
    market factors... so we need to postpone computing everything until the
    last minute, when everything is known.
    CAUTION: Many methods are not implemented and will throw an exception
    if invoked */
class CONVOLUTION_DLL FtDEffectiveCurve : public CObject,
                                          public virtual IDiscountCurveRisky
{
public:
    static CClassConstSP const TYPE;
    virtual ~FtDEffectiveCurve();

    ////Constructor with a risk free discount curve
    FtDEffectiveCurve(
        const DateTime&             valueDate,
        const YieldCurveConstSP     dsc,
        DateTimeArraySP             timeline,
        const vector<ICondLossDistributionsGenKeyArraySP>& keysByDate,
        IConditionalDefaultsModelSP condDefaultsModel,
        DoubleArraySP               namesLoss,
        const double                ftdNotional);

    //----------------------------
	// IDiscountCurveRisky methods
	//----------------------------
    /**Returns the probability that there are no default events from d1 to d2, 
       conditional that there are no default events up to d1*/
    virtual double survivalProb(const DateTime& d1, const DateTime& d2) const;

    /**Probability of no default in [this->baseDate(),dt], given no default at base date.*/
    virtual double survivalProb(const DateTime& dt) const;

    /**Returns the market recovery rate on defaulted assets given default at 
       time defaultDate. This allows the possibility of a term-structure of recovery rates. */
    virtual double getRecovery(const DateTime& defaultDate) const;    

    /**Recovery-rate for immediate default.*/
    virtual double getRecovery() const;

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
        ) const;

    /**Returns the value at paymentDate (and conditional on no default before then) 
     * of a contingent claim that pays recTyp (1, R, 1-R or, trivially, 0)
     * in case of default between startDate and endDate, and zero otherwise.
     * The default payment is made on recoveryDate, no matter when the default happens.
     * Whether integration is done continuously or discretely and what, 
     * if any approximations, * are used is up to the curve implementation. */
    virtual double protectionPV(const DateTime& paymentDate, 
                                const DateTime& startDt, 
                                const DateTime& endDt,
                                RecoveryType    recTyp,
                                const DateTime& recoveryDate) const;

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
    virtual double annuityPV(const CashFlowArray& payments,
                             const DateTime&      paymentDate,
                             RecoveryType         accruedRecTyp,
                             double               recoveryDelay=0,
                             DateTime             accrueStartDate = DateTime()) const;

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
        const CashFlowArray& payments,
        const DateTime&      paymentDate,
        RecoveryType         accruedRecTyp,
        const DateTime&      recoveryDate,
        DateTime             accrueStartDate = DateTime()
	) const;

    //------------------------------------------------------
    // Riskless PV methods
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
     * @return Discount factor between date1 & date2 */
    virtual double risklessPV(const DateTime& date1, 
                              const DateTime& date2) const;
    
    /** Compute price for settlement today of a zero-coupon bond
     * maturing at date. Note settlement is to TODAY, not to
     * the curve spot date. This is because some curves
     * may have ambiguous spot-dates - for example should a combined
     * credit and rates curve have spot date T+1 or T+2?
     * @param date To get discount factor/PV for
     * @return Discount factor between today & given date */
    virtual double risklessPV(const DateTime& date) const;
    
    /** Calculates present value to baseDate of supplied cash flows conditional
        on continued existence (i.e. no default for a risky curve)
        Cash-flows on or before baseDate are ignored. No
        ordering of the cashflows is assumed. */
    virtual double risklessPV(const CashFlowArray& cashFlows,
                              const DateTime&      baseDate) const;
    //-----------------------
    // IDiscountCurve methods
    //-----------------------

    /** @return Discounting curve's currency - typically this is the
        ISO code eg "USD" (although it is up to the client - you may
        assume however that yield curves in the same currency return the
        same value in getCcy(). Note it is NOT the name of the yield
        curve (which might be eg "GBP-LIBOR"). */
    virtual string getCcy() const;

    /** Compute discount factor between two dates. This is the price
     * at time date1 of a zero-coupon (discount) bond that pays 1
     * at time date2 if it still exists then. 
     * Note that, if this bond is risky, this method
     * in IDiscountCurveRisky means the PV of such a bond which knocks
     * out on default with no recovery at all, and the price is 
     * contingent on no default before date1.
     * @param date1 payment settlement date (conditional on existence at date1)
     * @param date2 payment/zero-coupon bond maturity date
     * @return Discount factor between date1 & date2 */
    virtual double pv(const DateTime& date1, 
                      const DateTime& date2) const;
    
    /** Compute price for settlement today of a zero-coupon bond
     * maturing at date. Note settlement is to TODAY, not to
     * the curve spot date. This is because some curves
     * may have ambiguous spot-dates - for example should a combined
     * credit and rates curve have spot date T+1 or T+2?
     * @param date To get discount factor/PV for
     * @return Discount factor between today & given date */
    virtual double pv(const DateTime& date) const;
    
    /** Calculates present value to baseDate of supplied cash flows conditional
        on continued existence (i.e. no default for a risky curve)
        Cash-flows on or before baseDate are ignored. No
        ordering of the cashflows is assumed. */
    virtual double pv(const CashFlowArray& cashFlows,
                      const DateTime&      baseDate) const;

    /** return the bootstrapped dates */
    virtual DateTimeArray zeroDates() const;

    /** Returns a key used to optimise repeated calculations of discount
        factors/forward rate. The calc method for this key returns the 
        natural logarithm of the discount factor (or equivalently the
        product of the forward rate (continuous, Act/365F) and the negative
        year fraction (Act/365F) betweeen the two dates.
        The default implementation has no performance improvements. */
    virtual IDiscountCurve::IKey* logOfDiscFactorKey() const;

    // accessor methods for logOfDiscFactorKey
    virtual IDiscountCurve::IKey* getDiscountKey() const;
    virtual DefaultRates::IKey* getRiskyKey() const;

    /** Returns a DefaultRates object which gives access to 
        useful functionality including "default rates", aka clean default
        spreads. The information needed to bootstrap the clean spreads is
        obtained from this object.
        The DefaultRate returned is immutable and therefore it will not be
        changed - which means that there is no need to clone it */
    virtual DefaultRatesSP getDefaultRates() const;

    //-----------------------------------------------
    // IGetMarket methods - exposed by IDiscountCurve
    //-----------------------------------------------

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    IObject* clone() const;

private:
    FtDEffectiveCurve(const FtDEffectiveCurve& rhs);
    FtDEffectiveCurve& operator=(const FtDEffectiveCurve& rhs);

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    FtDEffectiveCurve();

    ////Post construction basic validation
    void validate();

    //-------
    // Fields
    //-------
    DateTime                    valueDate;
    // 	string                      interpolationStyle; // Currently ignored (!)
    DateTimeArraySP             timeline;           // input timeline
    const vector<ICondLossDistributionsGenKeyArraySP> keysByDate;
    IConditionalDefaultsModelSP condDefaultsModel;
    IDiscountCurveSP            discount;           //may be risk-free or risky $unregistered
    DoubleArraySP               namesLoss;
    const double                ftdNotional;
};


DRLIB_END_NAMESPACE
#endif
