//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : Bond.hpp
//
//   Description : Virtual base class for defining bond cashflows and accrued
//                 interest
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : September 28, 2001
//
//
//----------------------------------------------------------------------------

#ifndef BOND_HPP
#define BOND_HPP

#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/AccrualCalendar.hpp"
#include "edginc/Settlement.hpp"

DRLIB_BEGIN_NAMESPACE
/** A Bond defines the cashflows and accrued interest for a bond. */

class PRODUCTS_DLL Bond:public CObject {
public:
    static CClassConstSP const TYPE;
    friend class BondHelper;

    virtual ~Bond();

    /** change the settlement of the bond if applicable */
    virtual void setSettlement(SettlementSP newSettle) =  0;

    /** get the accrued interest at a date */
    virtual double getAccruedAtDate(const DateTime &date) const = 0;

    /** returns only the cashFlows that occur after the start date */
    virtual CashFlowArraySP getCashFlows(const DateTime &startDate) const = 0;

    /** returns all the cashFlows */
    virtual CashFlowArraySP getCashFlows() const = 0;

    virtual CashFlowArrayConstSP getCashFlowRef() const;

    /** returns all coupons */
    virtual CashFlowArraySP getCoupons() const;

    /** returns only the coupons that occur after the start date */
    virtual CashFlowArraySP getCoupons(const DateTime &startDate) const;

    virtual DateTimeArraySP getExCouponDates() const;

    virtual DateTime getMaturityDate() const = 0;

    virtual DateTime getUnadjMaturityDate() const;

    virtual double getFaceValue() const = 0;

    virtual double getNotional(const DateTime &date) const;

    virtual DoubleArraySP getNotionalsAtDates(DateTimeArrayConstSP dates) const;

    virtual CashFlowArraySP getRedemptionPayments() const;

    virtual double getRedemption() const = 0;

    virtual double getRedemption(bool asPercent) const;

    virtual double getMaturityCoupon() const = 0;

    virtual void   getMarket(const IModel* model, const MarketData* market) {}

    /** value future bond cashFlows at date. Dirty value */
    virtual double presentValue(const DateTime &valueDate, YieldCurveConstSP discount) const;

    virtual double redemptionPV(const DateTime &baseDate, YieldCurveConstSP discount) const;

    virtual double duration(const DateTime &valueDate, YieldCurveConstSP discount) const;

    virtual double convexity(const DateTime &valueDate, YieldCurveConstSP discount) const;

    /** value of all coupons paid between fromDate and to Date on toDate */
    virtual double couponsPV(const DateTime&   fromDate,
                     const DateTime&   toDate,
                     YieldCurveConstSP discount);

    /** return the cashFlows where the pv of the coupons is paid on the ex date  */
    virtual CashFlowArraySP getExAdjCashFlows(const DateTime &startDate, YieldCurveConstSP discount) const;

    /** required for Yield to Maturity calculation */
    virtual double priceFromYield(double yield, bool clean, DateTime yDate) const = 0;

    /** Yield to Maturity calculation */
    double yieldToMaturity(double yPrice, bool clean, DateTime yDate);

    // PVBP calculation (based on YTM-0.5bp, YTM+0.5bp)
    double pvbp(double mktPrice, bool priceIsClean, DateTime fromDate);

    /** get the equivalent bond which matures on the put date. Used for yield to put. */
    virtual Bond* getBondToPut(const DateTime &putDate, double putLevel) const = 0;

    /** Yield to first put calculation */
    double yieldToFirstPut(double yPrice, bool clean, DateTime putDate, double putLevel, DateTime yDate);

    /** get coupon details at base date */
    void getCouponsDetails(const DateTime& baseDate, DateTime& startDate, DateTime& endDate,
                           double& couponRate, double& accruedInterest) const;

    virtual DateTime getAccrualStartDate() const = 0;

    virtual double   getBondCouponRate() const;

    // get all cashflows known today
    virtual CashFlowArraySP getKnownCashFlows() const;

    virtual DateTimeArraySP getPaymentDates() const;

    virtual CashFlowArraySP getCouponRates();

    virtual CashFlowArraySP getDayCountFractions();

    // when does bond settle?
    virtual DateTime settles (const DateTime& tradeDate) const = 0;

    // return any cashflow due in the next CPND_DATE_WINDOW calendar days
    // together with the coupon due date and the due date adjusted to next valid business day
    virtual AccrualCalendarArraySP couponDue(const DateTime& fromDate) const;

    // return ACCA_DATE_WINDOW calendar days of accrued interest data starting with fromDate
    virtual AccrualCalendarArraySP accrualCalendar(const DateTime& date) const;

protected:
    Bond(CClassConstSP clazz);
    static const double ONE_BASIS_POINT;

};

typedef smartConstPtr<Bond> BondConstSP;
typedef smartPtr<Bond> BondSP;
#ifndef QLIB_BOND_CPP
EXTERN_TEMPLATE(class PRODUCTS_DLL_SP smartPtr<Bond>);
#else
INSTANTIATE_TEMPLATE(class PRODUCTS_DLL smartPtr<Bond>);
#endif



DRLIB_END_NAMESPACE
#endif
