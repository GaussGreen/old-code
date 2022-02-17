//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ExplicitBond.hpp
//
//   Description : Explicit Bond implementation 
//
//   Author      : André Segger
//
//   Date        : 10 September 2003
//
//
//----------------------------------------------------------------------------

#ifndef EXPLICIT_BOND_HPP
#define EXPLICIT_BOND_HPP

#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Asset.hpp"
#include "edginc/Bond.hpp"
#include "edginc/Settlement.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/RateParallel.hpp"
#include "edginc/RatePointwise.hpp"

DRLIB_BEGIN_NAMESPACE
/** A Bond implementation where you input the cashflows */

class PRODUCTS_DLL ExplicitBond: public Bond,
                    public Theta::IShift
{
public:
    static CClassConstSP const TYPE;
    friend class ExplicitBondHelper;

    virtual void fieldsUpdated(const CFieldArray& fields);

    /** change the settlement of the bond if applicable */
    virtual void setSettlement(SettlementSP newSettle);

    void calculateLiborRates(/*bool isTheta*/);

    virtual CashFlowArraySP getOutstandingNotionals();

    virtual void   getMarket(const IModel* model, const MarketData* market);

    void validatePop2Object();

    virtual ~ExplicitBond() {}

    // methods derived from Bond

    /** get the accrued interest at a date */
    virtual double getAccruedAtDate(const DateTime &aiDate) const;
    
    /** returns only the cashFlows that occur after the start date */
    virtual CashFlowArraySP getCashFlows(const DateTime &startDate) const;

    /** returns all the cashFlows */
    virtual CashFlowArraySP getCashFlows() const;

    /** returns only the coupons that occur after the start date */
    virtual CashFlowArraySP getCoupons(const DateTime &startDate) const;

    /** returns all the coupons */
    virtual CashFlowArraySP getCoupons() const;

    virtual DateTimeArraySP getExCouponDates() const;
    
    virtual DateTime getMaturityDate() const;

    virtual DateTime getUnadjMaturityDate() const;
    
    virtual double getFaceValue() const;

    virtual double getNotional(const DateTime &date) const;

    virtual CashFlowArraySP getRedemptionPayments() const;

    virtual double getRedemption() const;

    virtual double getMaturityCoupon() const;

    /** required for Yield to Maturity calculation */
    virtual double priceFromYield(double yield, bool clean, DateTime yDate) const;

    /** get the equivalent bond which matures on the put date. Used for yield to put. */
    virtual Bond* getBondToPut(const DateTime &putDate, double putLevel) const;

    virtual DateTime getAccrualStartDate() const;

    virtual CashFlowArraySP getCouponRates();

    virtual CashFlowArraySP getDayCountFractions();

    // Theta
    bool sensShift(Theta* shift);
 
    // when does bond settle?
    virtual DateTime settles(const DateTime& tradeDate) const;

private:
    ExplicitBond();

    bool isZeroBond() const;

    DateTimeArraySP             resetDate;
    DateTimeArraySP             initialAccrualDate;
    DateTimeArraySP             finalAccrualDate;
    DateTimeArraySP             exCouponDate;
    DateTimeArraySP             paymentDate;
    DoubleArraySP               notional;
    DoubleArraySP               couponRate;
    CBoolArraySP                isFixed;
    DoubleArraySP               liborSpread;
    DoubleArraySP               couponCap;
    DoubleArraySP               couponFloor;
    DateTimeArraySP             redemptionDate;
    DoubleArraySP               redemption;
    StringArraySP               accrualDCC;
    StringArraySP               paymentDCC;
    MaturityPeriodSP            liborPeriod;
    string                      liborDCCString;
    YieldCurveWrapper           liborReferenceCurve;
    DateTime                    valueDate;
    int                         bondFrequency;
    SettlementSP                settle;           // Optional

    // internal fields
    DateTime                    datedDate;
    DoubleArraySP               rateToUse;
    CashFlowArraySP             cashFlows;
    DoubleArraySP               coupons;
    DoubleArraySP               forwards;
    DoubleArraySP               dayCountFracs;
    DayCountConventionSP        liborDCC;
    DayCountConventionArraySP   bondPaymentDCC;
    DayCountConventionArraySP   bondAccrualDCC;
};

typedef smartConstPtr<ExplicitBond> ExplicitBondConstSP;
typedef smartPtr<ExplicitBond>      ExplicitBondSP;

DRLIB_END_NAMESPACE
#endif
