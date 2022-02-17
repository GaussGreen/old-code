//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : BondFloatNtl.hpp
//
//   Description : Bond with notional that floats
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : June 18, 2002
//
//
//----------------------------------------------------------------------------

#ifndef FLOATING_BONDFLOATNTL_HPP
#define FLOATING_BONDFLOATNTL_HPP

#include "edginc/Bond.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/RhoParallel.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Settlement.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/RateParallel.hpp"
#include "edginc/RatePointwise.hpp"

DRLIB_BEGIN_NAMESPACE
/** A Bond implementation where you input the cashflows */

class PRODUCTS_DLL BondFloatNtl: public Bond,
                    public Theta::IShift
{
public:
    static CClassConstSP const TYPE;
    friend class BondFloatNtlHelper;
    friend class FloatNtlAccreted;

    virtual ~BondFloatNtl();
    
    /** re-initialize the transcient field if the yield curve is tweaked */
    virtual void fieldsUpdated(const CFieldArray& fields);

    /** change the settlement of the bond if applicable */
    virtual void setSettlement(SettlementSP newSettle);

    /** get the accrued interest at a date */
    virtual double getAccruedAtDate(const DateTime &aiDate) const;

    /** returns only the cashFlows that occur after the start date */
    virtual CashFlowArraySP getCashFlows(const DateTime &startDate) const;
    
    /** returns all the cashFlows */
    virtual CashFlowArraySP getCashFlows() const;

    virtual CashFlowArraySP getExAdjCashFlows(const DateTime &startDate, YieldCurveConstSP discount) const;

    virtual DateTime getMaturityDate() const;
    
    virtual double getFaceValue() const;

    virtual double getNotional(const DateTime &date) const;

    virtual double getRedemption() const;
    
    virtual double getMaturityCoupon() const;

    virtual DateTimeArraySP getExCouponDates() const;

    /** returns only the cashFlows that occur after the start date */
    virtual CashFlowArraySP getCoupons(const DateTime &startDate) const;

    /** returns all the coupons */
    virtual CashFlowArraySP getCoupons() const;

    virtual void validatePop2Object();

    virtual double priceFromYield(double yield, bool clean, DateTime yDate) const;

    virtual Bond* getBondToPut(const DateTime &putDate, double putLevel) const;

    void    getMarket(const IModel* model, const MarketData* market);

    void calculateLiborRates();

    void initialize();

    // get all cashflows known today
    virtual CashFlowArraySP getKnownCashFlows() const;

    // Theta
    bool sensShift(Theta* shift);

    virtual DateTime getAccrualStartDate() const;

    // when does bond settle?
    virtual DateTime settles(const DateTime& tradeDate) const;

private:
    BondFloatNtl();
    virtual void throwIfNotValid() const;
    virtual void validateInputs() const;
    
    DateTime             valueDate;
    DateTime             maturityDate;
    double               faceValue;
    DateTimeArraySP      fixingDates;
    DoubleArraySP        rates;
    MaturityPeriodSP     liborPeriod;
    string               liborDCCString;
    double               liborSpread;
    YieldCurveWrapper    liborReferenceCurve;
    double               redemptionPct;            // Optional
    DateTime             datedDate;

    // cap/floor parameters
    bool                        capFloatingCoupons;
    double                      floatingCouponCap;
    bool                        floorFloatingCoupons;
    double                      floatingCouponFloor;


    SettlementSP         settle;      // Optional

    // derived
    bool                 isValid;
    DoubleArraySP        ratesProcessed;
    DoubleArraySP        accretionRates;
    DoubleArraySP        notionalsProcessed;
    DayCountConventionSP liborDCC;
};

typedef smartConstPtr<BondFloatNtl> BondFloatNtlConstSP;
typedef smartPtr<BondFloatNtl> BondFloatNtlSP;

DRLIB_END_NAMESPACE
#endif
