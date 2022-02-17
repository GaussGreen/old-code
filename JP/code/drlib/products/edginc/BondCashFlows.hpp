//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : BondCashFlows.hpp
//
//   Description : Bond implementation where you input the cash flows
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : November 16, 2001
//
//
//----------------------------------------------------------------------------

#ifndef BONDCASHFLOWS_HPP
#define BONDCASHFLOWS_HPP

#include "edginc/Bond.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/Settlement.hpp"


DRLIB_BEGIN_NAMESPACE
/** A Bond implementation where you input the cashflows */

class PRODUCTS_DLL BondCashFlows:public Bond {
public:
    static CClassConstSP const TYPE;
    friend class BondCashFlowsHelper;
    friend class FloatingBond;
    
    virtual ~BondCashFlows();
    
    /** change the settlement of the bond if applicable */
    virtual void setSettlement(SettlementSP newSettle);

    void getMarket(const IModel* model, const MarketData* market);

    /** get the accrued interest at a date */
    virtual double getAccruedAtDate(const DateTime &aiDate) const;

    /** returns only the cashFlows that occur after the start date */
    virtual CashFlowArraySP getCashFlows(const DateTime &startDate) const;
    
    /** returns all the cashFlows */
    virtual CashFlowArraySP getCashFlows() const;

    virtual DateTime getMaturityDate() const;

    virtual DateTime getUnadjMaturityDate() const;

    virtual DateTimeArraySP getExCouponDates() const;
    
    virtual double getFaceValue() const;

    virtual double getRedemption() const;

    virtual double getMaturityCoupon() const;

    virtual void validatePop2Object();

    virtual double priceFromYield(double yield, bool clean, DateTime yDate) const;

    virtual Bond* getBondToPut(const DateTime &putDate, double putLevel) const;

    virtual DateTime getAccrualStartDate() const;

    // when does bond settle?
    virtual DateTime settles(const DateTime& tradeDate) const;

private:
    BondCashFlows();
    virtual void throwIfNotValid() const;
    virtual void validateInputs() const;
    
    double               faceValue;
    CashFlowArraySP      cashFlows;
    DateTime             datedDate;
    string               dayCountConvString;       // Optional
    double               redemptionPct;            // Optional
    bool                 lastCFIncludesRedemption; // Optional
    DateTimeArraySP      exCouponDates;             // Optional
    SettlementSP         settle;           // Optional

    // derived
    bool                 isValid;
    string               errorMessage;
    CashFlowArraySP      cashFlowsProcessed;
    DateTimeArraySP      exCouponDatesProcessed;
    DayCountConventionSP accruedDCC;
};

typedef smartConstPtr<BondCashFlows> BondCashFlowsConstSP;
typedef smartPtr<BondCashFlows> BondCashFlowsSP;

DRLIB_END_NAMESPACE
#endif
