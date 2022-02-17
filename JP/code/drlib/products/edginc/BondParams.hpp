//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : BondParams.hpp
//
//   Description : Full Bond implementation from parameters
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : December 12, 2001
//
//
//----------------------------------------------------------------------------

#ifndef BONDPARAMS_HPP
#define BONDPARAMS_HPP

#include "edginc/Bond.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Settlement.hpp"

DRLIB_BEGIN_NAMESPACE
/** Full Bond implementation from parameters */

class PRODUCTS_DLL BondParams:public Bond {
public:
    static CClassConstSP const TYPE;
    friend class BondParamsHelper;
    friend class CFAddin;
    friend class ConvBond;
    friend class OptOnConvBond;
    friend class BondFloatNtl;
    friend class FixedYieldParameters;
 
    virtual ~BondParams();

    /** change the settlement of the bond if applicable */
    virtual void setSettlement(SettlementSP newSettle);
     
    /** get the accrued interest at a date */
    virtual double getAccruedAtDate(const DateTime &aiDate) const;

    /** returns only the cashFlows that occur after the start date */
    virtual CashFlowArraySP getCashFlows(const DateTime &startDate) const;
    
    /** returns all the cashFlows */
    virtual CashFlowArraySP getCashFlows() const;

    /** return a reference to the cash flows */
    virtual CashFlowArrayConstSP getCashFlowRef() const;

    virtual DateTime getMaturityDate() const;

    virtual DateTime getUnadjMaturityDate() const;

    virtual DateTimeArraySP getExCouponDates() const;
    
    virtual double getFaceValue() const;

    virtual double getNotional(const DateTime &date) const;

    virtual double getRedemption() const;
    
    virtual double getMaturityCoupon() const;

    virtual void   validatePop2Object();

    void           getMarket(const IModel* model, const MarketData* market);

    virtual double priceFromYield(double yield, bool clean, DateTime yDate) const;

    virtual Bond* getBondToPut(const DateTime &putDate, double putLevel) const;

    virtual double   getBondCouponRate() const;

    virtual DateTime getAccrualStartDate() const;

    virtual CashFlowArraySP getRedemptionPayments() const;

    // when does bond settle?
    virtual DateTime settles(const DateTime& tradeDate) const;

private:
    BondParams();
    virtual void throwIfNotValid() const;
    virtual void initialize();
    virtual void validateInputs() const;
    
    // inputs
    double     faceValue;
    double     redemptionPct;       // Optional
    double     couponPct;
    int        frequency;
    DateTime   maturityDate;
    DateTime   datedDate;
    string     dayCountConvString; // Optional
    DateTime   firstCouponDate;    // Optional
    bool       oddLastShort;       // Optional
    bool       endOfMonthAdj;      // Optional
    bool       eomIgnoreLeapYear;  // Optional
    string     badDayConvString;   // Optional
    HolidayWrapper  hols;
    int        exDivDays;          // Optional
    string     exDivRule;          // Optional
    double     taxRate;            // Optional
    
    bool       isAnOID;            // Optional
    double     yieldForOID;        // Optional
    
    SettlementSP settle;           // Optional

    // derived
    DateTimeArraySP      exCouponDates;
    DateTimeArraySP      schedCouponDates; // scheduled pay dates unadjusted for bad days
    CashFlowArraySP      cashFlows;
    bool                 isValid;
    string               errorMessage;
    BadDayConventionSP   badDayConv;
    DayCountConventionSP accruedDCC;

    DateTime             pseudoStartDate;
    DateTime             pseudoMatDate;

    bool                 oddFirstLong;
    bool                 oddFirst;
    bool                 oddLastLong;
    bool                 oddLast;
    CashFlowArraySP      pseudoFrontFlows;
    CashFlowArraySP      pseudoBackFlows;
};

typedef smartConstPtr<BondParams> BondParamsConstSP;
typedef smartPtr<BondParams> BondParamsSP;

DRLIB_END_NAMESPACE
#endif
