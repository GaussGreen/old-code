//----------------------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR&D
//
//   Filename    : AbstractCashFlow.hpp
//
//   Description : Abstract base class for dealing with a variety of cashflow descriptions
//
//   Author      : Gordon Stephens
//
//   Date        : 8 June 2005
//-----------------------------------------------------------------------------------------

#ifndef AbstractCashFlow_HPP
#define AbstractCashFlow_HPP

#include "edginc/config.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/AccrualPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(AbstractCashFlow)

class MARKET_DLL AbstractCashFlow : public CObject,
                                    public virtual IGetMarket
{
public:
    static CClassConstSP const TYPE;

    /** Utility method allowing conversion into date/amount type structures */
    static CashFlowArraySP asCashFlows(const AbstractCashFlowArrayConstSP abstractCashflows,
                                       const IForwardRatePricerSP         model);

    virtual DateTime getPayDate()  const = 0;
    virtual double   getAmount(IForwardRatePricerSP model)   const = 0;
    virtual double   getNotional() const = 0;
    virtual void     setNotional(double newNotional) = 0;
    virtual double   getRate() const = 0; //implicitly a fixed rate
    virtual void     setRate(double newRate) = 0; //implicitly a fixed rate
    
    /** Returns cash flow amount using given rate override (if relevant) */
    virtual double getAmount(double rateOverride) const = 0;

    /** Returns known cash flow corresponding to a defaulted CDS, computing
     * accrued as required */
    virtual CashFlow getKnownCashFlowGivenDefault(
        const DateTime&      defaultDeterminationDate,
        const DateTime&      defaultPaymentDate,
        IForwardRatePricerSP model) const = 0;

    /** When to stop tweaking for Yield Curve type tweaks */
    virtual DateTime lastYCSensDate() const = 0;
    /** Populate market data structures */
    virtual void getMarket(const IModel* model, const MarketData* market) = 0;
    
    /** Says if this cash flow is actually an adhoc cash flow */
    virtual bool isAdhoc() const = 0;

    /** Says if this cash flow is risk free */
    virtual bool isRiskFree() const = 0;

    /** Returns the accrual period (start and end dates) corresponding to 
     * this cash flow. */
    virtual AccrualPeriodSP getAccrualPeriod() const = 0;

    /** Returns the fraction of this cashflow accrued in case of default.
        Adhoc cashflows do not accrue (since they have no accrue period) and
        an empty defaultDeterminationDate is interpreted as meaning that the 
        name did not default, so the complete cashflow is returned. */
    virtual CashFlow getAccruedCashFlow(
        const DateTime&      valuationDate,
        IForwardRatePricerSP model) const = 0;

    ~AbstractCashFlow();

protected:
    AbstractCashFlow(CClassConstSP clazz);

private:
    //infrastructure support
    static void load(CClassSP& clazz);
};

//typedef smartPtr<AbstractCashFlow> AbstractCashFlowSP;
//typedef smartConstPtr<AbstractCashFlow> AbstractCashFlowConstSP;
//typedef array<AbstractCashFlowSP, AbstractCashFlow>  AbstractCashFlowArray;
//typedef smartPtr<AbstractCashFlowArray> AbstractCashFlowArraySP;
//typedef smartConstPtr<AbstractCashFlowArray> AbstractCashFlowArrayConstSP;


// a base class for cashflows that holds common parameters
// and provides implementations of some of the basic AbstractCashFlow methods
class MARKET_DLL BaseCashFlow : public AbstractCashFlow
{
public:
    static CClassConstSP const TYPE;

    virtual DateTime getPayDate()  const;

    /** When to stop tweaking for Yield Curve type tweaks */
    virtual DateTime lastYCSensDate() const;
    
    /** Populate market data structures */
    virtual void getMarket(const IModel* model, const MarketData* market);

    ~BaseCashFlow();

protected:
    BaseCashFlow(CClassConstSP clazz);

    BaseCashFlow(CClassConstSP clazz,
                 DateTime      valueDate,
                 DateTime      payDate);
    //fields
    DateTime    valueDate; //from the market
    DateTime    payDate;   //when the cashflow is paid

private:
    //infrastructure support
    static void load(CClassSP& clazz);

};

DRLIB_END_NAMESPACE

#endif
