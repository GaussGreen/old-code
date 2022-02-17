//-------------------------------------------------------------
//
//   Group       : Credit Hybrids QR&D
//
//   Filename    : AdhocCashFlow.hpp
//
//   Description : An amount paid on a specific date
//
//   Author      : Gordon Stephens
//
//   Date        : 8 June 2005
//-------------------------------------------------------------

#ifndef ADHOCCASHFLOW_HPP
#define ADHOCCASHFLOW_HPP

#include "edginc/RiskFreeCashFlow.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL AdhocCashFlow : public RiskFreeCashFlow
{
public:
    static CClassConstSP const TYPE;

    //RiskFreeCashFlow methods
    virtual bool                 isZero() const;
    virtual bool                 isFixed() const;
    virtual double               getAmount(IForwardRatePricerSP model)   const;
    virtual double               getNotional() const;
    virtual void                 setNotional(double newNotional);
    virtual double               getRate() const; //implicitly a fixed rate
    virtual void                 setRate(double newRate); //implicitly a fixed rate
    virtual DateTime             getAccrueStart() const;
    virtual DateTime             getAccrueEnd() const;
    virtual DayCountConventionConstSP getAccrualDcc() const;
    virtual bool                 isAdhoc() const;

    /** Returns cash flow amount using given rate override (if relevant) */
    virtual double getAmount(double rateOverride) const;

    AdhocCashFlow();
    AdhocCashFlow(DateTime valueDate, DateTime payDate, double amount);
    AdhocCashFlow(DateTime payDate, double amount);

    ~AdhocCashFlow();
    virtual double               getAmount(const ZeroCurve& curve, 
                                           const IRVolBase* volModelIR) const;
    virtual void                 getCriticalDates(DateTimeArray& dates) const;
    virtual const DateTime       getMaturityDate() const;

    /** Returns the accrual period (start and end dates) corresponding to 
     *  this cash flow. */
    virtual AccrualPeriodSP getAccrualPeriod() const;


private:
    //infrastructure support
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    //fields
    double  amount; //the amount of the cashflow
};

DECLARE(AdhocCashFlow)

DRLIB_END_NAMESPACE

#endif
