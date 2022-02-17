//-------------------------------------------------------------
//
//   Group       : Credit Hybrids QR&D
//
//   Filename    : FixedCashFlow.hpp
//
//   Description : A fixed coupon that pays over a specified
//                 accrual period
//
//   Author      : Gordon Stephens
//
//   Date        : 8 June 2005
//-------------------------------------------------------------

#ifndef FIXEDCASHFLOW_HPP
#define FIXEDCASHFLOW_HPP

#include "edginc/RiskFreeCashFlow.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/CreditFeeLegCoupon.hpp"


DRLIB_BEGIN_NAMESPACE

class MARKET_DLL FixedCashFlow : 
                     public RiskFreeCashFlow,
                     public virtual IRestorableWithRespectTo<CreditFeeLegCoupon>
{
public:
    static CClassConstSP const TYPE;

    virtual bool                      isZero() const;
    virtual void                      getCriticalDates(DateTimeArray& dates) const;
    virtual const DateTime            getMaturityDate() const;
    virtual double                    getAmount(
                                        const ZeroCurve& curve, 
                                        const IRVolBase* volModelIR) const;

	virtual double               getAmount(IForwardRatePricerSP model) const;
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

    FixedCashFlow();
    FixedCashFlow(  DateTime             valueDate,
                    DateTime             payDate,
                    double               notional,
                    DateTime             accrueStart,
                    DateTime             accrueEnd,
                    DayCountConventionSP accrualDcc,
                    double               rate);

    ~FixedCashFlow();

    //CObject methods
    virtual void validatePop2Object();

    //// Implementing tweakable with respect to CreditFeelegCoupon 
    //// - a shift in coupon
    TweakOutcome sensShift(const PropertyTweak<CreditFeeLegCoupon>& shift);
    string sensName(const CreditFeeLegCoupon* shift) const;
    void sensRestore(const PropertyTweak<CreditFeeLegCoupon>& shift);

    /** Returns the accrual period (start and end dates) corresponding to 
     *  this cash flow. */
    virtual AccrualPeriodSP getAccrualPeriod() const;

private:
    //infrastructure support
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    //fields
    double                  notional;       // the size of the cashflow
    DateTime                accrueStart;    // the start of the accrual period
    DateTime                accrueEnd;      // the end of the accrual period
    DayCountConventionSP    accrualDcc;     // the accrual day count convention
    double                  rate;           // the fixed rate
};

DECLARE(FixedCashFlow)

DRLIB_END_NAMESPACE

#endif
