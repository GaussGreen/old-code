//-------------------------------------------------------------
//
//   Group       : Credit Hybrids QR&D
//
//   Filename    : RiskFreeCashFlow.hpp
//
//   Description : Abstract base class for dealing with a
//                 variety of risk free cashflow descriptions
//
//   Author      : Gordon Stephens
//
//   Date        : 8 June 2005
//-------------------------------------------------------------

#ifndef RiskFreeCashFlow_HPP
#define RiskFreeCashFlow_HPP

#include "edginc/config.hpp"
#include "edginc/AbstractCashFlow.hpp"
#include "edginc/DayCountConvention.hpp"


DRLIB_BEGIN_NAMESPACE


class ZeroCurve;
class IRVolBase;


class MARKET_DLL RiskFreeCashFlow : public BaseCashFlow
{
public:
    static CClassConstSP const TYPE;

    virtual double          getAmount(IForwardRatePricerSP model) const = 0;
    virtual double          getAmount(double rateOverride) const = 0;
    virtual double          getNotional() const = 0;
    virtual DateTime        getAccrueStart() const = 0;
    virtual DateTime        getAccrueEnd() const = 0;
    virtual DayCountConventionConstSP getAccrualDcc() const = 0;
    virtual bool            isAdhoc() const = 0;

    virtual bool            isRiskFree() const;

    /**
     * Enumerates critical dates; the critical points are the set of points 
     * that are required to correctly reprice the benchmark instruments.
     */
    virtual void            getCriticalDates(DateTimeArray& dates) const = 0;

    virtual const DateTime  getMaturityDate() const = 0;

    virtual bool            isZero() const = 0;
    virtual bool            isFixed() const;

    virtual double          getAmount(const ZeroCurve& curve, 
                                      const IRVolBase* volModelIR) const = 0;

    ~RiskFreeCashFlow();

    /** Returns known cash flow corresponding to a defaulted CDS. Since this
     * is a riskless cash flow, it will not be affected by the default */
    virtual CashFlow getKnownCashFlowGivenDefault(
        const DateTime&      defaultDeterminationDate,
        const DateTime&      accrualPaymentDate,
        IForwardRatePricerSP model) const;

    /** Returns the fraction of this cashflow accrued in case of default.
        Since these are riskfree cashflows, there is no concept of risky accrue
        period, so the complete cashflow is returned here. */
    virtual CashFlow getAccruedCashFlow(
        const DateTime&      valuationDate,
        IForwardRatePricerSP model) const;

protected:
    RiskFreeCashFlow();
    RiskFreeCashFlow(CClassConstSP clazz);
    RiskFreeCashFlow(CClassConstSP clazz,
                     DateTime valueDate,
                     DateTime payDate);

private:
    //infrastructure support
    static void load(CClassSP& clazz);
};

DECLARE(RiskFreeCashFlow)

DRLIB_END_NAMESPACE

#endif
