//-------------------------------------------------------------
//
//   Group       : Credit Hybrids QR&D
//
//   Filename    : RiskFreeCashFlow.cpp
//
//   Description : Abstract base class for dealing with a
//                 variety of risk free cashflow descriptions
//
//   Author      : Gordon Stephens
//
//   Date        : 8 June 2005
//-------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RiskFreeCashFlow.hpp"

DRLIB_BEGIN_NAMESPACE

RiskFreeCashFlow::RiskFreeCashFlow()
    : BaseCashFlow(TYPE)
{
}

RiskFreeCashFlow::RiskFreeCashFlow(CClassConstSP clazz)
    : BaseCashFlow(clazz)
{
}

RiskFreeCashFlow::RiskFreeCashFlow(CClassConstSP clazz,
                                   DateTime valueDate,
                                   DateTime payDate)
    : BaseCashFlow(clazz, valueDate, payDate)
{
}

RiskFreeCashFlow::~RiskFreeCashFlow()
{
}

bool RiskFreeCashFlow::isRiskFree() const
{
    return true;
}

/** Returns known cash flow corresponding to a defaulted CDS. Since this
 * is a riskless cash flow, it will not be affected by the default */
CashFlow RiskFreeCashFlow::getKnownCashFlowGivenDefault(
    const DateTime&      defaultDeterminationDate,
    const DateTime&      accrualPaymentDate,
    IForwardRatePricerSP model) const
{
    CashFlow cf(getPayDate(), getAmount(model));
    return cf;
}


/** Returns the fraction of this cashflow accrued in case of default.
    Since these are riskfree cashflows, there is no concept of risky accrue
    period, so the complete cashflow is returned here. */
CashFlow RiskFreeCashFlow::getAccruedCashFlow(
    const DateTime&      valuationDate,
    IForwardRatePricerSP model) const
{
    return getKnownCashFlowGivenDefault(valuationDate, valuationDate, model);
}

bool RiskFreeCashFlow::isFixed() const {
    return false;
}

/** Invoked when Class is 'loaded' */
void RiskFreeCashFlow::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(RiskFreeCashFlow, clazz);
    SUPERCLASS(BaseCashFlow);
}

CClassConstSP const RiskFreeCashFlow::TYPE = CClass::registerClassLoadMethod(
    "RiskFreeCashFlow", typeid(RiskFreeCashFlow), load);


DEFINE_TEMPLATE_TYPE(RiskFreeCashFlowArray);


DRLIB_END_NAMESPACE
