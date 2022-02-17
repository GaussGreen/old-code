//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : CCMBaseCorrelation.cpp
//
//   Description : Composite Copula Model with Base Correlation
//
//   Author      : Antoine Gregoire
//
//   Date        : May 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/CCM.hpp" // for static CCM::ccmHasStochasticRecoveries
#include "edginc/CCMBaseCorrelation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Same as CreditMetricsModel::createLossCalculator() method but
    returns one which is derived from
    CreditMetricsLossCalculatorBase. This method won't make sense
    for certain derived types - or rather its meaning changes.
    The implementation here creates an appropriate CCM loss
    calculator */
CreditMetricsLossCalculatorBase* CCMBaseCorrelation::createLossCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    // use utility method in parent
    return createCCMLossCalculator(timeline, tranche, cpty);
}


CreditMetricsLossCalculatorBase* CCMBaseCorrelation::createRecoveredNotionalCalculatorBase(
    const DateTimeArray&           timeline,    /* (I) */
    CreditTrancheLossConfigConstSP tranche,     /* (I) */
    CounterPartyCreditConstSP      cpty) const  /* (I) */
{
    // use utility method in parent
    return createCCMRecoveredNotionalCalculator(timeline, tranche, cpty);
}


/** Indicates whether this model supports stochastic recovery rates
    AND there are any engine parameters for any names specifying so. */
const bool CCMBaseCorrelation::hasStochasticRecoveries(
    CreditTrancheLossConfigConstSP tranche) const
{
    return CCM::ccmHasStochasticRecoveries(tranche);
}


/** Invoked when Class is 'loaded' */
void CCMBaseCorrelation::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CCMBaseCorrelation, clazz);
    SUPERCLASS(CreditMetricsBaseCorrelation);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

/** Private constructor */
CCMBaseCorrelation::CCMBaseCorrelation(const CClassConstSP& clazz) : CreditMetricsBaseCorrelation(clazz) {}

/** Default constructor */
IObject* CCMBaseCorrelation::defaultConstructor(){
    return new CCMBaseCorrelation(CCMBaseCorrelation::TYPE);
}

CClassConstSP const CCMBaseCorrelation::TYPE = CClass::registerClassLoadMethod(
    "CCMBaseCorrelation", typeid(CCMBaseCorrelation), CCMBaseCorrelation::load);

bool CCMBaseCorrelationLoad() {
    return (CCMBaseCorrelation::TYPE != 0);
}

DRLIB_END_NAMESPACE

