#include "edginc/config.hpp"
#include "edginc/CondLossDistributionsGenRisklessKey.hpp"
#include "edginc/DiscreteDistribution.hpp"

DRLIB_BEGIN_NAMESPACE

/** Constructor */
CondLossDistributionsGenRisklessKey::CondLossDistributionsGenRisklessKey():
    CObject(TYPE)
{
    // Builds condLossDistribution once for all
    DoubleArraySP values(new DoubleArray(1));
    DoubleArraySP probas(new DoubleArray(1));
    (*values)[0] = 0.0;
    (*probas)[0] = 1.0;
    condLossDistribution.reset(
        new DiscreteDistribution(values, probas));
}

/** [Implements ICondLossDistributionsGenKey] */
IDistribution1DConstSP CondLossDistributionsGenRisklessKey::conditionalLossDistribution(
    IMarketFactorValueConstSP marketFactorValue) const
{
    return condLossDistribution;            
}

double CondLossDistributionsGenRisklessKey::conditionalSurvProb(
    IMarketFactorValueConstSP marketFactorValue) const
{
    return 1.0; // (riskless)
}

void CondLossDistributionsGenRisklessKey::load(CClassSP& clazz)
{
    clazz->setPrivate();
    REGISTER(CondLossDistributionsGenRisklessKey, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICondLossDistributionsGenKey);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD_NO_DESC(condLossDistribution);
}

/** TYPE (for reflection) */        
CClassConstSP const CondLossDistributionsGenRisklessKey::TYPE =
CClass::registerClassLoadMethod(
    "CondLossDistributionsGenRisklessKey",
    typeid(CondLossDistributionsGenRisklessKey),
    CondLossDistributionsGenRisklessKey::load);

IObject* CondLossDistributionsGenRisklessKey::defaultConstructor()
{
    return new CondLossDistributionsGenRisklessKey();
}

DRLIB_END_NAMESPACE
