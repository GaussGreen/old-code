/**
 * @file IRiskQuantityFactory.cpp
 */

#include "edginc/config.hpp"
#include "edginc/IRiskQuantityFactory.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/RiskMapping.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/LazyRiskQuantityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
IRiskQuantityFactory::IRiskQuantityFactory(){}
IRiskQuantityFactory::~IRiskQuantityFactory(){}

CClassConstSP const IRiskQuantityFactory::TYPE = CClass::registerInterfaceLoadMethod(
    "IRiskQuantityFactory", typeid(IRiskQuantityFactory), 0);

DEFINE_TEMPLATE_TYPE(IRiskQuantityFactoryArray);

struct ExplicitRiskQuantityFactory: CObject, virtual IRiskQuantityFactory {

    static CClassConstSP const TYPE;
    static IObject* emptyShell();
    static void load(CClassSP& clazz);

    NamedRiskQuantityArraySP rqs; // $required

    ExplicitRiskQuantityFactory(
            NamedRiskQuantityArraySP rqs = NamedRiskQuantityArraySP()):
        CObject(TYPE),
        rqs(rqs)
    {}

    NamedRiskQuantityArraySP riskQuantities(
            MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const {
        return rqs;
    }

    LazyRiskQuantityFactoryArraySP lazies(MultiTweakGroupConstSP world) const {
        return LazyRiskQuantityFactoryArraySP();
    }
};

IRiskQuantityFactorySP fromRiskQuantities(NamedRiskQuantityArraySP rqs) {
    return IRiskQuantityFactorySP(new ExplicitRiskQuantityFactory(rqs));
}

IObject* ExplicitRiskQuantityFactory::emptyShell() {
  return new ExplicitRiskQuantityFactory();
}

void ExplicitRiskQuantityFactory::load(CClassSP& clazz) {
  REGISTER(ExplicitRiskQuantityFactory, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IRiskQuantityFactory);
  FIELD(rqs, "rqs");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const ExplicitRiskQuantityFactory::TYPE = CClass::registerClassLoadMethod(
  "ExplicitRiskQuantityFactory", typeid(ExplicitRiskQuantityFactory), ExplicitRiskQuantityFactory::load);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
