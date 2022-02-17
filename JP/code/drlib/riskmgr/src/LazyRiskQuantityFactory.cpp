/**
 * @file LazyRiskQuantityFactory.cpp
 */

#include "edginc/config.hpp"
#include "edginc/IHypothesis.hpp"
#include "edginc/IRiskQuantityFactory.hpp"
#include "edginc/LazyRiskQuantityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

LazyRiskQuantityFactory::LazyRiskQuantityFactory(IHypothesisConstSP hypothesis,
                                                 IRiskQuantityFactorySP greek):
    CObject(TYPE),
    hypothesis(hypothesis),
    greek(greek)
{}

LazyRiskQuantityFactory::~LazyRiskQuantityFactory() {}

static IObject* defaultLazyRiskQuantityFactory() {
    return new LazyRiskQuantityFactory(IHypothesisSP(), IRiskQuantityFactorySP());
}

void LazyRiskQuantityFactory::load(CClassSP& clazz) {
    REGISTER(LazyRiskQuantityFactory, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultLazyRiskQuantityFactory);
    FIELD(hypothesis, "hypothesis");
    FIELD(greek, "greek");
}

CClassConstSP const LazyRiskQuantityFactory::TYPE = CClass::registerClassLoadMethod(
    "LazyRiskQuantityFactory", typeid(LazyRiskQuantityFactory), load);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
