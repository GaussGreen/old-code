/**
 * @file SimpleTweakNameResolver.cpp
 */

#include "edginc/config.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/SimpleTweakNameResolver.hpp"

DRLIB_BEGIN_NAMESPACE

SimpleTweakNameResolver::SimpleTweakNameResolver(OutputNameConstSP name):
    CObject(TYPE),
    name(name)
{}

SimpleTweakNameResolver::~SimpleTweakNameResolver() {}

OutputNameConstSP SimpleTweakNameResolver::getMarketDataName() const {
    return name;
}

bool SimpleTweakNameResolver::nameMatches(const OutputName& name,
                                          IObjectConstSP obj) {
    return name.equals(MarketObjectConstSP::dynamicCast(obj)->getName());
}

static IObject* defaultSimpleTweakNameResolver() {
    return new SimpleTweakNameResolver(OutputNameConstSP());
}

void SimpleTweakNameResolver::load(CClassSP& clazz) {
    REGISTER(SimpleTweakNameResolver, clazz);
    SUPERCLASS(CObject);
    //IMPLEMENTS(ITweakNameResolver);
    EMPTY_SHELL_METHOD(defaultSimpleTweakNameResolver);
    FIELD(name, "name");
}

CClassConstSP const SimpleTweakNameResolver::TYPE = CClass::registerClassLoadMethod(
    "SimpleTweakNameResolver", typeid(SimpleTweakNameResolver), load);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
