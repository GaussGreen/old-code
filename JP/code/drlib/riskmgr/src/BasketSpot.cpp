/**
 * @file BasketSpot.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/BasketSpot.hpp"

DRLIB_BEGIN_NAMESPACE

BasketSpot::BasketSpot(): CObject(TYPE) {}
BasketSpot::~BasketSpot() {}

static void BasketSpot_load(CClassSP& clazz) {
    REGISTER(BasketSpot, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<BasketSpot>::iObject);
}

CClassConstSP const BasketSpot::TYPE = CClass::registerClassLoadMethod("BasketSpot", typeid(BasketSpot), BasketSpot_load);

RiskProperty_TYPES(BasketSpot)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
