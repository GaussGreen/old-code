/**
 * @file ParSpreadParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/ParSpreadParallel.hpp"

DRLIB_BEGIN_NAMESPACE

ParSpreadParallel::ParSpreadParallel(): CObject(TYPE) {}
ParSpreadParallel::~ParSpreadParallel() {}

static void ParSpreadParallel_load(CClassSP& clazz) {
    REGISTER(ParSpreadParallel, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<ParSpreadParallel>::iObject);
}

CClassConstSP const ParSpreadParallel::TYPE = CClass::registerClassLoadMethod("ParSpreadParallel", typeid(ParSpreadParallel), ParSpreadParallel_load);

RiskProperty_TYPES(ParSpreadParallel)

DRLIB_END_NAMESPACE
