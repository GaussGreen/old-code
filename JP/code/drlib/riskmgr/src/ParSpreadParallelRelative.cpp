/**
 * @file ParSpreadParallelRelative.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/ParSpreadParallelRelative.hpp"

DRLIB_BEGIN_NAMESPACE

ParSpreadParallelRelative::ParSpreadParallelRelative(): CObject(TYPE) {}
ParSpreadParallelRelative::~ParSpreadParallelRelative() {}

static void ParSpreadParallelRelative_load(CClassSP& clazz) {
    REGISTER(ParSpreadParallelRelative, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<ParSpreadParallelRelative>::iObject);
}

CClassConstSP const ParSpreadParallelRelative::TYPE = CClass::registerClassLoadMethod("ParSpreadParallelRelative", typeid(ParSpreadParallelRelative), ParSpreadParallelRelative_load);

RiskProperty_TYPES(ParSpreadParallelRelative)

DRLIB_END_NAMESPACE
