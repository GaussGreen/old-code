/**
 * @file VolParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_VOLPARALLEL_CPP
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/VolParallel.hpp"

DRLIB_BEGIN_NAMESPACE

VolParallel::VolParallel(): CObject(TYPE) {}
VolParallel::~VolParallel() {}

static void VolParallel_load(CClassSP& clazz) {
    REGISTER(VolParallel, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<VolParallel>::iObject);
}

CClassConstSP const VolParallel::TYPE = CClass::registerClassLoadMethod("VolParallel", typeid(VolParallel), VolParallel_load);

RiskProperty_TYPES(VolParallel)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<VolParallel>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<VolParallel>);

DRLIB_END_NAMESPACE
