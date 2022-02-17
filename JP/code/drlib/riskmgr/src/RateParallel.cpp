/**
 * @file RateParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_RATEPARALLEL_CPP
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/RateParallel.hpp"

DRLIB_BEGIN_NAMESPACE

RateParallel::RateParallel(): CObject(TYPE) {}
RateParallel::~RateParallel() {}

static void RateParallel_load(CClassSP& clazz) {
    REGISTER(RateParallel, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<RateParallel>::iObject);
}

CClassConstSP const RateParallel::TYPE = CClass::registerClassLoadMethod("RateParallel", typeid(RateParallel), RateParallel_load);

RiskProperty_TYPES(RateParallel)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<RateParallel>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<RateParallel>);

DRLIB_END_NAMESPACE
