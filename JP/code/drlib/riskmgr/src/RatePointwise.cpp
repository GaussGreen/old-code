/**
 * @file RatePointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_RATEPOINTWISE_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/RatePointwise.hpp"

DRLIB_BEGIN_NAMESPACE

RatePointwise::RatePointwise(): CObject(TYPE) {}
RatePointwise::~RatePointwise() {}

static void RatePointwise_load(CClassSP& clazz) {
    REGISTER(RatePointwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<RatePointwise>::iObject);
}

CClassConstSP const RatePointwise::TYPE = CClass::registerClassLoadMethod("RatePointwise", typeid(RatePointwise), RatePointwise_load);

RiskProperty_TYPES(RatePointwise)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<RatePointwise>);

INSTANTIATE_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<RatePointwise>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<RatePointwise>);

DRLIB_END_NAMESPACE
