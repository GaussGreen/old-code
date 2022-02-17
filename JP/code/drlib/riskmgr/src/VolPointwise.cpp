/**
 * @file VolPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_VOLPOINTWISE_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/VolPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

VolPointwise::VolPointwise(): CObject(TYPE) {}
VolPointwise::~VolPointwise() {}

static void VolPointwise_load(CClassSP& clazz) {
    REGISTER(VolPointwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<VolPointwise>::iObject);
}

CClassConstSP const VolPointwise::TYPE = CClass::registerClassLoadMethod("VolPointwise", typeid(VolPointwise), VolPointwise_load);

RiskProperty_TYPES(VolPointwise)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<VolPointwise>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<VolPointwise>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<VolPointwise>);
DRLIB_END_NAMESPACE
