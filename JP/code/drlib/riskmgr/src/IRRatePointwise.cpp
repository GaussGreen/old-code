/**
 * @file IRRatePointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_IRRATEPOINTWISE_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/IRRatePointwise.hpp"

DRLIB_BEGIN_NAMESPACE

IRRatePointwise::IRRatePointwise(): CObject(TYPE) {}
IRRatePointwise::~IRRatePointwise() {}

static void IRRatePointwise_load(CClassSP& clazz) {
    REGISTER(IRRatePointwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<IRRatePointwise>::iObject);
}

CClassConstSP const IRRatePointwise::TYPE = CClass::registerClassLoadMethod("IRRatePointwise", typeid(IRRatePointwise), IRRatePointwise_load);

RiskProperty_TYPES(IRRatePointwise)
INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<IRRatePointwise>);

DRLIB_END_NAMESPACE
