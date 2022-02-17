/**
 * @file LiquiditySpreadPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_LIQUIDITYSPREADPOINTWISE_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/LiquiditySpreadPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

LiquiditySpreadPointwise::LiquiditySpreadPointwise(): CObject(TYPE) {}
LiquiditySpreadPointwise::~LiquiditySpreadPointwise() {}

static void LiquiditySpreadPointwise_load(CClassSP& clazz) {
    REGISTER(LiquiditySpreadPointwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<LiquiditySpreadPointwise>::iObject);
}

CClassConstSP const LiquiditySpreadPointwise::TYPE = CClass::registerClassLoadMethod("LiquiditySpreadPointwise", typeid(LiquiditySpreadPointwise), LiquiditySpreadPointwise_load);

RiskProperty_TYPES(LiquiditySpreadPointwise)
INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<LiquiditySpreadPointwise>);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
