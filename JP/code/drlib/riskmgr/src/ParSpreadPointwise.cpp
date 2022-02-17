/**
 * @file ParSpreadPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_PARSPREADPOINTWISE_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/ParSpreadPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

ParSpreadPointwise::ParSpreadPointwise(): CObject(TYPE) {}
ParSpreadPointwise::~ParSpreadPointwise() {}

static void ParSpreadPointwise_load(CClassSP& clazz) {
    REGISTER(ParSpreadPointwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<ParSpreadPointwise>::iObject);
}

CClassConstSP const ParSpreadPointwise::TYPE = CClass::registerClassLoadMethod("ParSpreadPointwise", typeid(ParSpreadPointwise), ParSpreadPointwise_load);

RiskProperty_TYPES(ParSpreadPointwise)
INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<ParSpreadPointwise>);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
