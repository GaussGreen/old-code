/**
 * @file ParSpreadUpfronts.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_PARSPREADUPFRONTS_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/ParSpreadUpfronts.hpp"

DRLIB_BEGIN_NAMESPACE

ParSpreadUpfronts::ParSpreadUpfronts(): CObject(TYPE) {}
ParSpreadUpfronts::~ParSpreadUpfronts() {}

static void ParSpreadUpfronts_load(CClassSP& clazz) {
    REGISTER(ParSpreadUpfronts, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<ParSpreadUpfronts>::iObject);
}

CClassConstSP const ParSpreadUpfronts::TYPE = CClass::registerClassLoadMethod("ParSpreadUpfronts", typeid(ParSpreadUpfronts), ParSpreadUpfronts_load);

RiskProperty_TYPES(ParSpreadUpfronts)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<ParSpreadUpfronts>);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
