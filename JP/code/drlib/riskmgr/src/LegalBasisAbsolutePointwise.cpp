/**
 * @file LegalBasisAbsolutePointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_LEGALBASISABSOLUTEPOINTWISE_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/LegalBasisAbsolutePointwise.hpp"

DRLIB_BEGIN_NAMESPACE

LegalBasisAbsolutePointwise::LegalBasisAbsolutePointwise(): CObject(TYPE) {}
LegalBasisAbsolutePointwise::~LegalBasisAbsolutePointwise() {}

static void LegalBasisAbsolutePointwise_load(CClassSP& clazz) {
    REGISTER(LegalBasisAbsolutePointwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<LegalBasisAbsolutePointwise>::iObject);
}

CClassConstSP const LegalBasisAbsolutePointwise::TYPE = CClass::registerClassLoadMethod("LegalBasisAbsolutePointwise", typeid(LegalBasisAbsolutePointwise), LegalBasisAbsolutePointwise_load);

RiskProperty_TYPES(LegalBasisAbsolutePointwise)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<LegalBasisAbsolutePointwise>);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
