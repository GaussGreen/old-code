/**
 * @file LegalBasisRelativePointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_LEGALBASISRELATIVEPOINTWISE_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/LegalBasisRelativePointwise.hpp"

DRLIB_BEGIN_NAMESPACE

LegalBasisRelativePointwise::LegalBasisRelativePointwise(): CObject(TYPE) {}
LegalBasisRelativePointwise::~LegalBasisRelativePointwise() {}

static void LegalBasisRelativePointwise_load(CClassSP& clazz) {
    REGISTER(LegalBasisRelativePointwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<LegalBasisRelativePointwise>::iObject);
}

CClassConstSP const LegalBasisRelativePointwise::TYPE = CClass::registerClassLoadMethod("LegalBasisRelativePointwise", typeid(LegalBasisRelativePointwise), LegalBasisRelativePointwise_load);

RiskProperty_TYPES(LegalBasisRelativePointwise)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<LegalBasisRelativePointwise>);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
