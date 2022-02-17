/**
 * @file CreditIndexBasisPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_CREDITINDEXBASISPOINTWISE_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/CreditIndexBasisPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

CreditIndexBasisPointwise::CreditIndexBasisPointwise(): CObject(TYPE) {}
CreditIndexBasisPointwise::~CreditIndexBasisPointwise() {}

static void CreditIndexBasisPointwise_load(CClassSP& clazz) {
    REGISTER(CreditIndexBasisPointwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<CreditIndexBasisPointwise>::iObject);
}

CClassConstSP const CreditIndexBasisPointwise::TYPE = CClass::registerClassLoadMethod("CreditIndexBasisPointwise", typeid(CreditIndexBasisPointwise), CreditIndexBasisPointwise_load);

RiskProperty_TYPES(CreditIndexBasisPointwise)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<CreditIndexBasisPointwise>);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
