/**
 * @file CRVolPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/ExpiryPair.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/CRVolPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

CRVolPointwise::CRVolPointwise(): CObject(TYPE) {}
CRVolPointwise::~CRVolPointwise() {}

static void CRVolPointwise_load(CClassSP& clazz) {
    REGISTER(CRVolPointwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<CRVolPointwise>::iObject);
}

CClassConstSP const CRVolPointwise::TYPE = CClass::registerClassLoadMethod("CRVolPointwise", typeid(CRVolPointwise), CRVolPointwise_load);

RiskProperty_TYPES(CRVolPointwise)
 
DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
