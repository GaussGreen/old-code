/**
 * @file CRVolParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/CRVolParallel.hpp"

DRLIB_BEGIN_NAMESPACE

CRVolParallel::CRVolParallel(): CObject(TYPE) {}
CRVolParallel::~CRVolParallel() {}

static void CRVolParallel_load(CClassSP& clazz) {
    REGISTER(CRVolParallel, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<CRVolParallel>::iObject);
}

CClassConstSP const CRVolParallel::TYPE = CClass::registerClassLoadMethod("CRVolParallel", typeid(CRVolParallel), CRVolParallel_load);

RiskProperty_TYPES(CRVolParallel)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
