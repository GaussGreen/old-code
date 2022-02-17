/**
 * @file CRSpotVolParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/CRSpotVolParallel.hpp"

DRLIB_BEGIN_NAMESPACE

CRSpotVolParallel::CRSpotVolParallel(): CObject(TYPE) {}
CRSpotVolParallel::~CRSpotVolParallel() {}

static void CRSpotVolParallel_load(CClassSP& clazz) {
    REGISTER(CRSpotVolParallel, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<CRSpotVolParallel>::iObject);
}

CClassConstSP const CRSpotVolParallel::TYPE = CClass::registerClassLoadMethod("CRSpotVolParallel", typeid(CRSpotVolParallel), CRSpotVolParallel_load);

RiskProperty_TYPES(CRSpotVolParallel)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
