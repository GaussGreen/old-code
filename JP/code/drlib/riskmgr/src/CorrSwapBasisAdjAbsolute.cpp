/**
 * @file CorrSwapBasisAdjAbsolute.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/CorrSwapBasisAdjAbsolute.hpp"

DRLIB_BEGIN_NAMESPACE

CorrSwapBasisAdjAbsolute::CorrSwapBasisAdjAbsolute(): CObject(TYPE) {}
CorrSwapBasisAdjAbsolute::~CorrSwapBasisAdjAbsolute() {}

static void CorrSwapBasisAdjAbsolute_load(CClassSP& clazz) {
    REGISTER(CorrSwapBasisAdjAbsolute, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<CorrSwapBasisAdjAbsolute>::iObject);
}

CClassConstSP const CorrSwapBasisAdjAbsolute::TYPE = CClass::registerClassLoadMethod("CorrSwapBasisAdjAbsolute", typeid(CorrSwapBasisAdjAbsolute), CorrSwapBasisAdjAbsolute_load);

RiskProperty_TYPES(CorrSwapBasisAdjAbsolute)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
