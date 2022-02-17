/**
 * @file CorrSwapSamplingAdjAbsolute.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/CorrSwapSamplingAdjAbsolute.hpp"

DRLIB_BEGIN_NAMESPACE

CorrSwapSamplingAdjAbsolute::CorrSwapSamplingAdjAbsolute(): CObject(TYPE) {}
CorrSwapSamplingAdjAbsolute::~CorrSwapSamplingAdjAbsolute() {}

static void CorrSwapSamplingAdjAbsolute_load(CClassSP& clazz) {
    REGISTER(CorrSwapSamplingAdjAbsolute, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<CorrSwapSamplingAdjAbsolute>::iObject);
}

CClassConstSP const CorrSwapSamplingAdjAbsolute::TYPE = CClass::registerClassLoadMethod("CorrSwapSamplingAdjAbsolute", typeid(CorrSwapSamplingAdjAbsolute), CorrSwapSamplingAdjAbsolute_load);

RiskProperty_TYPES(CorrSwapSamplingAdjAbsolute)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
