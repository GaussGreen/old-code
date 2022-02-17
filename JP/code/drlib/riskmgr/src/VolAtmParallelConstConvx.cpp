/**
 * @file VolAtmParallelConstConvx.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/VolAtmParallelConstConvx.hpp"

DRLIB_BEGIN_NAMESPACE

VolAtmParallelConstConvx::VolAtmParallelConstConvx(): CObject(TYPE) {}
VolAtmParallelConstConvx::~VolAtmParallelConstConvx() {}

static void VolAtmParallelConstConvx_load(CClassSP& clazz) {
    REGISTER(VolAtmParallelConstConvx, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<VolAtmParallelConstConvx>::iObject);
}

CClassConstSP const VolAtmParallelConstConvx::TYPE = CClass::registerClassLoadMethod("VolAtmParallelConstConvx", typeid(VolAtmParallelConstConvx), VolAtmParallelConstConvx_load);

RiskProperty_TYPES(VolAtmParallelConstConvx)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
