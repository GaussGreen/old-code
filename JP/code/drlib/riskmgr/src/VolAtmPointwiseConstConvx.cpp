/**
 * @file VolAtmPointwiseConstConvx.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/VolAtmPointwiseConstConvx.hpp"

DRLIB_BEGIN_NAMESPACE

VolAtmPointwiseConstConvx::VolAtmPointwiseConstConvx(): CObject(TYPE) {}
VolAtmPointwiseConstConvx::~VolAtmPointwiseConstConvx() {}

static void VolAtmPointwiseConstConvx_load(CClassSP& clazz) {
    REGISTER(VolAtmPointwiseConstConvx, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<VolAtmPointwiseConstConvx>::iObject);
}

CClassConstSP const VolAtmPointwiseConstConvx::TYPE = CClass::registerClassLoadMethod("VolAtmPointwiseConstConvx", typeid(VolAtmPointwiseConstConvx), VolAtmPointwiseConstConvx_load);

RiskProperty_TYPES(VolAtmPointwiseConstConvx)

DRLIB_END_NAMESPACE
