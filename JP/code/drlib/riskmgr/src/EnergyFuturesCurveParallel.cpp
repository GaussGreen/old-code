/**
 * @file EnergyFuturesCurveParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/EnergyFuturesCurveParallel.hpp"

DRLIB_BEGIN_NAMESPACE

EnergyFuturesCurveParallel::EnergyFuturesCurveParallel(): CObject(TYPE) {}
EnergyFuturesCurveParallel::~EnergyFuturesCurveParallel() {}

static void EnergyFuturesCurveParallel_load(CClassSP& clazz) {
    REGISTER(EnergyFuturesCurveParallel, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<EnergyFuturesCurveParallel>::iObject);
}

CClassConstSP const EnergyFuturesCurveParallel::TYPE = CClass::registerClassLoadMethod("EnergyFuturesCurveParallel", typeid(EnergyFuturesCurveParallel), EnergyFuturesCurveParallel_load);

RiskProperty_TYPES(EnergyFuturesCurveParallel)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
