/**
 * @file EnergyFuturesCurvePointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_ENERGYFUTURESCURVEPOINTWISE_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/EnergyFuturesCurvePointwise.hpp"

DRLIB_BEGIN_NAMESPACE

EnergyFuturesCurvePointwise::EnergyFuturesCurvePointwise(): CObject(TYPE) {}
EnergyFuturesCurvePointwise::~EnergyFuturesCurvePointwise() {}

static void EnergyFuturesCurvePointwise_load(CClassSP& clazz) {
    REGISTER(EnergyFuturesCurvePointwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<EnergyFuturesCurvePointwise>::iObject);
}

CClassConstSP const EnergyFuturesCurvePointwise::TYPE = CClass::registerClassLoadMethod("EnergyFuturesCurvePointwise", typeid(EnergyFuturesCurvePointwise), EnergyFuturesCurvePointwise_load);

RiskProperty_TYPES(EnergyFuturesCurvePointwise)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<EnergyFuturesCurvePointwise>);

DRLIB_END_NAMESPACE

// Local Variables: 
// c-basic-offset:
// End:
