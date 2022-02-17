/**
 * @file Spot.cpp
 */

#include "edginc/config.hpp"
#define QLIB_SPOT_CPP
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/Spot.hpp"

DRLIB_BEGIN_NAMESPACE

Spot::Spot(): CObject(TYPE) {}
Spot::~Spot() {}

static void Spot_load(CClassSP& clazz) {
  REGISTER(Spot, clazz);
  SUPERCLASS(CObject);
  EMPTY_SHELL_METHOD(DefaultConstructor<Spot>::iObject);
}

CClassConstSP const Spot::TYPE = CClass::registerClassLoadMethod(
  "Spot", typeid(Spot), Spot_load);

RiskProperty_TYPES(Spot)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<Spot>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<Spot>);

DRLIB_END_NAMESPACE
