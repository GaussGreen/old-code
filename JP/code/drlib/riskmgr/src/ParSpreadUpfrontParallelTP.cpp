//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : ParSpreadUpfrontParallelTP.cpp
//
//   Description : ParSpread curve upfront parallel shift Tweak Property
//
//   Date        : Nov 2007
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/ParSpreadUpfrontParallelTP.hpp"

DRLIB_BEGIN_NAMESPACE

ParSpreadUpfrontParallelTP::ParSpreadUpfrontParallelTP(): CObject(TYPE) {}
ParSpreadUpfrontParallelTP::~ParSpreadUpfrontParallelTP() {}

static void ParSpreadUpfrontParallelTP_load(CClassSP& clazz) {
    REGISTER(ParSpreadUpfrontParallelTP, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<ParSpreadUpfrontParallelTP>::iObject);
}

CClassConstSP const ParSpreadUpfrontParallelTP::TYPE = CClass::registerClassLoadMethod("ParSpreadUpfrontParallelTP", typeid(ParSpreadUpfrontParallelTP), ParSpreadUpfrontParallelTP_load);

RiskProperty_TYPES(ParSpreadUpfrontParallelTP)

DRLIB_END_NAMESPACE

