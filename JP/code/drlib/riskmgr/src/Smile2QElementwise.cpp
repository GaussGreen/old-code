/**
 * @file Smile2QElementwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_Smile2QElementwise_CPP
#include "edginc/BoxedInt.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/Smile2QElementwise.hpp"

DRLIB_BEGIN_NAMESPACE

Smile2QElementwise::Smile2QElementwise(): CObject(TYPE) {}
Smile2QElementwise::~Smile2QElementwise() {}

static void Smile2QElementwise_load(CClassSP& clazz) {
    REGISTER(Smile2QElementwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<Smile2QElementwise>::iObject);
}

CClassConstSP const Smile2QElementwise::TYPE = CClass::registerClassLoadMethod("Smile2QElementwise", typeid(Smile2QElementwise), Smile2QElementwise_load);

RiskProperty_TYPES(Smile2QElementwise)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<Smile2QElementwise>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<Smile2QElementwise>);

/** Included in RiskMgrLib::linkInClasses() to force linkage into the
    Windows exe.  */
bool Smile2QElementwiseLinkIn() {
    return IRestorableWithRespectTo<Smile2QElementwise>::TYPE != NULL;
}

DRLIB_END_NAMESPACE
