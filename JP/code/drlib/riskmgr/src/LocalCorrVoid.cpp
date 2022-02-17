/**
 * @file LocalCorrVoid.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_LOCALCORRVOID_CPP
#include "edginc/Void.hpp"
#include "edginc/BoxedEnum.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/LocalCorrVoid.hpp"

DRLIB_BEGIN_NAMESPACE

LocalCorrVoid::LocalCorrVoid(Operation op): CObject(TYPE), op(op) {}

LocalCorrVoidConstSP LocalCorrVoid::SP(Operation op) {
    return LocalCorrVoidConstSP(new LocalCorrVoid(op));
}

LocalCorrVoid::~LocalCorrVoid() {}

static void LocalCorrVoid_load(CClassSP& clazz) {
    REGISTER(LocalCorrVoid, clazz);
    SUPERCLASS(CObject);
    FIELD(op, "Either shift or skew sensitivity");
    EMPTY_SHELL_METHOD(DefaultConstructor<LocalCorrVoid>::iObject);
}

CClassConstSP const LocalCorrVoid::TYPE = 
    CClass::registerClassLoadMethod("LocalCorrVoid", 
                                    typeid(LocalCorrVoid), 
                                    LocalCorrVoid_load);

RiskProperty_TYPES(LocalCorrVoid)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<LocalCorrVoid>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<LocalCorrVoid>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<LocalCorrVoid>);

START_PUBLIC_ENUM_DEFINITION(LocalCorrVoid::Operation, "How correlations are to be tweaked");
ENUM_VALUE_AND_NAME(LocalCorrVoid::PARALLEL_SHIFT, "Parallel", "Parallel shift up or down");
ENUM_VALUE_AND_NAME(LocalCorrVoid::SKEWED_SHIFT, "Skewed", "Skewed shift up or down");
END_ENUM_DEFINITION(LocalCorrVoid::Operation);

DRLIB_END_NAMESPACE
