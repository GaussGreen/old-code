/**
 * @file LocalCorrExpiry.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_LOCALCORREXPIRY_CPP
#include "edginc/Expiry.hpp"
#include "edginc/BoxedEnum.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/LocalCorrExpiry.hpp"

DRLIB_BEGIN_NAMESPACE

LocalCorrExpiry::LocalCorrExpiry(Operation op): CObject(TYPE), op(op) {}

LocalCorrExpiryConstSP LocalCorrExpiry::SP(Operation op) {
    return LocalCorrExpiryConstSP(new LocalCorrExpiry(op));
}

LocalCorrExpiry::~LocalCorrExpiry() {}

static void LocalCorrExpiry_load(CClassSP& clazz) {
    REGISTER(LocalCorrExpiry, clazz);
    SUPERCLASS(CObject);
    FIELD(op, "Either shift or skew sensitivity");
    EMPTY_SHELL_METHOD(DefaultConstructor<LocalCorrExpiry>::iObject);
}

CClassConstSP const LocalCorrExpiry::TYPE = 
    CClass::registerClassLoadMethod("LocalCorrExpiry", 
                                    typeid(LocalCorrExpiry), 
                                    LocalCorrExpiry_load);

RiskProperty_TYPES(LocalCorrExpiry)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<LocalCorrExpiry>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<LocalCorrExpiry>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<LocalCorrExpiry>);

START_PUBLIC_ENUM_DEFINITION(LocalCorrExpiry::Operation, "How correlations are to be tweaked");
ENUM_VALUE_AND_NAME(LocalCorrExpiry::PARALLEL_SHIFT, "Parallel", "Parallel shift up or down");
ENUM_VALUE_AND_NAME(LocalCorrExpiry::SKEWED_SHIFT, "Skewed", "Skewed shift up or down");
END_ENUM_DEFINITION(LocalCorrExpiry::Operation);

DRLIB_END_NAMESPACE
