/**
 * @file LocalCorrExpiryAndStrikewise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_LOCALCORREXPIRYANDSTRIKEWISE_CPP
#include "edginc/Expiry.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/LocalCorrExpiryAndStrikewise.hpp"

DRLIB_BEGIN_NAMESPACE

LocalCorrExpiryAndStrikewise::LocalCorrExpiryAndStrikewise(): CObject(TYPE) {}
LocalCorrExpiryAndStrikewise::~LocalCorrExpiryAndStrikewise() {}

static void LocalCorrExpiryAndStrikewise_load(CClassSP& clazz) {
    REGISTER(LocalCorrExpiryAndStrikewise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<LocalCorrExpiryAndStrikewise>::iObject);
}

CClassConstSP const LocalCorrExpiryAndStrikewise::TYPE = 
    CClass::registerClassLoadMethod("LocalCorrExpiryAndStrikewise", 
                                    typeid(LocalCorrExpiryAndStrikewise), 
                                    LocalCorrExpiryAndStrikewise_load);

RiskProperty_TYPES(LocalCorrExpiryAndStrikewise)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<LocalCorrExpiryAndStrikewise>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<LocalCorrExpiryAndStrikewise>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL RiskProperty<LocalCorrExpiryAndStrikewise>);
DRLIB_END_NAMESPACE
