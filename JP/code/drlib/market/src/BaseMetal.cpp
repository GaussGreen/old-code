//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BaseMetal.cpp
//
//   Description : Base metal asset
//
//   Author      : Andrew McCleery
//
//   Date        : 14 July 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BaseMetal.hpp"

DRLIB_BEGIN_NAMESPACE

BaseMetal::BaseMetal(): ContangoCommodity(TYPE) {}

class BaseMetalHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        clazz->setDescription("Base metal asset");
        REGISTER(BaseMetal, clazz);
        SUPERCLASS(ContangoCommodity);
        EMPTY_SHELL_METHOD(defaultBaseMetal);
    }

    static IObject* defaultBaseMetal(){
        return new BaseMetal();
    }
};

CClassConstSP const BaseMetal::TYPE = CClass::registerClassLoadMethod(
    "BaseMetal", typeid(BaseMetal), BaseMetalHelper::load);

DRLIB_END_NAMESPACE

