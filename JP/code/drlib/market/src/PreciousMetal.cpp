//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PreciousMetal.cpp
//
//   Description : Commodity asset for gold, silver etc.
//
//   Author      : Andrew J Swain
//
//   Date        : 30 September 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PreciousMetal.hpp"

DRLIB_BEGIN_NAMESPACE

PreciousMetal::PreciousMetal(): ContangoCommodity(TYPE) {}

class PreciousMetalHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        clazz->setDescription("PreciousMetal asset");
        REGISTER(PreciousMetal, clazz);
        SUPERCLASS(ContangoCommodity);
        EMPTY_SHELL_METHOD(defaultPreciousMetal);
    }

    static IObject* defaultPreciousMetal(){
        return new PreciousMetal();
    }
};

CClassConstSP const PreciousMetal::TYPE = CClass::registerClassLoadMethod(
    "PreciousMetal", typeid(PreciousMetal), PreciousMetalHelper::load);

DRLIB_END_NAMESPACE

