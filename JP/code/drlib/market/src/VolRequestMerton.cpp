//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequestMerton.cpp
//
//   Description : Merton vol request 
//
//   Author      : Oliver Brockhaus
//
//   Date        : April 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestMerton.hpp"

DRLIB_BEGIN_NAMESPACE

VolRequestMerton::VolRequestMerton(): 
    CVolRequest(TYPE){}

class VolRequestMertonHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VolRequestMerton, clazz);
        SUPERCLASS(CVolRequest);
        EMPTY_SHELL_METHOD(defaultVolRequestMerton);
    }

    static IObject* defaultVolRequestMerton(){
        return new VolRequestMerton();
    }
};

CClassConstSP const VolRequestMerton::TYPE = 
CClass::registerClassLoadMethod("VolRequestMerton", typeid(VolRequestMerton), 
                                VolRequestMertonHelper::load);


DRLIB_END_NAMESPACE

