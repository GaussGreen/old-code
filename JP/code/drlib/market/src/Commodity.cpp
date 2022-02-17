//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Commodity.cpp
//
//   Description : Commodity asset base class 
//
//   Author      : Andrew J Swain
//
//   Date        : 30 September 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Commodity.hpp"

DRLIB_BEGIN_NAMESPACE

Commodity::Commodity(CClassConstSP clazz): CAsset(clazz) {}

/** for reflection */
Commodity::Commodity(): CAsset(TYPE) {}


class CommodityHelper {
public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(Commodity, clazz);
        SUPERCLASS(Asset);
    }
};

CClassConstSP const Commodity::TYPE = CClass::registerClassLoadMethod(
    "Commodity", typeid(Commodity), CommodityHelper::load);

DRLIB_END_NAMESPACE
