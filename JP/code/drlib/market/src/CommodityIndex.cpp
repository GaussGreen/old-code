//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CommodityIndex.cpp
//
//   Description : Commodity index asset
//
//   Author      : Andrew J Swain
//
//   Date        : 30 September 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CommodityIndex.hpp"

DRLIB_BEGIN_NAMESPACE

CommodityIndex::CommodityIndex(): ContangoCommodity(TYPE) {}

class CommodityIndexHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        clazz->setDescription("Commodity index asset");
        REGISTER(CommodityIndex, clazz);
        SUPERCLASS(ContangoCommodity);
        EMPTY_SHELL_METHOD(defaultCommodityIndex);
    }

    static IObject* defaultCommodityIndex(){
        return new CommodityIndex();
    }
};

CClassConstSP const CommodityIndex::TYPE = CClass::registerClassLoadMethod(
    "CommodityIndex", typeid(CommodityIndex), CommodityIndexHelper::load);

DRLIB_END_NAMESPACE

