//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : AssetSwapStrike.cpp
//
//   Description : Base class for asset swap strike objects 
//
//   Author      : André Segger
//
//   Date        : 02 June 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AssetSwapStrike.hpp"

DRLIB_BEGIN_NAMESPACE


class AssetSwapStrikeHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AssetSwapStrike, clazz);
        SUPERCLASS(CObject);
        // no fields
    }
};

CClassConstSP const AssetSwapStrike::TYPE = CClass::registerClassLoadMethod(
    "AssetSwapStrike", typeid(AssetSwapStrike), AssetSwapStrikeHelper::load);


DRLIB_END_NAMESPACE
