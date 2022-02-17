//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GetMarket.cpp
//
//   Description : What an class must implement to get market data
//                 from the cache
//
//   Author      : Andrew J Swain
//
//   Date        : 24 September 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

IGetMarket::IGetMarket() {}

IGetMarket::~IGetMarket() {}

void IGetMarket::load(CClassSP& clazz){
    REGISTER_INTERFACE(IGetMarket, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IGetMarket::TYPE = CClass::registerInterfaceLoadMethod(
    "GetMarket", typeid(IGetMarket), load);

DRLIB_END_NAMESPACE
