//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : Map.cpp
//
//   Description : Interface which marks objects which have map like 
//                 properties
//
//   Author      : Mark A Robson
//
//   Date        : 13 Jan 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_HASHTABLE_CPP
#include "edginc/Map.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<IMap::IIterator>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<IMap>);

static void loadIMap(CClassSP& clazz){
    REGISTER_INTERFACE(IMap, clazz);
    EXTENDS(IObject);
}

static void loadIMapIterator(CClassSP& clazz){
    REGISTER_INTERFACE(IMap::IIterator, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IMap::TYPE = CClass::registerInterfaceLoadMethod(
    "IMap", typeid(IMap), loadIMap);

CClassConstSP const IMap::IIterator::TYPE = CClass::registerInterfaceLoadMethod(
    "IMap::IIterator", typeid(IMap::IIterator), loadIMapIterator);

DRLIB_END_NAMESPACE
