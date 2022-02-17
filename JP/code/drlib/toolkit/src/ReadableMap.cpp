//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : ReadableMap.cpp
//
//   Description : Interface which marks objects which have map like 
//                 properties including an ability to get elements directly
//                 out of the map
//
//   Author      : Mark A Robson
//
//   Date        : 13 Jan 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_READABLEMAP_CPP
#include "edginc/ReadableMap.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<IReadableMap>);

static void loadReadableMap(CClassSP& clazz){
    REGISTER_INTERFACE(IReadableMap, clazz);
    EXTENDS(IMap);
}

CClassConstSP const IReadableMap::TYPE = CClass::registerInterfaceLoadMethod(
    "IReadableMap", typeid(IReadableMap), loadReadableMap);

DRLIB_END_NAMESPACE
