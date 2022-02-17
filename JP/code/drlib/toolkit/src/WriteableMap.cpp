//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : WriteableMap.hpp
//
//   Description : Interface which marks objects which have map like 
//                 properties including an ability to put elements directly
//                 into the map
//
//   Author      : Mark A Robson
//
//   Date        : 13 Jan 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_WRITEABLE_CPP
#include "edginc/WriteableMap.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<IWriteableMap>);

static void loadWriteableMap(CClassSP& clazz){
    REGISTER_INTERFACE(IWriteableMap, clazz);
    EXTENDS(IMap);
}

CClassConstSP const IWriteableMap::TYPE = CClass::registerInterfaceLoadMethod(
    "IWriteableMap", typeid(IWriteableMap), loadWriteableMap);

DRLIB_END_NAMESPACE
