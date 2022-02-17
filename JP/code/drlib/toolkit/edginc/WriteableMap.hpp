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

#ifndef EDR_WRITEABLE_MAP_HPP
#define EDR_WRITEABLE_MAP_HPP
#include "edginc/Map.hpp"

DRLIB_BEGIN_NAMESPACE
/** Interface which marks objects which have map like properties  */
class TOOLKIT_DLL IWriteableMap: public virtual IMap{
public:
    static CClassConstSP const TYPE;
    // add an element to the map
    virtual void put(const string& key, const IObjectSP& obj) = 0;
};
typedef smartPtr<IWriteableMap> IWriteableMapSP;

#ifndef QLIB_WRITEABLE_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<IWriteableMap>);
#endif

DRLIB_END_NAMESPACE
#endif
