//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : ReadableMap.hpp
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

#ifndef EDR_READABLE_MAP_HPP
#define EDR_READABLE_MAP_HPP
#include "edginc/Map.hpp"

DRLIB_BEGIN_NAMESPACE
/** Interface which marks objects which have map like properties  */
class TOOLKIT_DLL IReadableMap: public virtual IMap{
public:
    static CClassConstSP const TYPE;
    // retrieve an element from the map
    virtual IObjectSP get(const string& key) const = 0;
    // Use IMap::IIterator to retrieve all the elements
};
typedef smartPtr<IReadableMap> IReadableMapSP;

#ifndef QLIB_READABLEMAP_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<IReadableMap>);
#endif

DRLIB_END_NAMESPACE
#endif
