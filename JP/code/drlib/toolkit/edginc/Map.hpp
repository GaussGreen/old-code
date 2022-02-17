//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : Map.hpp
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

#ifndef EDR_MAP_HPP
#define EDR_MAP_HPP
#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE
/** Interface which marks objects which have map like properties  */
class TOOLKIT_DLL IMap: public virtual IObject{
public:
    static CClassConstSP const TYPE;

    /** Iterator. Derived from IObject so easy to pass through DRInterface */
    class TOOLKIT_DLL IIterator: public virtual IObject{
    public:
        static CClassConstSP const TYPE;
        //// are there any elements left to iterate over
        virtual bool hasMoreElements() const = 0;
        //// get the key for the current element
        virtual const string& getKey() const = 0;
        //// get the current element (returns a reference)
        virtual IObjectSP getElement() const = 0;
        //// increment the iterator
        virtual void increment() = 0;
    };
    typedef smartPtr<IIterator> IIteratorSP;

    //// Builds an iterator
    virtual IIteratorSP createIterator() = 0;
    //// Is this object truly a map ie does toObject()/toMap() return this
    virtual bool isTrueMap() const = 0;
};
typedef smartPtr<IMap> IMapSP;

#ifndef QLIB_HASHTABLE_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<IMap::IIterator>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<IMap>);
#endif

DRLIB_END_NAMESPACE
#endif
