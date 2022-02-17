//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VirtualDestructorBase.hpp
//
//   Description : Class which has a virtual destructor and nothing else
//                 Base class to use if you want a common way to free derived
//                 types
//
//   Author      : Mark A Robson
//
//   Date        : 4 Mar 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VIRTUALBASEDESTRUCTOR_HPP
#define EDR_VIRTUALBASEDESTRUCTOR_HPP
#include "edginc/smartPtr.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

//// Base class to allows us to free the derived types using a single type as
//// well as implementing the core method required by smartPtr template.
//// Classes must derive virtually from this class (otherwise you get multiple
//// refCounts!)
class TOOLKIT_DLL VirtualDestructorBase{
    mutable int refCount;
public:
    // when we use copy constructor, we don't want to copy the number of references
    // FIXME: incompatible with "in-place new" ? (in this case we'd like to keep refCount)
    VirtualDestructorBase(const VirtualDestructorBase& /*rhs*/) : refCount(0) {}

    /** Leave refCount on lhs and rhs unchanged, so things like *sp1 = *sp2 works correctly */
    VirtualDestructorBase& operator= (const VirtualDestructorBase& /*rhs*/)
    {
        return *this;
    }

    virtual ~VirtualDestructorBase() {}; // was in Functor.cpp
    /** For smartPtr template */
    int& getRefCount() const { return refCount;} // for speed

    //// Ideally this wouldn't be here but we've started to smartPtr a lot
    //// with this class and VC71 insists (in dll mode) on instantiating the
    //// smartPtr clone method
    template<class T> static T* copy(const T* orig){
        throwError();  // in Functor.cpp - this really throws an exception
        return 0; // avoid compilation issues
    }
protected:
    VirtualDestructorBase() : refCount(0) {};  // was in Functor.cpp
private:
    static void throwError();  // was in Functor.cpp
};
// #if !defined(DEBUG) || defined(QLIB_OBJECT_CPP)
// OPT_INLINE int& VirtualDestructorBase::getRefCount() const{ return refCount;}
// #endif

DECLARE(VirtualDestructorBase);

DRLIB_END_NAMESPACE
#endif
