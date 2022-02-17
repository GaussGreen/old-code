//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Functor.cpp
//
//   Description : Templated class which have an execute method taking an
//                 IObject and return the templated type. Useful for wrapping
//                 C++ pointers to member functions.
//
//   Author      : Mark A Robson
//
//   Date        : 4 Mar 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Functor.hpp"

DRLIB_BEGIN_NAMESPACE
//VirtualDestructorBase::VirtualDestructorBase(): refCount(0){}
//VirtualDestructorBase::~VirtualDestructorBase(){}
void VirtualDestructorBase::throwError(){
    throw ModelException("VirtualDestructorBase::throwError",
                         "Clone is not supported on VirtualDestructorBase");
}

/** specialisation of Functor for return types actually used - make life
    easier for linker */
Functor<IObjectSP>::~Functor(){}
Functor<IObjectSP>::Functor(){}
Functor<string>::~Functor(){}
Functor<string>::Functor(){}
Functor<double>::~Functor(){}
Functor<double>::Functor(){}
Functor<int>::~Functor(){}
Functor<int>::Functor(){}
Functor<bool>::~Functor(){}
Functor<bool>::Functor(){}

DRLIB_END_NAMESPACE
