//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FunctorMember.hpp
//
//   Description : Templated class which wraps C++ pointers to member functions.
//                 The requirement being that these member functions take no
//                 no parameters (and return a templated type)
//
//   Author      : Mark A Robson
//
//   Date        : 4 Mar 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FUNCTORMEMBER_HPP
#define EDR_FUNCTORMEMBER_HPP

#include "edginc/Functor.hpp"

DRLIB_BEGIN_NAMESPACE

/** Functor for member functions. Here T is the class that contains the member
    function. ReturnType is what our Functor returns, whilst WrappedReturnType
    is what the member function returns */
template<class T, class ReturnType, class WrappedReturnType = ReturnType> 
class FunctorMember: public Functor<ReturnType> {
public:
    typedef WrappedReturnType (T::* Func)(void);//C++ pointer to member function
    // constructor
    FunctorMember(Func func): func(func) {}
    ~FunctorMember(){} // and destructor
    // Implementation of execute
    virtual ReturnType execute(IObject* obj) {
        // Cast to correct type and call member function on object
        T* t = DYNAMIC_CAST(T, obj); 
        return (t->*func)();
    }
private:
    Func func;      // Pointer to member function
};

DRLIB_END_NAMESPACE
#endif
