//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Functor.hpp
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

#ifndef EDR_FUNCTOR_HPP
#define EDR_FUNCTOR_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

//// Functor base class for each return type 
template<class ReturnType> class Functor: public virtual VirtualDestructorBase{
public:
    virtual ReturnType execute(IObject* obj) = 0;
    virtual ~Functor();
    Functor();
};

/** specialisation of Functor for return types actually used - make life
    easier for linker */
template <> class TOOLKIT_DLL Functor<IObjectSP>: public virtual VirtualDestructorBase{
public:
    virtual IObjectSP execute(IObject* obj) = 0;
    virtual ~Functor();
    Functor();
};
template <> class TOOLKIT_DLL Functor<string>: public virtual VirtualDestructorBase{
public:
    virtual string execute(IObject* obj) = 0;
    virtual ~Functor();
    Functor();
};
template <> class TOOLKIT_DLL Functor<double>: public virtual VirtualDestructorBase{
public:
    virtual double execute(IObject* obj) = 0;
    virtual ~Functor();
    Functor();
};

template <> class TOOLKIT_DLL Functor<int>: public virtual VirtualDestructorBase{
public:
    virtual int execute(IObject* obj) = 0;
    virtual ~Functor();
    Functor();
};

template <> class TOOLKIT_DLL Functor<bool>: public virtual VirtualDestructorBase{
public:
    virtual bool execute(IObject* obj) = 0;
    virtual ~Functor();
    Functor();
};

DRLIB_END_NAMESPACE
#endif
