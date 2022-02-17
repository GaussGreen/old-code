#ifndef FUNCTIONWRAPPER_HPP
#define FUNCTIONWRAPPER_HPP

DRLIB_BEGIN_NAMESPACE

#include "edginc/ModelException.hpp"

/** Policy class for non-null pointers */
class PtrNullCheck {
public:
    /** Checks a pointer for null */
    template<class T> static void check(const T* ptr) {
        static const string routine = "PtrChecker::check";
        if(!ptr) {
            throw ModelException(routine, "Null pointer detected");
        }
    }
};


/** No check */
class UTIL_DLL PtrNoCheck {
public:
    /** No check - will be inlined */
    template<class T> static void check(const T* ptr) { }
};


/** Wrapper around double obj->f(double x) where f() is a member function 
    of an object. Allows parameterization w.r.t. Object pointer (raw vs smart), 
    member function signature (const vs non const) and checking policy (null check
    vs no check). Based on Alexandrescu "Modern C++ design" p.117 - 118.
    May want to generalize the implemented operator() to arbitrary
    member function signatures and implement parameter binding. */
template<class _ObjPtr, class _MemFuncPtr, class CheckPolicy = PtrNullCheck>
class MemFuncWrapper {
public:
    /** Full constructor */
    MemFuncWrapper(const _ObjPtr& ptr, _MemFuncPtr func):
    ptr(ptr), func(func) { }

    /** Destructor */
    ~MemFuncWrapper() {}

    /** Call member function */
    double operator()(double x) const {
        // Check for null object
        CheckPolicy::check(ptr);
        
        double eval = ((*ptr).*func)(x);
        return eval;
    }

private:
    _ObjPtr      ptr;      //!< Pointer to object
    _MemFuncPtr  func;     //!< Pointer to member function
};

/** Wrapper around f(x) + c where c is a shift constant. The main application
    is numerical solution of equations. Solvers assume solution of the form f(x*) = 0 
    while author's may not want to pollute the code of f(x) by subtracting a target 
    value. */
template<class _Func>
class FuncShiftWrapper {
public:
    /** Full constructor */
    FuncShiftWrapper(_Func& func, double shift = 0.0): func(func), shift(shift) { }

    /** Destructor */
    ~FuncShiftWrapper() {}
    
    /** Resets constant offset */
    void resetShift(double shift) {
        this->shift = shift;
    }

    /** Call member function and shift */
    double operator()(double x) const {
        double eval = func(x) + shift;
        return eval;
    }

private:
    _Func  func;        //!< Underlying functor object
    double shift;       //!< Scalar shift
};

DRLIB_END_NAMESPACE

#endif
