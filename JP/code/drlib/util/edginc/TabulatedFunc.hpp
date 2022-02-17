//---------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TabulatedFunc.hpp
//
//   Description : 
//
//   Date        : May 2002
//
//
//----------------------------------------------------------------------------

#ifndef TABULATEDFUNC_HPP
#define TABULATEDFUNC_HPP

#include <string>
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

/** 2 arrays of doubles {x} and {f(x)}
    Given x, return f(x) with perhaps interpolation */

class UTIL_DLL TabulatedFunc : public CObject {
public:
    static CClassConstSP const TYPE;

    /** how you can interpolate */
    static const string INTERP_NONE;
    static const string INTERP_LINEAR;
    static const string INTERP_STAIRS;
    
    TabulatedFunc(const DoubleArray& xArray,
                  const DoubleArray& fArray,
                  const string&      interp);

    virtual ~TabulatedFunc();

    /** how long is the TabulatedFunc ? */
    int length() const;

    /** v. basic interpolator */
    double interpolate(double x) const;

    /** validation */
    virtual void validatePop2Object();

    /** return interp type */
    string getInterp() const{return interp;}

    /** return date list (deep copy) */
    DoubleArray getOrdinates() const{return xArray;}

    /** return value list (deep copy) */
    DoubleArray getValues() const{return fArray;} 

private:
    friend class TabulatedFuncHelper;

    TabulatedFunc();
    TabulatedFunc(const TabulatedFunc &rhs);
    TabulatedFunc& operator=(const TabulatedFunc& rhs);

    typedef enum InterpType {
        interpNone = 0,
        interpLinear,
        interpStairs
    } InterpType;

    mutable DoubleArray   xArray; // because locate() requires non-const param!
    DoubleArray   fArray;
    string        interp;
    InterpType    interpType; // $unregistered
    
};

typedef smartConstPtr<TabulatedFunc> TabulatedFuncConstSP;
typedef smartPtr<TabulatedFunc> TabulatedFuncSP;
#ifndef QLIB_TABULATEDFUNC_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<TabulatedFunc>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<TabulatedFunc>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<TabulatedFunc>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<TabulatedFunc>);
#endif

DRLIB_END_NAMESPACE
#endif


