//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : IMSLFuncNDProvider.hpp
//
//   Description :
//
//   Date        : 4 Sep 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Function.hpp"
#include "edginc/DECLARE.hpp"

#ifndef EDR_IMSLFUNCTIONNDPROVIDER_HPP
#define EDR_IMSLFUNCTIONNDPROVIDER_HPP

DRLIB_BEGIN_NAMESPACE

typedef double (IMSLFunctionND)(int , double*);

/** Class that gives out static function pointers for use with eg IMSL 
    ie use this to turn C++ types with methods into a C stlye function pointer.
    This is not thread-safe but could be made so with some simple locking.
*/
class IMSLFunctionNDProvider{
    int           index;
public:
    IMSLFunctionNDProvider(const FunctionNDDouble* object);

    ~IMSLFunctionNDProvider();

    IMSLFunctionND* function();

private:
    static int getFirstNullEntry();
    
    //// wrapped functions
    static double wrappedMethod1(int n, double* x);
    static double wrappedMethod2(int n, double* x);
    static double wrappedMethod3(int n, double* x);
    static double wrappedMethod4(int n, double* x);

    static IMSLFunctionND* globalFuncs[];
    static int numGlobalObjects;
    static const FunctionNDDouble* globalObjects[];
};

DRLIB_END_NAMESPACE

#endif
