//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : IMSLFunctionNDProvider.cpp
//
//   Description : 
//
//   Date        : 4 Sep 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/IMSLFunctionNDProvider.hpp"
#include <math.h>

DRLIB_BEGIN_NAMESPACE

/** Class that gives out static function pointers for use with eg IMSL 
    ie use this to turn C++ types with methods into a C stlye function pointer.
    This is not thread-safe but could be made so with some simple locking.
*/
IMSLFunctionNDProvider::IMSLFunctionNDProvider(const FunctionNDDouble* object):
        index(getFirstNullEntry())
{
        globalObjects[index] = object; // set global static variable
}

IMSLFunctionNDProvider::~IMSLFunctionNDProvider()
{
        // clear globalObject
        globalObjects[index] = 0;
}

IMSLFunctionND* IMSLFunctionNDProvider::function()
{
        return globalFuncs[index];
}


int IMSLFunctionNDProvider::getFirstNullEntry()
{
        // search for first null entry
        for (int i = 0; i < numGlobalObjects; i++){
            if (!globalObjects[i]){
                return i;
            }
        }
        throw ModelException(
            "IMSLFunctionNDProvider (in IntFuncNDGeneral.cpp)",
            "Run out of spare static variables");
}
    
//// wrapped functions
double IMSLFunctionNDProvider::wrappedMethod1(int n, double* x)
{
        CDoubleArray xx = CDoubleArray(n);
        for (int i=0; i<n; ++i)
        {
            xx[i] = x[i];
        }
        return (*globalObjects[0])(xx);
}

double IMSLFunctionNDProvider::wrappedMethod2(int n, double* x)
{
        CDoubleArray xx = CDoubleArray(n);
        for (int i=0; i<n; ++i)
        {
            xx[i] = x[i];
        }
        return (*globalObjects[1])(xx);
}

double IMSLFunctionNDProvider::wrappedMethod3(int n, double* x)
{
        CDoubleArray xx = CDoubleArray(n);
        for (int i=0; i<n; ++i)
        {
            xx[i] = x[i];
        }
        return (*globalObjects[2])(xx);
}

double IMSLFunctionNDProvider::wrappedMethod4(int n, double* x)
{
        CDoubleArray xx = CDoubleArray(n);
        for (int i=0; i<n; ++i)
        {
            xx[i] = x[i];
        }
        return (*globalObjects[3])(xx);
}

//// set up statics
IMSLFunctionND* IMSLFunctionNDProvider::globalFuncs[] = 
    {
        wrappedMethod1,
        wrappedMethod2,
        wrappedMethod3,
        wrappedMethod4
    };
int IMSLFunctionNDProvider::numGlobalObjects =
    sizeof(globalFuncs)/sizeof(IMSLFunctionND*);
const FunctionNDDouble* IMSLFunctionNDProvider::globalObjects[] = {0,0,0,0};

DRLIB_END_NAMESPACE

