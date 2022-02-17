//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : IMSLFunc2DProvider.hpp
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

typedef double (IMSLFunction2D)(double, double);
typedef double (IMSLBoundFunction)(double);

/** Class that gives out static function pointers for use with eg IMSL 
    ie use this to turn C++ types with methods into a C stlye function pointer.
    This is not thread-safe but could be made so with some simple locking.
*/
class IMSLFunction2DProvider{
    int           index;
public:
    IMSLFunction2DProvider(const FunctionNDDouble* object);

    ~IMSLFunction2DProvider();

    IMSLFunction2D* function();
    IMSLBoundFunction *lowerBound();
    IMSLBoundFunction *upperBound();

private:
    static int getFirstNullEntry();
    
    //// wrapped functions
    static double wrappedMethod1(double x, double y);
    static double wrappedMethod2(double x, double y);
    static double wrappedMethod3(double x, double y);
    static double wrappedMethod4(double x, double y);

    //// wrapped functions for the lower bound
    static double wrappedMethodLowerB1(double x);
    static double wrappedMethodLowerB2(double x);
    static double wrappedMethodLowerB3(double x);
    static double wrappedMethodLowerB4(double x);

    //// wrapped functions for the upper bound
    static double wrappedMethodUpperB1(double x);
    static double wrappedMethodUpperB2(double x);
    static double wrappedMethodUpperB3(double x);
    static double wrappedMethodUpperB4(double x);

    static IMSLFunction2D*    globalFuncs[];
    static IMSLBoundFunction* globalFuncsLowerB[];
    static IMSLBoundFunction* globalFuncsUpperB[];
    static int numGlobalObjects;
    static const FunctionNDDouble* globalObjects[];
};

DRLIB_END_NAMESPACE

#endif
