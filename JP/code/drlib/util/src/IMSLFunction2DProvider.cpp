//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : IMSLFunction2DProvider.cpp
//
//   Description : 
//
//   Date        : 4 Sep 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/IMSLFunction2DProvider.hpp"
#include <math.h>

DRLIB_BEGIN_NAMESPACE

/** Class that gives out static function pointers for use with eg IMSL 
    ie use this to turn C++ types with methods into a C stlye function pointer.
    This is not thread-safe but could be made so with some simple locking.
*/
IMSLFunction2DProvider::IMSLFunction2DProvider(const FunctionNDDouble* object):
        index(getFirstNullEntry())
{
        globalObjects[index] = object; // set global static variable
}

IMSLFunction2DProvider::~IMSLFunction2DProvider()
{
        // clear globalObject
        globalObjects[index] = 0;
}

IMSLFunction2D* IMSLFunction2DProvider::function()
{
        return globalFuncs[index];
}

IMSLBoundFunction *IMSLFunction2DProvider::lowerBound()
{
        return globalFuncsLowerB[index];
}

IMSLBoundFunction *IMSLFunction2DProvider::upperBound()
{
        return globalFuncsUpperB[index];
}

int IMSLFunction2DProvider::getFirstNullEntry()
{
        // search for first null entry
        for (int i = 0; i < numGlobalObjects; i++){
            if (!globalObjects[i]){
                return i;
            }
        }
        throw ModelException(
            "IMSLFunction2DProvider (in IntFuncNDGeneral.cpp)",
            "Run out of spare static variables");
}
    
//// wrapped functions
double IMSLFunction2DProvider::wrappedMethod1(double x, double y)
{
        CDoubleArray xx = CDoubleArray(2);
        xx[0] = x;
        xx[1] = y;
        return (*globalObjects[0])(xx);
}

double IMSLFunction2DProvider::wrappedMethod2(double x, double y)
{
        CDoubleArray xx = CDoubleArray(2);
        xx[0] = x;
        xx[1] = y;
        return (*globalObjects[1])(xx);
}

double IMSLFunction2DProvider::wrappedMethod3(double x, double y)
{

        CDoubleArray xx = CDoubleArray(2);
        xx[0] = x;
        xx[1] = y;
        return (*globalObjects[2])(xx);
}

double IMSLFunction2DProvider::wrappedMethod4(double x, double y)
{
        CDoubleArray xx = CDoubleArray(2);
        xx[0] = x;
        xx[1] = y;
        return (*globalObjects[3])(xx);
}

double IMSLFunction2DProvider::wrappedMethodLowerB1(double x)
{
    static const string method = "IMSLFunction2DProvider::wrappedMethodLowerB1";
    if ((*globalObjects[0]).getNbVars() != 2)
    {
        throw ModelException(method,
            "internal error: function provided is not a 2d function");            
    }
    const RangeArray& interval = (*globalObjects[0]).getIntervals();
    Range::checkIsNonEmpty(*(interval[1]));
    const Boundary& b = (*(interval[1])).getLower();
    if (b.isInfinite())
    {
        throw ModelException(method,
                    "infinite intervals are not supported; got " + (*(interval[1])).toString());        
    }
    return b.getValue();
}

double IMSLFunction2DProvider::wrappedMethodLowerB2(double x)
{
    static const string method = "IMSLFunction2DProvider::wrappedMethodLowerB2";
    if ((*globalObjects[1]).getNbVars() != 2)
    {
        throw ModelException(method,
            "internal error: function provided is not a 2d function");            
    }
    const RangeArray& interval = (*globalObjects[1]).getIntervals();
    Range::checkIsNonEmpty(*(interval[1]));
    const Boundary& b = (*(interval[1])).getLower();
    
     if ( !b.isClosedBracket() || b.isInfinite())
     {
         throw ModelException(method,
                             "(Semi-)Open and / or infinite intervals are not supported; got " + (*(interval[1])).toString());
     }
     return b.getValue();
}

double IMSLFunction2DProvider::wrappedMethodLowerB3(double x)
{
    static const string method = "IMSLFunction2DProvider::wrappedMethodLowerB3";
    if ((*globalObjects[2]).getNbVars() != 2)
    {
        throw ModelException(method,
            "internal error: function provided is not a 2d function");            
    }
    const RangeArray& interval = (*globalObjects[2]).getIntervals();
    Range::checkIsNonEmpty(*(interval[1]));
    const Boundary& b = (*(interval[1])).getLower();
    if (b.isInfinite())
    {
        throw ModelException(method,
                    "infinite intervals are not supported; got " + (*(interval[1])).toString());        
    }
    return b.getValue();
}

double IMSLFunction2DProvider::wrappedMethodLowerB4(double x)
{
    static const string method = "IMSLFunction2DProvider::wrappedMethodLowerB4";
    if ((*globalObjects[3]).getNbVars() != 2)
    {
        throw ModelException(method,
            "internal error: function provided is not a 2d function");            
    }
    const RangeArray& interval = (*globalObjects[3]).getIntervals();
    Range::checkIsNonEmpty(*(interval[1]));
    const Boundary& b = (*(interval[1])).getLower();
    if (b.isInfinite())
    {
        throw ModelException(method,
                    "infinite intervals are not supported; got " + (*(interval[1])).toString());        
    }
    return b.getValue();
}

double IMSLFunction2DProvider::wrappedMethodUpperB1(double x)
{
    static const string method = "IMSLFunction2DProvider::wrappedMethodUpperB1";
    if ((*globalObjects[0]).getNbVars() != 2)
    {
        throw ModelException(method,
            "internal error: function provided is not a 2d function");            
    }
    const RangeArray& interval = (*globalObjects[0]).getIntervals();
    Range::checkIsNonEmpty(*(interval[1]));
    const Boundary& b = (*(interval[1])).getUpper();
    if (b.isInfinite())
    {
        throw ModelException(method,
                    "infinite intervals are not supported; got " + (*(interval[1])).toString());        
    }
    return b.getValue();
}

double IMSLFunction2DProvider::wrappedMethodUpperB2(double x)
{
    static const string method = "IMSLFunction2DProvider::wrappedMethodUpperB2";
    if ((*globalObjects[1]).getNbVars() != 2)
    {
        throw ModelException(method,
            "internal error: function provided is not a 2d function");            
    }
    const RangeArray& interval = (*globalObjects[1]).getIntervals();
    Range::checkIsNonEmpty(*(interval[1]));
    const Boundary& b = (*(interval[1])).getUpper();
    if (b.isInfinite())
    {
        throw ModelException(method,
                    "infinite intervals are not supported; got " + (*(interval[1])).toString());
    }
    return b.getValue();
}

double IMSLFunction2DProvider::wrappedMethodUpperB3(double x)
{
    static const string method = "IMSLFunction2DProvider::wrappedMethodUpperB3";
    if ((*globalObjects[2]).getNbVars() != 2)
    {
        throw ModelException(method,
            "internal error: function provided is not a 2d function");            
    }
    const RangeArray& interval = (*globalObjects[2]).getIntervals();
    Range::checkIsNonEmpty(*(interval[1]));
    const Boundary& b = (*(interval[1])).getUpper();
    if (b.isInfinite())
    {
        throw ModelException(method,
                    "infinite intervals are not supported; got " + (*(interval[1])).toString());        
    }
    return b.getValue();
}

double IMSLFunction2DProvider::wrappedMethodUpperB4(double x)
{
    static const string method = "IMSLFunction2DProvider::wrappedMethodUpperB4";
    if ((*globalObjects[3]).getNbVars() != 2)
    {
        throw ModelException(method,
            "internal error: function provided is not a 2d function");            
    }
    const RangeArray& interval = (*globalObjects[3]).getIntervals();
    Range::checkIsNonEmpty(*(interval[1]));
    const Boundary& b = (*(interval[1])).getUpper();
    if (b.isInfinite())
    {
        throw ModelException(method,
                    "infinite intervals are not supported; got " + (*(interval[1])).toString());        
    }
    return b.getValue();
}

//// set up statics
IMSLFunction2D* IMSLFunction2DProvider::globalFuncs[] = 
{
    wrappedMethod1,
    wrappedMethod2,
    wrappedMethod3,
    wrappedMethod4
};

IMSLBoundFunction* IMSLFunction2DProvider::globalFuncsLowerB[] =
{
    wrappedMethodLowerB1,
    wrappedMethodLowerB2,
    wrappedMethodLowerB3,
    wrappedMethodLowerB4
};

IMSLBoundFunction* IMSLFunction2DProvider::globalFuncsUpperB[] =
{
    wrappedMethodUpperB1,
    wrappedMethodUpperB2,
    wrappedMethodUpperB3,
    wrappedMethodUpperB4
};

int IMSLFunction2DProvider::numGlobalObjects =
    sizeof(globalFuncs)/sizeof(IMSLFunction2D*);
const FunctionNDDouble* IMSLFunction2DProvider::globalObjects[] = {0,0,0,0};

DRLIB_END_NAMESPACE

