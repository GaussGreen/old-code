//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 05-Oct-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_FUNCTIONOPERATIONS_HPP
#define QLIB_FUNCTIONOPERATIONS_HPP

#include "edginc/Function.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Collection of utility (static) operations on MFunctionND and derived classes
 * */
class UTIL_DLL FunctionOperations
{
public:
    /** Build g(x) = f1(x) * f2(x) defined on the input ranges */
    static MFunctionNDSP multiply(const MFunctionND& f1, const MFunctionND& f2, const RangeArray& ranges);

    /** Build g(x) = f1(x) + f2(x) defined on the input ranges */
    static MFunctionNDSP add(const MFunctionND& f1, const MFunctionND& f2, const RangeArray& ranges);

    // Build g(x) = f(x+shift)
    static MFunctionNDSP shift(const CDoubleArray& shift, const MFunctionND& f);

    /**
     * If f:R^N->R^M such that f(x) = (f_0(x), f_1(x), ..., f_M-1(x)),
     * builds g(x) = f_space(x)
     * 
     * NB: "space" belongs to [0, ..., M-1]
     * */
    static MFunctionNDSP project(const MFunctionND& f, int space);

    /**
     * If f:R^N->R^M with f(x) = f(x_0, x_1, ..., x_N-1),
     * builds g:R^(N-1)->R^M such that
     * g(x_1, ...,x_space-1, x_space+1, ...,x_N-1) = f(x_1, ...,x_space-1, y, x_space+1, ...,x_N-1)
     * 
     * NB: "space" belongs to [0, ..., M-1]
     * */
    static MFunctionNDSP partialFunction(const MFunctionND& f, int space, double y);
};

DRLIB_END_NAMESPACE

#endif /*QLIB_FUNCTIONOPERATIONS_HPP*/
