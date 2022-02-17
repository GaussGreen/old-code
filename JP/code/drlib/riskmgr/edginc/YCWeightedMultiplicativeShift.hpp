//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : YCWeightedMultiplicativeShift.hpp
//
//   Description : Scenario shift: Yield Curve benchmarks can be shifted by
//                 different amounts. If the user defined shiftSizes lie between
//                 Yield Curve benchmarks then we interpolate. Extrapolate flat off
//                 the ends.
//
//   Author      : Jose Hilera
//
//   Date        : 15-June-2005
//
//
//----------------------------------------------------------------------------

#ifndef _YCWEIGHTMULTIPLICATIVESHIFT__HPP
#define _YCWEIGHTMULTIPLICATIVESHIFT__HPP

#include "edginc/YCWeightedShift.hpp"


DRLIB_BEGIN_NAMESPACE

/** Scenario shift: Weighted multiplicative shift to a Yield Curve.
    Benchmarks can be shifted by different amounts. If the user defined 
    shiftSizes lie between Yield Curve benchmarks then we interpolate. 
    Extrapolate flat off the ends. */
class RISKMGR_DLL YCWeightedMultiplicativeShift: public YCWeightedShift {
public:    
    static CClassConstSP const TYPE;
 
    /* Additive shift to the DoubleArray passed as an argument */
     virtual void shiftArray(ExpiryArraySP rate_expiries, 
                             DoubleArray* rates,
                             DateTime& today);

private:
    friend class YCWeightedMultiplicativeShiftHelper;
    // for reflection
    YCWeightedMultiplicativeShift();
};

DRLIB_END_NAMESPACE

#endif
