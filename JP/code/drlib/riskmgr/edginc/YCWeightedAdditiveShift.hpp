//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : YCWeightedAdditiveShift.hpp
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

#ifndef _YCWEIGHTADDITIVESHIFT__HPP
#define _YCWEIGHTADDITIVESHIFT__HPP

#include "edginc/YCWeightedShift.hpp"


DRLIB_BEGIN_NAMESPACE

/** Scenario shift: Weighted additive shift to a Yield Curve.
    Benchmarks can be shifted by different amounts. If the user defined 
    shiftSizes lie between Yield Curve benchmarks then we interpolate. 
    Extrapolate flat off the ends. */
class RISKMGR_DLL YCWeightedAdditiveShift: public YCWeightedShift {
public:    
    static CClassConstSP const TYPE;
 
    /* Additive shift to the DoubleArray passed as an argument */
     virtual void shiftArray(ExpiryArraySP rate_expiries, 
                             DoubleArray* rates,
                             DateTime& today);

protected:
    YCWeightedAdditiveShift(const CClassConstSP& clazz);

private:
    friend class YCWeightedAdditiveShiftHelper;
    // for reflection
    YCWeightedAdditiveShift();
};

DRLIB_END_NAMESPACE

#endif
