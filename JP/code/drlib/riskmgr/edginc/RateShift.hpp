//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RateShift.hpp
//
//   Description : Scenario shift: Yield Curve benchmarks can be shifted by
//                 different amounts. If the user defined shiftSizes lie between
//                 Yield Curve benchmarks then we interpolate. Extrapolate flat off
//                 the ends.
//                 
//                 ** NOTE **
//                 ==========
//                 This class is to be retired and should not be used.
//                 Use class YCWeightedAdditiveShift instead.
//
//   Author      : Stephen Hope
//
//   Date        : 7 May 2003
//
//
//----------------------------------------------------------------------------

#ifndef _RATESHIFT__HPP
#define _RATESHIFT__HPP

#include "edginc/YCWeightedAdditiveShift.hpp"

DRLIB_BEGIN_NAMESPACE

/** Scenario shift: Yield Curve benchmarks can be shifted by
    different amounts. If the user defined shiftSizes lie between
    Yield Curve benchmarks then we interpolate. Extrapolate flat off
    the ends. */

class RISKMGR_DLL RateShift: public YCWeightedAdditiveShift {
public:    
    static CClassConstSP const TYPE;

private:
    friend class RateShiftHelper;
    // for reflection
    RateShift();
};

DRLIB_END_NAMESPACE

#endif
