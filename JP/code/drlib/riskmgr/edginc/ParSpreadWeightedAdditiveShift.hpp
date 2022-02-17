//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParSpreadWeightedAdditiveShift.hpp
//
//   Description : Scenario shift: Par Spreads benchmarks can be shifted by
//                 different amounts. If the user defined shiftSizes lie between
//                 ParSpread benchmarks then we interpolate. Extrapolate flat off
//                 the ends.
//
//   Author      : Jose Hilera
//
//   Date        : 15-June-2005
//
//
//----------------------------------------------------------------------------

#ifndef _PARSPREADWEIGHTADDITIVESHIFT__HPP
#define _PARSPREADWEIGHTADDITIVESHIFT__HPP

#include "edginc/ParSpreadWeightedShift.hpp"


DRLIB_BEGIN_NAMESPACE

/** Scenario shift: Weighted additive shift to a Par Spread Curve.
    Benchmarks can be shifted by different amounts. If the user defined 
    shiftSizes lie between Par Spread benchmarks then we interpolate. 
    Extrapolate flat off the ends. */
class RISKMGR_DLL ParSpreadWeightedAdditiveShift: public ParSpreadWeightedShift {
public:    
    static CClassConstSP const TYPE;
 
    /* Additive shift to the DoubleArray passed as an argument */
     virtual void shiftArray(ExpiryArraySP rate_expiries, 
                             DoubleArray* rates,
                             DateTime& today);

private:
    friend class ParSpreadWeightedAdditiveShiftHelper;
    // for reflection
    ParSpreadWeightedAdditiveShift();
};

DRLIB_END_NAMESPACE

#endif
