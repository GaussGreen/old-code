//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParSpreadWeightedMultiplicativeShift.hpp
//
//   Description : Scenario shift: Par Spread Curve benchmarks can be shifted by
//                 different amounts. If the user defined shiftSizes lie between
//                 Par Spread benchmarks then we interpolate. Extrapolate flat 
//                 off the ends.
//
//   Author      : Jose Hilera
//
//   Date        : 15-June-2005
//
//
//----------------------------------------------------------------------------

#ifndef _PARSPREADWEIGHTMULTIPLICATIVESHIFT__HPP
#define _PARSPREADWEIGHTMULTIPLICATIVESHIFT__HPP

#include "edginc/ParSpreadWeightedShift.hpp"


DRLIB_BEGIN_NAMESPACE

/** Scenario shift: Weighted multiplicative shift to a Par Spread Curve.
    Benchmarks can be shifted by different amounts. If the user defined 
    shiftSizes lie between Par Spread benchmarks then we interpolate. 
    Extrapolate flat off the ends. */
class RISKMGR_DLL ParSpreadWeightedMultiplicativeShift: public ParSpreadWeightedShift {
public:    
    static CClassConstSP const TYPE;
    const static string NAME;
 
    /* Additive shift to the DoubleArray passed as an argument */
     virtual void shiftArray(ExpiryArraySP rate_expiries, 
                             DoubleArray* rates,
                             DateTime& today);

private:
    friend class ParSpreadWeightedMultiplicativeShiftHelper;
    // for reflection
    ParSpreadWeightedMultiplicativeShift();
};

DRLIB_END_NAMESPACE

#endif
