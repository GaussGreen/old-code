//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMLossUnit.hpp
//
//   Description : Calculates loss unit
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------
#ifndef EDR_CCMLOSSUNIT_HPP
#define EDR_CCMLOSSUNIT_HPP

#include "edginc/AtomicArray.hpp"
DRLIB_BEGIN_NAMESPACE

class CONVOLUTION_DLL CCMLossUnit{
public:
    /**
     * Calculate loss unit (bin size for loss distribution discretization)
     * lossAmt and lossAmtCata must already be allocated with size n
     * if a name is defaulted, it is still in the lossAmt list with amt 0.
     */
    static double calcLossUnit(
        const DoubleArray& notionals,    /* (I) name notionals */
        int                maxSlice,     /* (I) max nb of slice              */
        int                subdivision); /* (I) subdivion of the GCD */
    
    
    /* Calcualtes bounds for the density arrays */ 
    static void calcLossDensityBounds(
        const DoubleArray& amt,           /* (I) amounts as a double, 
                                             usually loss Cata */
        int&               maxLongLoss,   /* (O) sum(long loss amt) */
        int&               maxShortLoss); /* (O) sum(short loss amount) */
};

DRLIB_END_NAMESPACE
#endif
