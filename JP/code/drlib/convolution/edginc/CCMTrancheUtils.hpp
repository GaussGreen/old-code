//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMTrancheUtils.hpp
//
//   Date        : Aug 2004
//
//----------------------------------------------------------------------------
#ifndef EDR_CCMTRANCHEUTILS_HPP
#define EDR_CCMTRANCHEUTILS_HPP

#include "edginc/AtomicArray.hpp"
DRLIB_BEGIN_NAMESPACE

/** Class containing a couple of little utility methods */
class CONVOLUTION_DLL CCMTrancheUtils{
public:
    /**
     * Calculate tranche expected loss, which is given by 
     * E(min(max(L(0,t)+L(t,T) - K1, 0), K2-K1))
     *
     * note:
     * for wrappers, L(0,t) is the historical loss from obsStart to obsEnd
     * density is the loss distribution of L(t,T), which includes all the names
     * that were not marked as defaulted.
     */
    static double calcExpectedTrancheLoss(
        double             k1,       /* (I) lower strike in $ */
        double             k2,       /* (I) upper strike in $ */
        double             histLoss, /* (I) historical loss L(0,t) in $ */
        double             lossUnit, /* (I) loss unit in $ */
        const DoubleArray& density,  /* (I) L(t,T) loss density */
        int                maxS,     /* (I) max index on the Short side */
        int                maxL);    /* (I) max index on the Long side */
};

DRLIB_END_NAMESPACE

#endif
