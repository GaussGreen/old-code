//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMTrancheUtils.cpp
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CCMTrancheUtils.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Calculate tranche expected loss, which is given by 
 * E(min(max(L(0,t)+L(t,T) - K1, 0), K2-K1))
 *
 * note:
 * for wrappers, L(0,t) is the historical loss from obsStart to obsEnd
 * density is the loss distribution of L(t,T), which includes all the names
 * that were not marked as defaulted.
 */
double CCMTrancheUtils::calcExpectedTrancheLoss(
    double             k1,       /* (I) lower strike in $ */
    double             k2,       /* (I) upper strike in $ */
    double             histLoss, /* (I) historical loss L(0,t) in $ */
    double             lossUnit, /* (I) loss unit in $ */
    const DoubleArray& density,  /* (I) L(t,T) loss density */
    int                maxS,     /* (I) max index on the Short side */
    int                maxL){    /* (I) max index on the Long side */
    if(k2<k1)
    {
        return 0.;
    }
    
    double tmp = 0.;

    /* For performance reasons (we can call ccmCollar over 200,000 times in ONE 
     * typical testcase) consider two cases here: if k1>0 we can call collar() 
     * directly instead of ccmCollar */
    if (k1 >= 0 ) {
    
#ifdef PRINT_LOSS_DISTRIBUTION  
        FILE *fp = fopen("C:\\debug2.txt","a");
#endif
        for(int i= -maxS; i <= maxL; ++i) {
#ifdef PRINT_LOSS_DISTRIBUTION
            fprintf(fp, "%.10lf\t", density[i+maxS]);
#endif
            tmp += Maths::shiftedCollar(histLoss + i*lossUnit, k1, k2) * density[i+maxS];
        }
#ifdef PRINT_LOSS_DISTRIBUTION
        fprintf(fp, "\n");
        fclose(fp);
#endif
    } else {
        for(int i= -maxS; i <= maxL; ++i) {
            tmp += Maths::creditCollar(histLoss + i*lossUnit, k1, k2) * density[i+maxS];
        }
    }
    return tmp;
}

DRLIB_END_NAMESPACE

