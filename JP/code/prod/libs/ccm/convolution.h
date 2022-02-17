#include "payoff_tp.h"

/** calculate loss density */
int calcLossDensity(
    const int    *lossAmt,   /* (I) credit portfolio (loss amt)       */
    int     n,         /* (I) number of names                   */
    int     maxLoss,   /* (I) sum(loss amt)                     */
    const COPULA *c,         /* (I) copula 0=gauss, 1=skew            */
    const double *p,         /* (I) array[n] of uncond survival proba */
    double *density    /* (O) array[maxLoss+1] of loss density  */
    );

