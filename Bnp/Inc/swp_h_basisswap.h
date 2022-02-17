

#include "utallhdr.h"
#include <NUM_H_ALLHDR.H>

/* ========================================================================== */

#ifndef BASISSWAP_H
#define BASISSWAP_H

Err CPBasisSwap(
    long StartDate, double Maturity,
    char *cshortRefRateCode, /* Rate corresponding to the maturity of the CP */
    char *clongRefRateCode,  /* Rate corresponding to the Libor leg */
    char *szYieldCurveName,
    double *
        ShortVol, /* Vol used to take into account the delayed payment effect */
    double *ConstSpread, double *SpikeLibor, double *SpikeCP,
    double *Spread, /* Spread CP-Libor for each open day */
    double *Dates,  /* Open days during the life of the swap */
    double *TotalMargins,
    double
        *MarginsSpreads /* Part of the margin coming from the spread CP-Libor */
);

Err srt_f_CPCap(long StartDate, long EndDate, double Strike, double chi,
                char *cCPRefRateCode, /* Rate corresponding to the CP */
                char *cRefRateCode,   /* Rate corresponding to the Libor 3M */
                char *szYieldCurveName, char *szVolCurveName, double *price,
                SrtCallPutType CallPut);

#endif