/* =======================================================================================
   FILENAME :       srt_h_futuresfct.h

   AUTHOR :         C. Drozo
   ======================================================================================= */

#ifndef SRT_H_FUTURESFCT_H
#define SRT_H_FUTURESFCT_H

#include "srt_h_futures.h"
/*---------------------------------------------------------*/
Err future_Rutkowski(long index, double* dPeriod, double* dParams, double* futuresRate);

Err future_HJM(long index, double* dPeriod, double* dParams, double* futuresRate);

Err future_HL(long index, double* dPeriod, double* futuresRate);

Err future_Bloom(long index, double* dPeriod, double* dParams, double* futuresRate);

/*---------------------------------------------------------*/
char* eFutFunc_Rut(
    double dIndex, double* dParams, double* dConvexity, double* dDerivative, int nParams);

char* eFutFunc_HJM(
    double dIndex, double* dParams, double* dConvexity, double* dDerivative, int nParams);

char* eFutFunc_Bloom(
    double dIndex, double* dParams, double* dConvexity, double* dDerivative, int nParams);
/*---------------------------------------------------------*/
char* eBuildVols(long lToday);

#endif  // SRT_H_FUTURESFCT_H
