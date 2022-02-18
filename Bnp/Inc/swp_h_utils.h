/*--------------------------------------------------------------
        FILE: swp_h_utils.h
        PURPOSE: Useful swap utilities
        AUTHOR: Dimitri Mayevski
        DATE: 3/10/2002
  --------------------------------------------------------------*/

#ifndef __SWP_H_UTILS_H__
#define __SWP_H_UTILS_H__

#include "uterror.h"

enum ESwapElements
{
    SE_NOTIONALS = 1,
    SE_FIXED     = 2,
    SE_SPREADS   = 4
};

Err SwapElements_Init(
    long     start,
    long     end,
    char*    freqStr,
    char*    basisStr,
    char*    refrateStr,
    char*    recpayStr,
    double   fixed,
    int      elements,
    long     today,
    double** times,
    double** cflows,
    int*     len);

typedef struct _SCashFlows
{
    int      nex;
    double*  ex;
    int*     nmat;
    double **mat, **cpn;
} SCashFlows;

Err CashFlows_Init(SCashFlows* g, int nex);
Err CashFlows_Free(SCashFlows* g);

#endif /* #ifndef __SWP_H_UTILS_H__ */