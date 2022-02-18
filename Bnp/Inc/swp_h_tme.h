#include <NUM_H_ALLHDR.H>

#include "utallhdr.h"

#ifndef SWP_H_CMS_H
#define SWP_H_CMS_H

/*--------------------------------------------------------------*/

Err swp_f_tme(
    long    TMEDate,
    double  dTECSpread,
    char*   szSwapYieldCurveName,
    char*   szBondYieldCurveName,
    char*   szBondVolCurveName,
    char*   cSwapRefRateCode,
    char*   cBondRefRateCode,
    int     iNumStrikesInVol,
    double* Strikes,
    int     bAdjForSpread,
    Err (*GetVol)(
        double      dStart,
        double      dEnd,
        double      dStrike,
        SRT_Boolean bAdjForSpread,
        double      dForward,
        double      dSpread,
        double*     pdBsVol),
    double* dAnswer);

#endif