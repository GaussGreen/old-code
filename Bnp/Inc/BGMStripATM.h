

#include <num_h_allhdr.h>

#include "utallhdr.h"

/* ========================================================================== */

#ifndef BGMStripATM_H
#define BGMStripATM_H

/* Coded in BGMStripATM.c */
Err srt_f_BGMSensitivityVolMatrix(
    long     MaxNumPeriod,
    long     MaturityInPeriod,
    long     nCorrelRow,
    double** ppdHistCorrel,
    char*    szRefRate,
    char*    szYieldCurveName,
    char*    szVolCurveName,
    char*    szCorrelMode,
    long     iStart, /* of the swaption to bump */
    long     iEnd,   /* of the swaption to bump */
    double   bump,
    double** answer);

/* Coded in BGMStripATM.c */
Err srt_f_BGMOrthogonalSensitivity(
    long     MaxNumPeriod,
    long     MaturityInPeriod,
    long     nCorrelRow,
    double** ppdHistCorrel,
    char*    szRefRate,
    char*    szYieldCurveName,
    char*    szVolCurveName,
    char*    szCorrelMode,
    long     nInstrLiquid,
    long*    iStartLiquid,
    long*    iEndLiquid,
    long     IndextoBump,
    double   bump,
    double** answer);

/* Coded in BGMStripATM.c */
Err srt_f_DeformationsBGMSABR(
    int       MaxNumPeriod,
    int       MaturityInPeriod,
    double**  Correl,
    char*     szRefRate,
    char*     szYieldCurveName,
    char*     szVolCurveName,
    char*     szCorrelMode,
    long*     IMatCube,
    int       NumMatCube,
    long*     IUndCube,
    int       NumUndCube,
    long*     IStartLiquidATM,
    long*     IEndLiquidATM,
    int       NumLiquidATM,
    long*     IStartLiquidSABR,
    long*     IEndLiquidSABR,
    int       NumLiquidSABR,
    double**  TSLibor,
    double*** ATMSens,
    double*** alphaSens,
    double*** rhoSens,
    double    bumpATM,
    double    bumpalpha,
    double    bumprho);

/* Coded in BGMStripATM.c */
Err srt_f_CalibrateATMBGMSABR(
    char*     szRefRate,
    char*     szYieldCurveName,
    char*     szVolCurveName,
    char*     szCorrelMode,
    int       MaxNumPeriod,
    int       MaturityInPeriod,
    int       NumLiquidATM,
    int       NumLiquidSABR,
    int       NumMatCube,
    int       NumUndCube,
    long*     IMatCube,
    long*     IUndCube,
    long*     IStartLiquidATM,
    long*     IEndLiquidATM,
    double**  Correl,
    double*** ATMSens,
    double    bumpATM,
    double    bumpalpha,
    double    bumprho,
    double**  OutputATM,
    double**  OutputAlpha,
    double**  OutputRho);

/* Coded in BGMStripATM.c */
Err srt_f_getinfofromund(
    char*      dom_nme,
    int*       MaxNumPeriod,
    int*       MaturityInPeriod,
    int*       NumLiquidATM,
    int*       NumLiquidSABR,
    int*       NumMatCube,
    int*       NumUndCube,
    long**     IMatCube,
    long**     IUndCube,
    long**     IStartLiquidATM,
    long**     IEndLiquidATM,
    double***  Correl,
    double**** ATMSens,
    double*    bumpATM,
    double*    bumpalpha,
    double*    bumprho);

#endif
