#ifndef ESL_MARKET_DOT_H
#define ESL_MARKET_DOT_H

#ifdef  __cplusplus
extern "C" {
#endif

#include "esl_types.h"

int EslMktVolCalibration(
    MKTVOL_DATA*        mktvol_data,  // (O) Volatility data
    char const*         Index,        // (I) Index to calibrate
    T_CURVE const*      t_curve,      // (I) Term structure data
    BASEVOL_DATA const* baseVol_data, // (I) Base vol
    SWAPVOL_DATA const* swapVol_data);// (I) Swaption matrix

/* read basevol.dat and swapvol.dat */
int EslReadVolsW(
    BASEVOL_DATA*   baseVol,
    SWAPVOL_DATA*   swapVol,
    char const*     baseVolFilename,
    char const*     swapVolFilename);

#ifdef  __cplusplus
}
#endif

#endif



