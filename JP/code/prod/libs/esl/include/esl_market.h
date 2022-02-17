#ifndef ESL_MARKET_DOT_H
#define ESL_MARKET_DOT_H

#include "esl_types.h"

#ifdef  __cplusplus
extern "C" {
#endif

int EslMktVolCalibration(
    MKTVOL_DATA*        mktvol_data,  /* (O) Volatility data */
    char const*         Index,        /* (I) Index to calibrate */
    T_CURVE const*      t_curve,      /* (I) Term structure data */
    BASEVOL_DATA const* baseVol_data, /* (I) Base vol */
    SWAPVOL_DATA const* swapVol_data);/* (I) Swaption matrix */

int EslMktVolCalibrationAndRecording(
    MKTVOL_DATA*        mktvol_data,  /* (O) Volatility data */
    BASEVOL_EXPOSURE_DATA* selected_baseVols, /* (O) Selected base vol dates */
    SWAPVOL_EXPOSURE_DATA* selected_swapVols, /* (O) Selected swap vol dates */
    char const*         Index,        /* (I) Index to calibrate */
    T_CURVE const*      t_curve,      /* (I) Term structure data */
    BASEVOL_DATA const* baseVol_data, /* (I) Base vol */
    SWAPVOL_DATA const* swapVol_data);/* (I) Swaption matrix */

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



