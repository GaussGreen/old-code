
#ifndef __DIAGCALIBDLMSV_H
#define __DIAGCALIBDLMSV_H

#include "CPDCalib.h"

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
Err LGMSVcalibration_dlm(

    int fix_lambda, /*	0: calib lambda to cap      , 1: fix lambda calib
                                        to diagonal */
    int nlam,       /*	Lambda TS: may NOT be changed in the process */
    double lam_time[], double lam[],

    void *AllPrimaryInst, void *AllSecondaryInst, void *Model, void *Params,

    Err (*CalibrationFunction)(void *AllPrimaryInst, void *Model, void *Params),

    Err (*PricingFunction)(void *AllSecondaryInst, void *Model, void *Params,
                           int iIndex, double *dPrice),

    Err (*BumpingModelFunction)(void *Model, int iIndex, double NewValue),

    int nb_exe_dates, double *target_prices, double *target_vegas,

    DIAG_CALIB_LM_PARAMS lm_params);

#endif