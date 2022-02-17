#ifndef Fx5FGrfnH
#define Fx5FGrfnH

#include "grf_h_mdlcomm.h"

Err grfn_payoff_4_5dfx_mc(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double spot_fx, double R1D, double R2D, double R1F, double R2F, double Z,
    /* Results */
    int num_col, double *res, int *stop_path);

#endif
