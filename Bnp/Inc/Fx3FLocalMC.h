
#ifndef Fx3FLocalMCH
#define Fx3FLocalMCH

#include "Fx3FQuadGRFN.h"

/*	Main function */
/*	------------- */
Err mcLocal_main_3dfx(
    /*	Time data */
    long    npaths,
    long    nsteps,
    int     num_col,
    double* time,
    double* date,
    double* dom_ifr,
    double* dom_std,
    double* dom_phi,
    double* for_ifr,
    double* for_std,
    double* for_phi,
    double* fx_fwd,
    double* fx_std,

    LOCAL_MODEL_FUNC*  model_func,
    LOCAL_MODEL_PARAM* model_param,

    /*	Product data */
    void** func_parm_tab,
    int*   has_evt,
    /*	Model data */
    double  dom_lam,
    double  for_lam,
    double* corr_dom_for,
    double* corr_dom_fx,
    double* corr_for_fx,
    /*	Market data */
    double spot_fx,
    char*  dom_yc,
    char*  for_yc,
    /* do PECS adjustment */
    int do_pecs,
    /*	Payoff function */
    Err (*payoff_func)(
        /* Event */
        double evt_date,
        double evt_time,
        void*  func_parm,
        /* Market data */
        double spot_fx,
        void*  dom_yc,
        double dom_lam,
        double dom_phi,
        void*  for_yc,
        double for_lam,
        double for_phi,
        double Xdom,
        double Yfor,
        double Zfx,
        /* Results */
        int     num_col,
        double* res,
        int*    stop_path),
    /*	Results */
    double** res);

/*	Main function */
/*	------------- */
Err mcLocal_main_3dfx2(
    /*	Time data */
    long    npaths,
    long    nsteps,
    int     num_col,
    double* time,
    double* date,
    double* dom_ifr,
    double* dom_std,
    double* dom_phi,
    double* for_ifr,
    double* for_std,
    double* for_phi,
    double* fx_fwd,
    double* fx_std,
    double (*vol_ln_func)(double t, double Spot, double Fwd, double* params),
    double* params,
    /*	Product data */
    void** func_parm_tab,
    int*   has_evt,
    /*	Model data */
    double  dom_lam,
    double  for_lam,
    double* corr_dom_for,
    double* corr_dom_fx,
    double* corr_for_fx,
    /*	Market data */
    double spot_fx,
    char*  dom_yc,
    char*  for_yc,
    /* do PECS adjustment */
    int do_pecs,
    /*	Payoff function */
    Err (*payoff_func)(
        /* Event */
        double evt_date,
        double evt_time,
        void*  func_parm,
        /* Market data */
        double spot_fx,
        void*  dom_yc,
        double dom_lam,
        double dom_phi,
        void*  for_yc,
        double for_lam,
        double for_phi,
        double Xdom,
        double Yfor,
        double Zfx,
        /* Results */
        int     num_col,
        double* res,
        int*    stop_path),
    /*	Results */
    double** res);

#endif
