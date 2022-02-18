
#ifndef Fx3FMCH
#define Fx3FMCH

#include "MCEBOptimisation.h"

/*	Main function */
/*	------------- */
Err mc_main_3dfx(
    /*	Time data */
    long    npaths,
    int     num_col,
    double* time,
    double* date,
    long    nb_dates,
    double* dom_ifr, /*	Distributions */
    double* dom_fwd,
    double* dom_std,
    double* dom_phi,
    double* dom_beta,
    double* dom_bond_pay,
    double* dom_beta_pay,
    double* for_ifr,
    double* for_fwd,
    double* for_std,
    double* for_phi,
    double* for_beta,
    double* fx_fwd,
    double* fx_std,
    double* dom_for_cov,
    double* dom_fx_cov,
    double* for_fx_cov,
    /*	Product data */
    void** func_parm_tab,
    /*	Model data */
    double dom_lam,
    double for_lam,
    double corr_dom_for,
    double corr_dom_fx,
    double corr_for_fx,
    /*	Market data */
    double spot_fx,
    char*  dom_yc,
    char*  for_yc,
    /*  Do PECS adjustment */
    int do_pecs,
    /* for Optimisation of exercise boundary */
    int        do_optimisation,
    int*       optimise,
    MCEBPARAMS params,
    /*	Initialisation function to be called at the beggining of each path
                    or NULL if none */
    void (*init_func)(),
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
    double** res,
    double** optim_res);

#endif
