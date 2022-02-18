#ifndef __3F_QUANTO_TREE_H
#define __3F_QUANTO_TREE_H

Err tree_main_3fquanto(
    /*	Time data */
    long    nstp,
    double* time,
    double* date,
    int*    vol_change, /*	1 if one of the volatilities has changed, 0 otherwise */

    /*	Term structures & Distributions */
    double  r1_lam,
    double* r1_sig,
    double* r1_fwd,
    double* r1_var,

    double  r2_lam,
    double* r2_sig,
    double* r2_fwd,
    double* r2_var,

    double* r1_r2_covar,

    double* lgm2F_ifr,
    char*   lgm2F_yc,
    int     lgm2F_is_dom_for,

    double  r3_lam,
    double* r3_sig,
    double* r3_fwd,
    double* r3_var,

    double* lgm1F_ifr,
    char*   lgm1F_yc,

    double* fx_sig,

    /*	Product data */
    void** func_parm_tab,
    int*   eval_evt,

    /*	Correlations */
    double* corr_r1_r2,
    double* corr_r1_r3,
    double* corr_r1_fx,
    double* corr_r2_r3,
    double* corr_r2_fx,
    double* corr_r3_fx,

    /*	Payoff function */
    Err (*payoff_func)(
        /* Event */
        double evt_date,
        double evt_time,
        void*  func_parm,
        /* Market data */
        void*  lgm2F_yc,
        double r1_lambda,
        double r2_lambda,
        double r1_phi,
        double r2_phi,
        double r1_r2_phi,
        double rho,
        void*  lgm1F_yc,
        double r3_lambda,
        double r3_phi,
        int    lgm2F_is_dom_for,
        /* Nodes data */
        long n1,
        long n2,
        long n3,
        /* i: d1, j: d2, k: d3, l = {0: r1, 1: r2, 2: r3} */
        double**** sv,
        /* Vector of results to be updated */
        long       nprod,
        double**** prod_val),
    /*	Result */
    int     nprod,
    int     discount,
    double* res);

#endif
