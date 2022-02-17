
#ifndef Fx3FBetaDLMTreeH
#define Fx3FBetaDLMTreeH

Err tree_main_3dBetaDLM_QBeta(
    /*	Time data */
    long nstp, double *time, double *date,

    /*	Model data */
    int *vol_change,

    double dom_lam, double for_lam,

    double *sig_dom, double *sig_for, double *sig_fx, double *corr_dom_for,
    double *corr_dom_fx, double *corr_for_fx,

    double *dr_const_dom, double *dr_coef_dom, double *dr_const_for,
    double *dr_coef1_for, double *dr_coef2_for, double *dr_coef3_for,
    double *dr_const_fx, double *dr_coef1_fx, double *dr_coef2_fx,
    double *dr_coef3_fx,

    /*	Distributions and Constants */
    double *dom_fwd, double *dom_var, double *for_fwd, double *for_var,
    double *fx_fwd, double *fx_var,

    double *glob_corr_dom_for, double *glob_corr_dom_fx,
    double *glob_corr_for_fx,

    /*	Market data */
    double spot_fx, char *dom_yc, char *for_yc, double *dom_ifr,

    /*	Product data */
    void **func_parm_tab, int *eval_evt, double *bar_lvl, int *bar_col,
    int *is_bar,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,
                       /* Market Data */
                       double spot_fx, char *dom_yc, char *for_yc,
                       /* Nodes data */
                       long n1, long n2, long n3, double ****sv,
                       /* Results */
                       long nprod, double ****prod_val,
                       /* Barrier details */
                       int is_bar, int bar_k, int **bar_idx, int bar_col,
                       double bar_lvl),
    /*	Result */
    int nprod, int discount, double *res);

/*	Main function */
/*	------------- */
Err tree_main_3dBetaDLM_QTStar(
    /*	Time data */
    long nstp, double *time, double *date,
    /*	Model data */
    int *vol_change, /*	1 if one of the volatlities has changed  ,
                                             0 otherwise */
    double *sig_dom, /*	Term structures */
    double *mu_quanto_const, double *mu_quanto_lin, double *sig_for,
    double *sig_fx, double *dom_fwd, double *dom_var, double *for_fwd,
    double *for_var, double *fx_fwd, double *fx_var, double dom_lam,
    double for_lam, double *corr_dom_for, double *corr_dom_fx,
    double *corr_for_fx, double *glob_corr_dom_for, double *glob_corr_dom_fx,
    double *glob_corr_for_fx,

    /*	Product data */
    void **func_parm_tab, int *eval_evt, double *bar_lvl, int *bar_col,
    int *is_bar,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,
                       /* Nodes data */
                       long n1, long n2, long n3, double ****sv,
                       /* Vector of results to be updated */
                       long nprod, double ****prod_val,
                       /* Barrier details */
                       int is_bar, int bar_k, int **bar_idx, int bar_col,
                       double bar_lvl),
    /*	Result */
    int nprod, double *res);

#endif
