
#ifndef Fx3FQuadTreeH
#define Fx3FQuadTreeH

#include "Fx3FQuadGrfn.h"

/*	Main function */
/*	------------- */

Err treeQuad_main_3dfx(
    /*	Time data */
    long nstp, double *time, double *date,
    int *vol_change, /*	1 if one of the volatlities has changed      ,
                                         0 otherwise */
    double *sig_dom, /*	Term structures */
    double *sig_for, double *sig_fx,

    double *dom_ifr, /*	Distributions */
    double *dom_fwd, double *dom_var, double *for_ifr, double *for_fwd,
    double *for_var, double *fx_fwd_approx, double *fx_fwd0,
    double *fx_var_approx,
    /*	Product data */
    void **func_parm_tab, int *eval_evt, double *bar_lvl, int *bar_col,
    int *is_bar,
    /*	Model data */
    double dom_lam, double for_lam, double corr_dom_for, double corr_dom_fx,
    double corr_for_fx,
    /*	Market data */
    double spot_fx, char *dom_yc, char *for_yc,

    /*	Transformation functions */
    LOCAL_MODEL_FUNC *model_func, LOCAL_MODEL_PARAM *model_param,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, long step,
                       void *func_parm,
                       /* Market data */
                       double fwd_fx, void *dom_yc, double dom_lam,
                       double dom_phi, void *for_yc, double for_lam,
                       double for_phi,
                       /* Nodes data */
                       long n1, long n2, long n3,
                       /* i: d1      , j: d2      , k: d3      ,
                                       l = {0: xDom      , 1: xFor      , 2: log
                          (Fx/Fx0)} */
                       double ****sv,
                       /* function to go get back fx */
                       LOCAL_MODEL_FUNC *model_func,
                       LOCAL_MODEL_PARAM *model_param,

                       /* Vector of results to be updated */
                       long nprod, double ****prod_val,
                       /* Barrier details */
                       int is_bar, int bar_k, int **bar_idx, int bar_col,
                       double bar_lvl),
    /*	Result */
    int nprod, int discount, double *res);

#endif
