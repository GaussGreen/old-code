
#ifndef Fx3FBetaGrfnH
#define Fx3FBetaGrfnH

#include "Fx3FGrfn.h"
#include "grf_h_mdlcomm.h"

Err grfn_payoff_4_3dfxBeta_tree(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double spot_fx, void *dom_yc, double dom_lam, double dom_phi, void *for_yc,
    double for_lam, double for_phi,
    /* Nodes data */
    long n1, long n2, long n3,
    /* i: d1  , j: d2  , k: d3  , l = {0: xDom  , 1: xFor  , 2: log (Fx/Fx0)} */
    double beta, double ****sv,
    /* Vector of results to be updated */
    long nprod, double ****prod_val,
    /* Barrier details */
    int is_bar, int bar_k, int **bar_idx, int bar_col, double bar_lvl);

Err grfn_payoff_4_3dfx_Betamc(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double spot_fx, void *dom_yc, double dom_lam, double dom_phi, void *for_yc,
    double for_lam, double for_phi, double Xdom, double Yfor, double Zfx,
    /* Results */
    int num_col, double *res, int *stop_path);

#endif
