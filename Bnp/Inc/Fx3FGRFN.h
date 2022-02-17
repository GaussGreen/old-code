#ifndef Fx3FGrfnH
#define Fx3FGrfnH

#include "grf_h_mdlcomm.h"

typedef struct {
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  int fx_idx;
  int dom_idx; /* -1 if none */
  int for_idx; /* -1 if none */

  long num_fx_df;
  double *fx_df_tms;
  long *fx_df_dts;

  long num_dom_df;
  double *dom_df_tms;
  long *dom_df_dts;

  long num_for_df;
  double *for_df_tms;
  long *for_df_dts;

} grfn_parm_tree, *GRFNPARMTREE;

Err grfn_payoff_4_3dfx_tree(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double spot_fx, void *dom_yc, double dom_lam, double dom_phi, void *for_yc,
    double for_lam, double for_phi,
    /* Nodes data */
    long n1, long n2, long n3,
    /* i: d1  , j: d2  , k: d3  , l = {0: xDom  , 1: xFor  , 2: log (Fx/Fx0)} */
    double ****sv,
    /* Vector of results to be updated */
    long nprod, double ****prod_val,
    /* Barrier details */
    int is_bar, int bar_k, int **bar_idx, int bar_col, double bar_lvl);

typedef struct {
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  int fx_idx;
  int dom_idx; /* -1 if none */
  int for_idx; /* -1 if none */

  long num_fx_df;
  double *fx_df_tms;
  long *fx_df_dts;
  double *fx_dff;
  double *fx_gam;
  double *fx_gam2;
  int do_fx;

  long num_dom_df;
  double *dom_df_tms;
  long *dom_df_dts;
  double *dom_dff;
  double *dom_gam;
  double *dom_gam2;
  int do_dom;

  long num_for_df;
  double *for_df_tms;
  long *for_df_dts;
  double *for_dff;
  double *for_gam;
  double *for_gam2;
  int do_for;

} grfn_parm_mc, *GRFNPARMMC;

Err grfn_payoff_4_3dfx_mc(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double spot_fx, void *dom_yc, double dom_lam, double dom_phi, void *for_yc,
    double for_lam, double for_phi, double Xdom, double Yfor, double Zfx,
    /* Results */
    int num_col, double *res, int *stop_path);

#endif
