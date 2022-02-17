#ifndef Fx3FBetaDLMGrfnH
#define Fx3FBetaDLMGrfnH

#include "grf_h_mdlcomm.h"

typedef struct {
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  int dom_idx; /* -1 if none */
  int for_idx; /* -1 if none */
  int fx_idx;  /* -1 if none */

  /* Bond Payment reconstruction */
  double bond_pay_const;
  double bond_pay_lin;

  /* Spot Fx reconstruction */
  double fx_coef_const;
  double fx_coef_lin;
  double fx_coef_quad;
  double fx_coef_dom;
  double fx_coef_for;

  /* Fx discount factors */
  int do_fx;
  long num_fx_df;
  double *fx_df_tms;
  long *fx_df_dts;
  double *df_fx_coef_const;
  double *df_fx_coef_lin;

  /* Domestic discount factors */
  int do_dom;
  long num_dom_df;
  double *dom_df_tms;
  long *dom_df_dts;
  double *df_dom_coef_const;
  double *df_dom_coef_lin;

  /* Foreign discount factors */
  int do_for;
  long num_for_df;
  double *for_df_tms;
  long *for_df_dts;
  double *df_for_coef_const;
  double *df_for_coef_lin;

} grfn_parm_beta_MC_DLM, *GRFNPARMBETAMCDLM;

typedef struct {
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  int dom_idx; /* -1 if none */
  int for_idx; /* -1 if none */
  int fx_idx;  /* -1 if none */

  /* Preclalculations */
  long today;

  double lam_dom;
  double lam_for;
  double phi_dom;
  double phi_for;
  double beta_dom;
  double beta_for;

  double B0;
  double C0;
  double varFFX;

  /* Fx discount factors */
  long num_fx_df;
  double *fx_df_tms;
  long *fx_df_dts;

  /* Domestic discount factors */
  long num_dom_df;
  double *dom_df_tms;
  long *dom_df_dts;

  /* Foreign discount factors */
  long num_for_df;
  double *for_df_tms;
  long *for_df_dts;

} grfn_parm_beta_Tree_DLM, *GRFNPARMBETATREEDLM;

Err grfn_payoff_4_3dfxBetaDLM_mc(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Model data */
    double Xdom, double Yfor, double Zfx,
    /* Results */
    int num_col, double *res, int *stop_path);

Err grfn_payoff_4_3dfxBetaDLMQTStar_Tree(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Nodes data */
    long n1, long n2, long n3, double ****sv,
    /* Results */
    long nprod, double ****prod_val,
    /* Barrier details */
    int is_bar, int bar_k, int **bar_idx, int bar_col, double bar_lvl);

Err grfn_payoff_4_3dfxBetaDLMQBeta_Tree(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market Data */
    double spot_fx, char *dom_yc, char *for_yc,
    /* Nodes data */
    long n1, long n2, long n3, double ****sv,
    /* Results */
    long nprod, double ****prod_val,
    /* Barrier details */
    int is_bar, int bar_k, int **bar_idx, int bar_col, double bar_lvl);

#endif
