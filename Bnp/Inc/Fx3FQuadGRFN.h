
#ifndef Fx3FQuadGrfnH
#define Fx3FQuadGrfnH

#include "Fx3FGrfn.h"
#include "grf_h_mdlcomm.h"

typedef struct {
  long nbSteps;
  double *time;

  double *Zmin;
  double *Zmax;
  double *Umin;
  double *Umax;

  /* Floor volatility function */
  double *volD_param1;
  double *volderivD_param1;

  double *ZU_D_param1;
  double *ZU_D_param2;

  /* Middle volatility function */
  double volM_param1;
  double volM_param2;
  double volM_param3;
  double volderivM_param1;
  double volderivM_param2;

  double ZU_M_param1;
  double ZU_M_param2;
  double ZU_M_param3;
  double ZU_M_param4;

  /* Cap volatility function */
  double *volU_param1;
  double *volderivU_param1;

  double *ZU_U_param1;
  double *ZU_U_param2;

  /* Custom Volatility parameters */
  long nbPoints;
  double *Z_grid;
  double *U_grid;
  double *Vol_grid;
  double *Vol_deriv_grid;

  double *ZS_param1;
  double ZS_param2;

} LOCAL_MODEL_PARAM;

typedef struct {
  /* function to go from Z to U */
  double (*Z_to_U)(long i, double Z, LOCAL_MODEL_PARAM *model_param);

  double (*U_to_Z)(long i, double U, LOCAL_MODEL_PARAM *model_param);

  double (*Z_to_S)(long i, double Z, double Fwd,
                   LOCAL_MODEL_PARAM *model_param);

  double (*S_to_Z)(long i, double S, double Fwd,
                   LOCAL_MODEL_PARAM *model_param);

  double (*vol_func)(long i, double Z, LOCAL_MODEL_PARAM *model_param);

  double (*vol_func_deriv)(long i, double Z, LOCAL_MODEL_PARAM *model_param);

  double (*vol_func_deriv_t)(long i, double Z, LOCAL_MODEL_PARAM *model_param);

  double (*vol_ln_func)(long i, double Z, LOCAL_MODEL_PARAM *model_param);

} LOCAL_MODEL_FUNC;

Err allocate_local_model(long nbSteps, LOCAL_MODEL_FUNC **model_func,
                         LOCAL_MODEL_PARAM **model_param);

Err fill_local_model(LOCAL_MODEL_FUNC *model_func,
                     LOCAL_MODEL_PARAM *model_param, long nbSteps, double *time,
                     double alpha, double beta, double gamma, double sig0,
                     double *fwd_fx);

Err allocate_local_model_cust(long nbPoints, long nbSteps,
                              LOCAL_MODEL_FUNC **model_func,
                              LOCAL_MODEL_PARAM **model_param);
Err fill_local_model_cust(LOCAL_MODEL_FUNC *model_func,
                          LOCAL_MODEL_PARAM *model_param, long nbPoints,
                          double *time, long nbSteps, double alpha, double beta,
                          double gamma, double sig0, double *fwd_fx);

void free_local_model_param(LOCAL_MODEL_FUNC *model_func,
                            LOCAL_MODEL_PARAM *model_param);

void free_local_model_param_cust(LOCAL_MODEL_FUNC *model_func,
                                 LOCAL_MODEL_PARAM *model_param);

double Z_to_U_deltaPo_Quad(long i, double Z, LOCAL_MODEL_PARAM *model_param);
double U_to_Z_deltaPo_Quad(long i, double U, LOCAL_MODEL_PARAM *model_param);
double Z_to_U_deltaNeg_Quad(long i, double Z, LOCAL_MODEL_PARAM *model_param);
double U_to_Z_deltaNeg_Quad(long i, double U, LOCAL_MODEL_PARAM *model_param);

double Z_to_S_Quad(long i, double Z, double Fwd,
                   LOCAL_MODEL_PARAM *model_param);
double S_to_Z_Quad(long i, double S, double Fwd,
                   LOCAL_MODEL_PARAM *model_param);

double vol_func_Quad(long i, double Z, LOCAL_MODEL_PARAM *model_param);
double vol_func_deriv_Quad(long i, double Z, LOCAL_MODEL_PARAM *model_param);
double vol_func_deriv_t_Quad(long i, double Z, LOCAL_MODEL_PARAM *model_param);
double vol_lnATM_func_Quad(double t, double S, double Fwd, double *param);

double vol_ln_func_Quad_Mc(double std, double S, double F, double *params);

double Z_to_U_Quad_Tanh(long i, double Z, LOCAL_MODEL_PARAM *model_param);
double U_to_Z_Quad_Tanh(long i, double U, LOCAL_MODEL_PARAM *model_param);

double Z_to_S_Quad_Tanh(long i, double Z, double Fwd,
                        LOCAL_MODEL_PARAM *model_param);
double S_to_Z_Quad_Tanh(long i, double S, double Fwd,
                        LOCAL_MODEL_PARAM *model_param);

double vol_func_Quad_Tanh(long i, double Z, LOCAL_MODEL_PARAM *model_param);
double vol_func_deriv_Quad_Tanh(long i, double Z,
                                LOCAL_MODEL_PARAM *model_param);
double vol_func_deriv_t_Quad_Tanh(long i, double Z,
                                  LOCAL_MODEL_PARAM *model_param);

Err grfn_payoff_4_3dfxQuad_tree(
    /* Event */
    double evt_date, double evt_time, long step, void *func_parm,
    /* Market data */
    double fwd_fx, void *dom_yc, double dom_lam, double dom_phi, void *for_yc,
    double for_lam, double for_phi,
    /* Nodes data */
    long n1, long n2, long n3,
    /* i: d1      , j: d2      , k: d3      , l = {0: xDom      , 1: xFor      ,
       2: log (Fx/Fx0)} */
    double ****sv, LOCAL_MODEL_FUNC *model_func, LOCAL_MODEL_PARAM *model_param,
    /* Vector of results to be updated */
    long nprod, double ****prod_val,
    /* Barrier details */
    int is_bar, int bar_k, int **bar_idx, int bar_col, double bar_lvl);

Err grfn_payoff_4_3dfx_Quadmc(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double spot_fx, void *dom_yc, double dom_lam, double dom_phi, void *for_yc,
    double for_lam, double for_phi, double Xdom, double Yfor, double Zfx,
    /* Results */
    int num_col, double *res, int *stop_path);

Err FxCall_payoff_4_3dfxQuad_tree(
    /* Event */
    double evt_date, double evt_time, long step, void *func_parm,
    /* Market data */
    double fwd_fx, void *dom_yc, double dom_lam, double dom_phi, void *for_yc,
    double for_lam, double for_phi,
    /* Nodes data */
    long n1, long n2, long n3,
    /* i: d1      , j: d2      , k: d3      , l = {0: xDom      , 1: xFor      ,
       2: log (Fx/Fx0)} */
    double ****sv, LOCAL_MODEL_FUNC *model_func, LOCAL_MODEL_PARAM *model_param,
    /* Vector of results to be updated */
    long nprod, double ****prod_val,
    /* Barrier details */
    int is_bar, int bar_k, int **bar_idx, int bar_col, double bar_lvl);

#endif
