#ifndef SABRCALIBRRBT
#define SABRCALIBRRBT

#include "DiagcalibGen.h"
#include "srt_h_all.h"

Err calib_sabr_rr_bt_given_beta(
    double           fwd,
    double           mat,
    double           vol_atm,
    double           strike_down,
    double           vol_down,
    double           strike_up,
    double           vol_up,
    double*          sigma_beta,
    double*          alpha,
    double           fixed_beta,
    double*          rho,
    SrtDiffusionType vol_type,
    int              solve_on_vol,
    double           precision,
    int              nb_iter_max,
    double*          calib_err);

/* All the structures needed */
typedef struct
{
    double vol_atm;
    double price_atm;

    double strike_up;
    double target_vol_up;
    double target_price_up;
    double vol_up;
    double price_up;

    double strike_down;
    double target_vol_down;
    double target_price_down;
    double vol_down;
    double price_down;

    SrtDiffusionType vol_type;

} sabr_rrbt_inst, *SABR_RRBT_INST;

typedef struct
{
    double fwd;
    double mat;
    double vol_atm;
    double sigma_beta;
    double alpha;
    double beta;
    double rho;

    double rhoalpha;

    SrtDiffusionType vol_type;

} sabr_rrbt_model, *SABR_RRBT_MODEL;

typedef struct
{
    int             solve_on_vol;
    int             solve_on_rhoalpha;
    CALIBFUNCTIONS  AllFunctions_ForAlpha;
    CALIBGEN_PARAMS CalibParams_ForAlpha;
    CALIBFUNCTIONS  AllFunctions_ForRho;
    CALIBGEN_PARAMS CalibParams_ForRho;

} sabr_rrbt_params, *SABR_RRBT_PARAMS;

Err transform_sabr_beta(
    double           mat,
    double           fwd,
    double           init_sigma,
    double           init_alpha,
    double           init_beta,
    double           init_rho,
    SrtDiffusionType vol_type,
    double*          new_sigmabeta,
    double*          new_alpha,
    double           new_beta,
    double*          new_rho,
    double           nb_std,
    int              solve_on_vol,
    double           precision,
    int              nb_iter_max,
    double*          calib_err);

#endif