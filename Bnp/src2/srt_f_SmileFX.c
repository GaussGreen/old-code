#include "math.h"
#include "srt_h_BrownianMotion.h"
#include "srt_h_FeynmannKac.h"
#include "srt_h_all.h"
#include "srt_h_fwdcurve.h"

Err srt_f_SABW_Opt(
    Date   eval_date,
    Date   exp_date,
    double option_strike,

    char* rec_pay_str,
    char* greek_str,

    char* yield_crv_name,
    char* grow_crv_name,

    /* LOCAL VOLATILITY PARAMETERS */
    double basevol,
    double beta,
    double gamma,
    double omega,

    /* PDE PARAMETERS */
    long   min_node,
    long   min_num_mesh,
    double X0,
    double X_min,
    double X_max,
    double h_max,

    double* option_pv)
{
    Err    err = NULL;
    long   NumTimeStep, i, j, k, NumMesh;
    double yr_to_exp_date, TimeStep = 1.0 / 24, dX, PDEBoundary[2], *X = NULL, *exp_X = NULL,
                           *TerminalPayoff = NULL, **D = NULL, **V = NULL, *U = NULL,
                           **div_mat = NULL;
    double time, std, init_fwd_val;
    double df_exp, drift;

    SrtPdePDEType      PDEType         = LINEAR;
    SrtPdeBoundaryCond PDEBoundaryCond = DIRICHLET;
    SrtCurvePtr        yldcrv;
    SrtCurvePtr        fwdcrv;
    SrtGreekType       greek_type;
    SrtReceiverType    RecPay;

    yr_to_exp_date = (double)(exp_date - eval_date) * YEARS_IN_DAY;

    TimeStep    = min(TimeStep, yr_to_exp_date / min_node);
    NumTimeStep = ((long)(yr_to_exp_date / TimeStep));

    /* COMPUTE THE SPACE STEP */
    dX = (X_max - X_min) / min_num_mesh;
    dX = min(dX, h_max);

    /* ADAPT dX SO THAT X0 IS A MESH POINT */
    k       = (long)((X0 - X_min) / dX) + 1;
    dX      = (X0 - X_min) / k;
    NumMesh = (long)((X_max - X_min) / dX);

    X              = dvector(1, NumMesh - 1);
    exp_X          = dvector(1, NumMesh - 1);
    U              = dvector(1, NumMesh - 1);
    TerminalPayoff = dvector(1, NumMesh - 1);
    V              = dmatrix(0, NumTimeStep, 1, NumMesh);
    D              = dmatrix(0, NumTimeStep, 1, NumMesh);
    div_mat        = dmatrix(0, NumTimeStep, 1, NumMesh);

    fwdcrv = lookup_curve(grow_crv_name);
    yldcrv = lookup_curve(yield_crv_name);

    df_exp       = swp_f_df(eval_date, exp_date, yldcrv);
    init_fwd_val = exp(X0) * srt_f_forward_from_fwdcrv(yr_to_exp_date, fwdcrv, NULL) / df_exp;
    drift        = log(init_fwd_val / exp(X0)) / yr_to_exp_date;

    err = interp_rec_pay(rec_pay_str, &RecPay);
    if (err)
        return err;

    /* FILL X, D, V , TerminalPayoff, AND PDEBoundary */
    for (j = 1; j < NumMesh; j++)
    {
        X[j]     = X_min + j * dX;
        exp_X[j] = exp(X[j]);
        if (RecPay == SRT_RECEIVER)
            TerminalPayoff[j] = max(option_strike - exp(X[j]), 0);
        else if (RecPay == SRT_PAYER)
            TerminalPayoff[j] = max(exp(X[j]) - option_strike, 0);
        else if (RecPay == SRT_STRADDLE)
            TerminalPayoff[j] =
                max(option_strike - exp(X[j]), 0) + max(exp(X[j]) - option_strike, 0);

        for (i = 0; i <= NumTimeStep; i++)
        {
            if (i == 0)
                time = yr_to_exp_date;
            else
                time = max(TimeStep, yr_to_exp_date - (i - 1) * TimeStep);

            std = log(exp_X[j] / exp(X0)) / (basevol * sqrt(time));

            V[i][j] = basevol * (1 + beta * omega * tanh(std / omega) +
                                 gamma * omega * omega * tanh(std / omega) * tanh(std / omega));
            D[i][j] = drift - V[i][j] * V[i][j] / 2;
        }
    }

    PDEBoundary[0] = max(option_strike - exp(X_min), 0);
    PDEBoundary[1] = max(option_strike - exp(X_max), 0);

    err = srt_f_feynmann_kac_pde(
        NumMesh,
        NumTimeStep,
        dX,
        TimeStep,

        X,
        D,
        V,
        PDEType,
        PDEBoundaryCond,
        PDEBoundary,
        TerminalPayoff,
        NULL,
        div_mat,
        &U);
    if (err)
        return err;

    err = interp_greeks(greek_str, &greek_type);
    if (err)
        return err;

    if (greek_type == PREMIUM)
        (*option_pv) = df_exp * U[k];
    else if (greek_type == DELTA)
        (*option_pv) = df_exp * (U[k + 1] - U[k]) / dX;
    else if (greek_type == GAMMA)
        (*option_pv) = df_exp * (U[k + 1] - 2 * U[k] + U[k - 1]) / (dX * dX);

    if (X)
        free_dvector(X, 1, NumMesh - 1);
    X = NULL;
    if (exp_X)
        free_dvector(exp_X, 1, NumMesh - 1);
    exp_X = NULL;
    if (U)
        free_dvector(U, 1, NumMesh - 1);
    U = NULL;
    if (TerminalPayoff)
        free_dvector(TerminalPayoff, 1, NumMesh - 1);
    TerminalPayoff = NULL;
    if (V)
        free_matrix(V, 0, NumTimeStep, 1, NumMesh);
    V = NULL;
    if (D)
        free_matrix(D, 0, NumTimeStep, 1, NumMesh);
    D = NULL;
    if (div_mat)
        free_dmatrix(div_mat, 0, NumTimeStep, 1, NumMesh);
    div_mat = NULL;

    return err;
}
