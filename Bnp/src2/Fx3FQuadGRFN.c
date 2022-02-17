#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

#define CAP_STDEV 2.0
#define FLOOR_STDEV 2.0

/* transformation functions */
/* if delta is positive */

double Z_to_U_deltaPo_Quad(long i, double Z, LOCAL_MODEL_PARAM *model_param) {

  if (Z > model_param->Zmax[i]) {
    /* Cap on Z */
    return Z / model_param->ZU_U_param1[i] + model_param->ZU_U_param2[i];
  } else if (Z >= model_param->Zmin[i]) {
    return log((Z - model_param->ZU_M_param1) /
               (Z - model_param->ZU_M_param2)) /
               model_param->ZU_M_param3 -
           model_param->ZU_M_param4;
  } else {
    /* Floor on Z */
    return (log(Z) / model_param->ZU_D_param1[i] + model_param->ZU_D_param2[i]);
    // return Z / model_param -> ZU_D_param1[i] + model_param -> ZU_D_param2[i];
  }
}

double U_to_Z_deltaPo_Quad(long i, double U, LOCAL_MODEL_PARAM *model_param) {
  static double x;

  if (U > model_param->Umax[i]) {
    /* Cap on U */
    return (U - model_param->ZU_U_param2[i]) * model_param->ZU_U_param1[i];
  } else if (U >= model_param->Umin[i]) {
    x = exp((U + model_param->ZU_M_param4) * model_param->ZU_M_param3);
    x = (model_param->ZU_M_param2 * x - model_param->ZU_M_param1) / (x - 1.0);
    return x;
  } else {
    /* Floor on U */
    return exp((U - model_param->ZU_D_param2[i]) * model_param->ZU_D_param1[i]);
    // return (U - model_param -> ZU_D_param2[i]) * model_param ->
    // ZU_D_param1[i];
  }
}

/* if delta is negative */

double Z_to_U_deltaNeg_Quad(long i, double Z, LOCAL_MODEL_PARAM *model_param) {
  if (Z > model_param->Zmax[i]) {
    /* Cap on Z */
    return Z / model_param->ZU_U_param1[i] + model_param->ZU_U_param2[i];
  } else if (Z >= model_param->Zmin[i]) {
    return model_param->ZU_M_param1 *
               atan(model_param->ZU_M_param2 * (Z + model_param->ZU_M_param3)) -
           model_param->ZU_M_param4;
  } else {
    /* Floor on Z */
    return log(Z) / model_param->ZU_D_param1[i] + model_param->ZU_D_param2[i];
    // return Z / model_param -> ZU_D_param1[i] + model_param -> ZU_D_param2[i];
  }
}

double U_to_Z_deltaNeg_Quad(long i, double U, LOCAL_MODEL_PARAM *model_param) {
  if (U > model_param->Umax[i]) {
    /* Cap on U */
    return (U - model_param->ZU_U_param2[i]) * model_param->ZU_U_param1[i];
  } else if (U >= model_param->Umin[i]) {
    return tan((U + model_param->ZU_M_param4) / model_param->ZU_M_param1) /
               model_param->ZU_M_param2 -
           model_param->ZU_M_param3;
  } else {
    /* Floor on U */
    return exp((U - model_param->ZU_D_param2[i]) * model_param->ZU_D_param1[i]);
    // return (U - model_param -> ZU_D_param2[i]) * model_param ->
    // ZU_D_param1[i];
  }
}

double Z_to_S_Quad(long i, double Z, double Fwd,
                   LOCAL_MODEL_PARAM *model_param) {
  return Z * Fwd;
}

double S_to_Z_Quad(long i, double S, double Fwd,
                   LOCAL_MODEL_PARAM *model_param) {
  return S / Fwd;
}

double vol_func_Quad(long i, double Z, LOCAL_MODEL_PARAM *model_param) {
  if (Z > model_param->Zmax[i]) {
    return model_param->volU_param1[i];
  } else if (Z >= model_param->Zmin[i]) {
    return (model_param->volM_param1 +
            Z * (model_param->volM_param2 + model_param->volM_param3 * Z));
  } else {
    return model_param->volD_param1[i] * Z;
    // return model_param -> volD_param1[i];
  }
}

double vol_func_deriv_Quad(long i, double Z, LOCAL_MODEL_PARAM *model_param) {
  if (Z > model_param->Zmax[i]) {
    return model_param->volderivU_param1[i];
  } else if (Z >= model_param->Zmin[i]) {
    return (model_param->volderivM_param1 + model_param->volderivM_param2 * Z);
  } else {
    return model_param->volderivD_param1[i];
    // return model_param -> volderivD_param1[i];
  }
}

double vol_func_deriv_t_Quad(long i, double Z, LOCAL_MODEL_PARAM *model_param) {
  return 0.0;
}

double vol_lnATM_func_Quad(double t, double S, double Fwd, double *param) {
  static double Z;

  Z = S / Fwd;
  return (param[0] + Z * (param[1] + param[2] * Z) / Z);
}

double(vol_ln_func_Quad_Mc)(double std, double S, double F, double *params) {
  static double x;

  x = F / S;

  if (x < 2) {
    x = 2;
  } else if (x > -2) {
    x = -2;
  }
  return params[0] * x + (1.0 - x) * (params[1] + params[2] * (1.0 / x - 1.0));
}

Err allocate_local_model(long nbSteps, LOCAL_MODEL_FUNC **model_func,
                         LOCAL_MODEL_PARAM **model_param) {
  Err err = NULL;

  (*model_func) = calloc(1, sizeof(LOCAL_MODEL_FUNC));
  (*model_param) = calloc(1, sizeof(LOCAL_MODEL_PARAM));

  if (!(*model_param) || !(*model_func)) {
    free_local_model_param(*model_func, *model_param);
    err = "Memory allocation error in allocate_local_model_param (1)";
    return err;
  }

  (*model_param)->nbSteps = nbSteps;

  (*model_param)->Zmin = calloc(nbSteps, sizeof(double));
  (*model_param)->Zmax = calloc(nbSteps, sizeof(double));
  (*model_param)->Umin = calloc(nbSteps, sizeof(double));
  (*model_param)->Umax = calloc(nbSteps, sizeof(double));

  (*model_param)->volD_param1 = calloc(nbSteps, sizeof(double));
  (*model_param)->volderivD_param1 = calloc(nbSteps, sizeof(double));
  (*model_param)->ZU_D_param1 = calloc(nbSteps, sizeof(double));
  (*model_param)->ZU_D_param2 = calloc(nbSteps, sizeof(double));

  (*model_param)->volU_param1 = calloc(nbSteps, sizeof(double));
  (*model_param)->volderivU_param1 = calloc(nbSteps, sizeof(double));
  (*model_param)->ZU_U_param1 = calloc(nbSteps, sizeof(double));
  (*model_param)->ZU_U_param2 = calloc(nbSteps, sizeof(double));

  if (!(*model_param)->Zmin || !(*model_param)->Zmax || !(*model_param)->Umin ||
      !(*model_param)->Umax ||

      !(*model_param)->volD_param1 || !(*model_param)->volderivD_param1 ||
      !(*model_param)->ZU_D_param1 || !(*model_param)->ZU_D_param2 ||

      !(*model_param)->volU_param1 || !(*model_param)->volderivU_param1 ||
      !(*model_param)->ZU_U_param1 || !(*model_param)->ZU_U_param2) {
    err = "Memory allocation error in allocate_local_model_param (1)";
    free_local_model_param(*model_func, *model_param);
    return err;
  }

  return err;
}

Err fill_local_model(LOCAL_MODEL_FUNC *model_func,
                     LOCAL_MODEL_PARAM *model_param, long nbSteps, double *time,
                     double alpha, double beta, double gamma, double sig0,
                     double *fwd_fx) {
  Err err = NULL;
  double delta, sq_delta;
  long i;

  model_param->Zmax[0] = 100;
  model_param->Zmin[0] = 1.0E-10;

  model_func->Z_to_S = &(Z_to_S_Quad);
  model_func->S_to_Z = &(S_to_Z_Quad);
  model_func->vol_func = &(vol_func_Quad);
  model_func->vol_func_deriv = &(vol_func_deriv_Quad);
  model_func->vol_func_deriv_t = &(vol_func_deriv_t_Quad);

  model_param->time = time;

  /* First we do the middle function */

  delta = beta * beta - 4.0 * alpha * gamma;

  if (delta >= 0) {
    sq_delta = sqrt(delta);
    model_param->ZU_M_param3 = sq_delta;
    model_param->ZU_M_param1 = (-beta + sq_delta) / (2.0 * gamma);
    model_param->ZU_M_param2 = (-beta - sq_delta) / (2.0 * gamma);

    model_func->Z_to_U = &(Z_to_U_deltaPo_Quad);
    model_func->U_to_Z = &(U_to_Z_deltaPo_Quad);

    model_param->ZU_M_param4 = 0;
    model_param->ZU_M_param4 = Z_to_U_deltaPo_Quad(0, 1.0, model_param);
  } else {
    sq_delta = sqrt(-delta);
    model_param->ZU_M_param1 = 2.0 / sq_delta;
    model_param->ZU_M_param2 = gamma * model_param->ZU_M_param1;
    model_param->ZU_M_param3 = beta / (2.0 * gamma);

    model_func->Z_to_U = &(Z_to_U_deltaNeg_Quad);
    model_func->U_to_Z = &(U_to_Z_deltaNeg_Quad);

    model_param->ZU_M_param4 = 0;
    model_param->ZU_M_param4 = Z_to_U_deltaNeg_Quad(0, 1.0, model_param);
  }

  model_param->volM_param1 = alpha;
  model_param->volM_param2 = beta;
  model_param->volM_param3 = gamma;
  model_param->volderivM_param1 = beta;
  model_param->volderivM_param2 = 2.0 * gamma;

  /* now we do the others */
  for (i = 1; i < nbSteps; i++) {
    // model_param -> Zmax[i]	 = 1.0 + CAP_STDEV * sig0 * sqrt(time[i]);
    // model_param -> Zmin[i]	 = 1.0 - FLOOR_STDEV * sig0 * sqrt(time[i]);

    // model_param -> Zmax[i]	 = 1.0 + CAP_STDEV * sig0 *
    // sqrt(time[nbSteps-1]); model_param -> Zmin[i]	 = 1.0 - FLOOR_STDEV *
    // sig0 * sqrt(time[nbSteps-1]);

    model_param->Zmax[i] = 1E30;
    model_param->Zmin[i] = 1.0e-15;

    if (model_param->Zmin[i] < 0.0) {
      model_param->Zmin[i] = 1.0E-08;
    }

    // model_param -> Zmax[i]	 = 500.0;
    // model_param -> Zmin[i]	 = 1.0E-10;

    model_param->Umax[i] =
        model_func->Z_to_U(i, model_param->Zmax[i], model_param);
    model_param->Umin[i] =
        model_func->Z_to_U(i, model_param->Zmin[i], model_param);

    /* Param vol Up */
    model_param->ZU_U_param1[i] =
        model_func->vol_func(i, model_param->Zmax[i], model_param);

    model_param->ZU_U_param2[i] =
        model_func->Z_to_U(i, model_param->Zmax[i], model_param) -
        model_param->Zmax[i] / model_param->ZU_U_param1[i];

    model_param->volU_param1[i] = model_param->ZU_U_param1[i];
    model_param->volderivU_param1[i] = 0.0;

    /* Param vol Down */

    /* Case Ln */
    model_param->ZU_D_param1[i] =
        model_func->vol_func(i, model_param->Zmin[i], model_param) /
        model_param->Zmin[i];

    model_param->ZU_D_param2[i] =
        model_func->Z_to_U(i, model_param->Zmin[i], model_param) -
        log(model_param->Zmin[i]) / model_param->ZU_D_param1[i];

    model_param->volD_param1[i] = model_param->ZU_D_param1[i];
    model_param->volderivD_param1[i] = model_param->ZU_D_param1[i];

    /* case Floor
    model_param -> ZU_D_param1[i] = model_func -> vol_func(i  , model_param ->
    Zmin[i]  , model_param);

    model_param -> ZU_D_param2[i] = model_func -> Z_to_U(i  , model_param ->
    Zmin[i]  , model_param)
                                                                    -
    model_param -> Zmin[i] / model_param -> ZU_D_param1[i];

    model_param -> volD_param1[i] = model_param -> ZU_D_param1[i];
    model_param -> volderivD_param1[i] = 0.0;
    */
  }

  return err;
}

void free_local_model_param(LOCAL_MODEL_FUNC *model_func,
                            LOCAL_MODEL_PARAM *model_param) {

  if (model_param) {
    if (model_param->Zmin)
      free(model_param->Zmin);
    if (model_param->Zmax)
      free(model_param->Zmax);
    if (model_param->Umin)
      free(model_param->Umin);
    if (model_param->Umax)
      free(model_param->Umax);

    if (model_param->volD_param1)
      free(model_param->volD_param1);
    if (model_param->volderivD_param1)
      free(model_param->volderivD_param1);
    if (model_param->ZU_D_param1)
      free(model_param->ZU_D_param1);
    if (model_param->ZU_D_param2)
      free(model_param->ZU_D_param2);

    if (model_param->volU_param1)
      free(model_param->volU_param1);
    if (model_param->volderivU_param1)
      free(model_param->volderivU_param1);
    if (model_param->ZU_U_param1)
      free(model_param->ZU_U_param1);
    if (model_param->ZU_U_param2)
      free(model_param->ZU_U_param2);

    free(model_param);
    model_param = NULL;
  }

  if (model_func) {
    free(model_func);
    model_func = NULL;
  }
}

Err grfn_payoff_4_3dfxQuad_tree(
    /* Event */
    double evt_date, double evt_time, long step, void *func_parm,
    /* Market data */
    double fwd_fx, void *dom_yc, double dom_lam, double dom_phi, void *for_yc,
    double for_lam, double for_phi,
    /* Nodes data */
    long n1, long n2, long n3,
    /* i: d1  , j: d2  , k: d3  , l = {0: xDom  , 1: xFor  , 2: log (Fx/Fx0)} */
    double ****sv, LOCAL_MODEL_FUNC *model_func, LOCAL_MODEL_PARAM *model_param,
    /* Vector of results to be updated */
    long nprod, double ****prod_val,
    /* Barrier details */
    int is_bar, int bar_k, int **bar_idx, int bar_col, double bar_lvl) {
  GRFNPARMTREE total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;
  double *dom_dff = NULL, *dom_gam = NULL, *dom_gam2 = NULL, *fx_dff = NULL,
         *fx_gam = NULL, *fx_gam2 = NULL, *for_dff = NULL, *for_gam = NULL,
         *for_gam2 = NULL;
  int i, j, k, l;
  double temp;
  unsigned short do_dom, do_for, do_fx;
  double *temp1 = NULL, *temp2 = NULL;
  int start_idx;
  Err err = NULL;

  /* Get the event */
  total = (GRFNPARMTREE)func_parm;
  global = total->global;
  local = total->local;

  /* Precalculate DFF  , gamma and 0.5 * gamma * gamma */

  if (total->num_dom_df > 0 && total->dom_idx != -1) {
    dom_dff = dvector(0, total->num_dom_df - 1);
    dom_gam = dvector(0, total->num_dom_df - 1);
    dom_gam2 = dvector(0, total->num_dom_df - 1);

    if (!dom_dff || !dom_gam || !dom_gam2) {
      err = "Memory allocation error (1) in grfn_payoff_4_3dfx_tree";
      goto FREE_RETURN;
    }

    for (i = 0; i < total->num_dom_df; i++) {
      dom_dff[i] = swp_f_df(evt_date, total->dom_df_dts[i], (char *)dom_yc);
      dom_gam[i] = (1.0 - exp(-dom_lam * total->dom_df_tms[i])) / dom_lam;
      dom_gam2[i] = 0.5 * dom_gam[i] * dom_gam[i] * dom_phi;
    }

    do_dom = 1;
  } else {
    do_dom = 0;
  }

  if (total->num_fx_df > 0 && total->fx_idx != -1) {
    fx_dff = dvector(0, total->num_fx_df - 1);
    fx_gam = dvector(0, total->num_fx_df - 1);
    fx_gam2 = dvector(0, total->num_fx_df - 1);

    if (!fx_dff || !fx_gam || !fx_gam2) {
      err = "Memory allocation error (2) in grfn_payoff_4_3dfx_tree";
      goto FREE_RETURN;
    }

    for (i = 0; i < total->num_fx_df; i++) {
      fx_dff[i] = swp_f_df(evt_date, total->fx_df_dts[i], (char *)dom_yc);
      fx_gam[i] = (1.0 - exp(-dom_lam * total->fx_df_tms[i])) / dom_lam;
      fx_gam2[i] = 0.5 * fx_gam[i] * fx_gam[i] * dom_phi;
    }

    do_fx = 1;
  } else {
    do_fx = 0;
  }

  if (total->num_for_df > 0 && total->for_idx != -1) {
    for_dff = dvector(0, total->num_for_df - 1);
    for_gam = dvector(0, total->num_for_df - 1);
    for_gam2 = dvector(0, total->num_for_df - 1);

    if (!for_dff || !for_gam || !for_gam2) {
      err = "Memory allocation error (3) in grfn_payoff_4_3dfx_tree";
      goto FREE_RETURN;
    }

    for (i = 0; i < total->num_for_df; i++) {
      for_dff[i] = swp_f_df(evt_date, total->for_df_dts[i], (char *)for_yc);
      for_gam[i] = (1.0 - exp(-for_lam * total->for_df_tms[i])) / for_lam;
      for_gam2[i] = 0.5 * for_gam[i] * for_gam[i] * for_phi;
    }

    do_for = 1;
  } else {
    do_for = 0;
  }

  /* Calculate digital jumps and store them in the last column */
  if (is_bar) {
    temp1 = (double *)calloc(nprod, sizeof(double));
    temp2 = (double *)calloc(nprod, sizeof(double));

    if (!temp1 || !temp2) {
      err = "Memory allocation error (3) in grfn_payoff_4_3dfx_tree";
      goto FREE_RETURN;
    }

    switch (bar_k) {
    case 0:

      for (j = 0; j < n2; j++) {
        for (k = 0; k < n3; k++) {
          i = bar_idx[j][k];

          if (i >= 0) {
            if (do_dom) {
              for (l = 0; l < total->num_dom_df; l++) {
                local->evt->df[total->dom_idx][l] =
                    dom_dff[l] *
                    exp(-dom_gam[l] * sv[i][j][k][0] - dom_gam2[l]);
              }
            }
            if (do_fx) {
              for (l = 0; l < total->num_fx_df; l++) {
                local->evt->df[total->fx_idx][l] =
                    fx_dff[l] * exp(-fx_gam[l] * sv[i][j][k][0] - fx_gam2[l]);
              }
            }
            if (do_for) {
              for (l = 0; l < total->num_for_df; l++) {
                local->evt->df[total->for_idx][l] =
                    for_dff[l] *
                    exp(-for_gam[l] * sv[i][j][k][1] - for_gam2[l]);
              }
            }

            local->smp.und[total->fx_idx].sv[SPOT] =
                model_func->Z_to_S(
                    step, model_func->U_to_Z(step, bar_lvl, model_param),
                    fwd_fx, model_param) -
                1.0e-04;
            memcpy(temp1, prod_val[i][j][k], nprod * sizeof(double));
            err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, temp1,
                                 &temp);

            if (err) {
              goto FREE_RETURN;
            }
            local->smp.und[total->fx_idx].sv[SPOT] =
                model_func->Z_to_S(
                    step, model_func->U_to_Z(step, bar_lvl, model_param),
                    fwd_fx, model_param) +
                1.0e-04;
            memcpy(temp2, prod_val[i][j][k], nprod * sizeof(double));
            err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, temp2,
                                 &temp);

            if (err) {
              goto FREE_RETURN;
            }

            for (l = 0; l < n1; l++) {
              prod_val[l][j][k][nprod] = temp2[bar_col] - temp1[bar_col];
            }
          }
        }
      }

      break;

    case 1:

      for (i = 0; i < n1; i++) {
        for (k = 0; k < n3; k++) {
          j = bar_idx[i][k];

          if (j >= 0) {
            if (do_dom) {
              for (l = 0; l < total->num_dom_df; l++) {
                local->evt->df[total->dom_idx][l] =
                    dom_dff[l] *
                    exp(-dom_gam[l] * sv[i][j][k][0] - dom_gam2[l]);
              }
            }
            if (do_fx) {
              for (l = 0; l < total->num_fx_df; l++) {
                local->evt->df[total->fx_idx][l] =
                    fx_dff[l] * exp(-fx_gam[l] * sv[i][j][k][0] - fx_gam2[l]);
              }
            }
            if (do_for) {
              for (l = 0; l < total->num_for_df; l++) {
                local->evt->df[total->for_idx][l] =
                    for_dff[l] *
                    exp(-for_gam[l] * sv[i][j][k][1] - for_gam2[l]);
              }
            }

            local->smp.und[total->fx_idx].sv[SPOT] =
                model_func->Z_to_S(
                    step, model_func->U_to_Z(step, bar_lvl, model_param),
                    fwd_fx, model_param) -
                1.0e-04;
            memcpy(temp1, prod_val[i][j][k], nprod * sizeof(double));
            err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, temp1,
                                 &temp);

            if (err) {
              goto FREE_RETURN;
            }
            local->smp.und[total->fx_idx].sv[SPOT] =
                model_func->Z_to_S(
                    step, model_func->U_to_Z(step, bar_lvl, model_param),
                    fwd_fx, model_param) +
                1.0e-04;
            memcpy(temp2, prod_val[i][j][k], nprod * sizeof(double));
            err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, temp2,
                                 &temp);

            if (err) {
              goto FREE_RETURN;
            }

            for (l = 0; l < n2; l++) {
              prod_val[i][l][k][nprod] = temp2[bar_col] - temp1[bar_col];
            }
          }
        }
      }

      break;

    case 2:

      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          k = bar_idx[i][j];

          if (k >= 0) {
            if (do_dom) {
              for (l = 0; l < total->num_dom_df; l++) {
                local->evt->df[total->dom_idx][l] =
                    dom_dff[l] *
                    exp(-dom_gam[l] * sv[i][j][k][0] - dom_gam2[l]);
              }
            }
            if (do_fx) {
              for (l = 0; l < total->num_fx_df; l++) {
                local->evt->df[total->fx_idx][l] =
                    fx_dff[l] * exp(-fx_gam[l] * sv[i][j][k][0] - fx_gam2[l]);
              }
            }
            if (do_for) {
              for (l = 0; l < total->num_for_df; l++) {
                local->evt->df[total->for_idx][l] =
                    for_dff[l] *
                    exp(-for_gam[l] * sv[i][j][k][1] - for_gam2[l]);
              }
            }

            local->smp.und[total->fx_idx].sv[SPOT] =
                model_func->Z_to_S(
                    step, model_func->U_to_Z(step, bar_lvl, model_param),
                    fwd_fx, model_param) -
                1.0e-04;
            memcpy(temp1, prod_val[i][j][k], nprod * sizeof(double));
            err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, temp1,
                                 &temp);

            if (err) {
              goto FREE_RETURN;
            }
            local->smp.und[total->fx_idx].sv[SPOT] =
                model_func->Z_to_S(
                    step, model_func->U_to_Z(step, bar_lvl, model_param),
                    fwd_fx, model_param) +
                1.0e-04;
            memcpy(temp2, prod_val[i][j][k], nprod * sizeof(double));
            err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, temp2,
                                 &temp);

            if (err) {
              goto FREE_RETURN;
            }

            for (l = 0; l < n3; l++) {
              prod_val[i][j][l][nprod] = temp2[bar_col] - temp1[bar_col];
            }
          }
        }
      }

      break;

    default:
      break;
    }
  }

  /* Eval payoff */
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) {
      for (k = 0; k < n3; k++) {
        local->smp.und[total->fx_idx].sv[SPOT] =
            model_func->Z_to_S(
                step, model_func->U_to_Z(step, sv[i][j][k][2], model_param),
                fwd_fx, model_param) +
            1.0E-08;

        if (do_dom) {
          for (l = 0; l < total->num_dom_df; l++) {
            local->evt->df[total->dom_idx][l] =
                dom_dff[l] * exp(-dom_gam[l] * sv[i][j][k][0] - dom_gam2[l]);
          }
        }
        if (do_fx) {
          for (l = 0; l < total->num_fx_df; l++) {
            local->evt->df[total->fx_idx][l] =
                fx_dff[l] * exp(-fx_gam[l] * sv[i][j][k][0] - fx_gam2[l]);
          }
        }
        if (do_for) {
          for (l = 0; l < total->num_for_df; l++) {
            local->evt->df[total->for_idx][l] =
                for_dff[l] * exp(-for_gam[l] * sv[i][j][k][1] - for_gam2[l]);
          }
        }
        err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL,
                             prod_val[i][j][k], &temp);

        if (err) {
          goto FREE_RETURN;
        }
      }
    }
  }

  /* Remove digital jumps */
  if (is_bar) {
    switch (bar_k) {
    case 0:

      for (j = 0; j < n2; j++) {
        for (k = 0; k < n3; k++) {
          start_idx = bar_idx[j][k];
          if (start_idx >= 0) {
            if (sv[start_idx][j][k][2] < bar_lvl - 1.0e-08) {
              start_idx++;
            }

            for (i = start_idx; i < n1; i++) {
              prod_val[i][j][k][bar_col] -= prod_val[i][j][k][nprod];
            }
          }
        }
      }

      break;

    case 1:

      for (i = 0; i < n1; i++) {
        for (k = 0; k < n3; k++) {
          start_idx = bar_idx[i][k];
          if (start_idx >= 0) {
            if (sv[i][start_idx][k][2] < bar_lvl - 1.0e-08) {
              start_idx++;
            }

            for (j = start_idx; j < n2; j++) {
              prod_val[i][j][k][bar_col] -= prod_val[i][j][k][nprod];
            }
          }
        }
      }

      break;

    case 2:

      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          start_idx = bar_idx[i][j];
          if (start_idx >= 0) {
            if (sv[i][j][start_idx][2] < bar_lvl - 1.0e-08) {
              start_idx++;
            }

            for (k = start_idx; k < n3; k++) {
              prod_val[i][j][k][bar_col] -= prod_val[i][j][k][nprod];
            }
          }
        }
      }

      break;

    default:
      break;
    }
  }

FREE_RETURN:

  if (temp1) {
    free(temp1);
  }

  if (temp2) {
    free(temp2);
  }

  if (dom_dff) {
    free_dvector(dom_dff, 0, total->num_dom_df - 1);
  }

  if (dom_gam) {
    free_dvector(dom_gam, 0, total->num_dom_df - 1);
  }

  if (dom_gam2) {
    free_dvector(dom_gam2, 0, total->num_dom_df - 1);
  }

  if (for_dff) {
    free_dvector(for_dff, 0, total->num_for_df - 1);
  }

  if (for_gam) {
    free_dvector(for_gam, 0, total->num_for_df - 1);
  }

  if (for_gam2) {
    free_dvector(for_gam2, 0, total->num_for_df - 1);
  }

  if (fx_dff) {
    free_dvector(fx_dff, 0, total->num_fx_df - 1);
  }

  if (fx_gam) {
    free_dvector(fx_gam, 0, total->num_fx_df - 1);
  }

  if (fx_gam2) {
    free_dvector(fx_gam2, 0, total->num_fx_df - 1);
  }

  return err;
}

Err grfn_payoff_4_3dfx_Quadmc(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double spot_fx, void *dom_yc, double dom_lam, double dom_phi, void *for_yc,
    double for_lam, double for_phi, double Xdom, double Yfor, double Zfx,
    /* Results */
    int num_col, double *res, int *stop_path) {
  GRFNPARMMC total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;
  int l;
  double temp;
  Err err;

  if (func_parm == NULL) {
    memset(res, 0, num_col * sizeof(double));
    *stop_path = 0;
    return NULL;
  } else {
    err = NULL;
  }

  memset(res, 0, num_col * sizeof(double));

  /* Get the event */
  total = (GRFNPARMMC)func_parm;
  global = total->global;
  local = total->local;

  /* Calc market data */
  local->smp.und[total->fx_idx].sv[SPOT] = exp(Zfx);

  if (total->do_dom) {
    for (l = 0; l < total->num_dom_df; l++) {
      local->evt->df[total->dom_idx][l] =
          total->dom_dff[l] *
          exp(-total->dom_gam[l] * Xdom - total->dom_gam2[l]);
    }
  }
  if (total->do_fx) {
    for (l = 0; l < total->num_fx_df; l++) {
      local->evt->df[total->fx_idx][l] =
          total->fx_dff[l] * exp(-total->fx_gam[l] * Xdom - total->fx_gam2[l]);
    }
  }
  if (total->do_for) {
    for (l = 0; l < total->num_for_df; l++) {
      local->evt->df[total->for_idx][l] =
          total->for_dff[l] *
          exp(-total->for_gam[l] * Yfor - total->for_gam2[l]);
    }
  }

  err = FIRSTEvalEvent(global, local, num_col, 2, NULL, NULL, res, &temp);

  *stop_path = 0;
  return err;
}

Err FxCall_payoff_4_3dfxQuad_tree(
    /* Event */
    double evt_date, double evt_time, long step, void *func_parm,
    /* Market data */
    double fwd_fx, void *dom_yc, double dom_lam, double dom_phi, void *for_yc,
    double for_lam, double for_phi,
    /* Nodes data */
    long n1, long n2, long n3,
    /* i: d1  , j: d2  , k: d3  , l = {0: xDom  , 1: xFor  , 2: log (Fx/Fx0)} */
    double ****sv, LOCAL_MODEL_FUNC *model_func, LOCAL_MODEL_PARAM *model_param,
    /* Vector of results to be updated */
    long nprod, double ****prod_val,
    /* Barrier details */
    int is_bar, int bar_k, int **bar_idx, int bar_col, double bar_lvl) {

  long i, j, k, n;
  double Fx;
  double *strikes;

  strikes = (double *)(func_parm);

  /* Eval payoff */
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) {
      for (k = 0; k < n3; k++) {
        Fx = model_func->Z_to_S(
            step, model_func->U_to_Z(step, sv[i][j][k][2], model_param), fwd_fx,
            model_param);

        for (n = 0; n < nprod; n++) {
          if (Fx > strikes[n]) {
            prod_val[i][j][k][n] = Fx - strikes[n];
          } else {
            prod_val[i][j][k][n] = 0;
          }
        }
      }
    }
  }

  return NULL;
}