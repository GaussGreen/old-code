#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

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
                exp((log((1.0 - beta) * bar_lvl + 1.0)) / (1.0 - beta)) -
                1.0e-04;
            memcpy(temp1, prod_val[i][j][k], nprod * sizeof(double));
            err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, temp1,
                                 &temp);

            if (err) {
              goto FREE_RETURN;
            }
            local->smp.und[total->fx_idx].sv[SPOT] =
                exp((log((1.0 - beta) * bar_lvl + 1.0)) / (1.0 - beta)) +
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
                exp((log((1.0 - beta) * bar_lvl + 1.0)) / (1.0 - beta)) -
                1.0e-04;
            memcpy(temp1, prod_val[i][j][k], nprod * sizeof(double));
            err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, temp1,
                                 &temp);

            if (err) {
              goto FREE_RETURN;
            }
            local->smp.und[total->fx_idx].sv[SPOT] =
                exp((log((1.0 - beta) * bar_lvl + 1.0)) / (1.0 - beta)) +
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
                exp((log((1.0 - beta) * bar_lvl + 1.0)) / (1.0 - beta)) -
                1.0e-04;
            memcpy(temp1, prod_val[i][j][k], nprod * sizeof(double));
            err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, temp1,
                                 &temp);

            if (err) {
              goto FREE_RETURN;
            }
            local->smp.und[total->fx_idx].sv[SPOT] =
                exp((log((1.0 - beta) * bar_lvl + 1.0)) / (1.0 - beta)) +
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
            exp((log((1.0 - beta) * sv[i][j][k][2] + 1.0)) / (1.0 - beta)) +
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

Err grfn_payoff_4_3dfx_Betamc(
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