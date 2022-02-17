
/* ==========================================================================
   FILE_NAME:	3FQuantoGRFN.c

   PURPOSE:		Functions for the GRFN payoff evaluation in a 3FQuanto
   model

   DATE:		26 Sep 2003

   AUTHOR:		J.M.L.
   ========================================================================== */

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_all3FQuanto.h"

Err grfn_payoff_4_3fquanto_tree(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    void *lgm2F_yc, double r1_lambda, double r2_lambda, double r1_phi,
    double r2_phi, double r1_r2_phi, double rho, void *lgm1F_yc,
    double r3_lambda, double r3_phi, int lgm2F_is_dom_for,
    /* Nodes data */
    long n1, long n2, long n3,
    /* i: d1  , j: d2  , k: d3  , l = {0: r1  , 1: r2  , 2: r3} */
    double ****sv,
    /* Vector of results to be updated */
    long nprod, double ****prod_val) {
  GRFNPARMTREE total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;
  double *dom_dff = NULL,
         *dom_1f_gam = NULL, // First factor for DOM
             *dom_1f_gam2 = NULL,
         *dom_2f_gam = NULL, // Second factor for DOM (if DOM is LGM2F)
             *dom_2f_gam2 = NULL, *dom_cross_gam = NULL,

         *fx_dff = NULL,
         *fx_1f_gam = NULL, // First factor for FX
             *fx_1f_gam2 = NULL,
         *fx_2f_gam = NULL, // Second factor for FX (if FX = DOM is LGM2F)
             *fx_2f_gam2 = NULL, *fx_cross_gam = NULL,

         *for_dff = NULL,
         *for_1f_gam = NULL, // First factor for FOR
             *for_1f_gam2 = NULL,
         *for_2f_gam = NULL, // Second factor for FOR (if FOR is LGM2F)
             *for_2f_gam2 = NULL, *for_cross_gam = NULL;

  int i, j, k, l;
  double temp;
  unsigned short do_dom, do_for, do_fx;
  Err err = NULL;

  /* Get the event */
  total = (GRFNPARMTREE)func_parm;
  global = total->global;
  local = total->local;

  /* Precalculate DFF  , gamma and 0.5 * gamma * gamma and cross gammas */

  // For DOMESTIC DFs
  if (total->num_dom_df > 0 && total->dom_idx != -1) {
    dom_dff = dvector(0, total->num_dom_df - 1);
    dom_1f_gam = dvector(0, total->num_dom_df - 1);
    dom_1f_gam2 = dvector(0, total->num_dom_df - 1);

    if (!dom_dff || !dom_1f_gam || !dom_1f_gam2) {
      err = "Memory allocation error (1.1) in grfn_payoff_4_3fquanto_tree";
      goto FREE_RETURN;
    }

    if (lgm2F_is_dom_for == 0) {
      dom_2f_gam = dvector(0, total->num_dom_df - 1);
      dom_2f_gam2 = dvector(0, total->num_dom_df - 1);
      dom_cross_gam = dvector(0, total->num_dom_df - 1);
      if (!dom_2f_gam || !dom_2f_gam2 || !dom_cross_gam) {
        err = "Memory allocation error (1.2) in grfn_payoff_4_3fquanto_tree";
        goto FREE_RETURN;
      }
    }

    for (i = 0; i < total->num_dom_df; i++) {
      if (lgm2F_is_dom_for == 0) {
        // LGM2F is DOMESTIC
        dom_dff[i] = swp_f_df(evt_date, total->dom_df_dts[i], (char *)lgm2F_yc);
        dom_1f_gam[i] =
            (1.0 - exp(-r1_lambda * total->dom_df_tms[i])) / r1_lambda;
        dom_1f_gam2[i] = 0.5 * dom_1f_gam[i] * dom_1f_gam[i] * r1_phi;
        dom_2f_gam[i] =
            (1.0 - exp(-r2_lambda * total->dom_df_tms[i])) / r2_lambda;
        dom_2f_gam2[i] = 0.5 * dom_2f_gam[i] * dom_2f_gam[i] * r2_phi;
        dom_cross_gam[i] = rho * dom_1f_gam[i] * dom_2f_gam[i] * r1_r2_phi;
      } else {
        // LGM2F is FOREIGN
        dom_dff[i] = swp_f_df(evt_date, total->dom_df_dts[i], (char *)lgm1F_yc);
        dom_1f_gam[i] =
            (1.0 - exp(-r3_lambda * total->dom_df_tms[i])) / r3_lambda;
        dom_1f_gam2[i] = 0.5 * dom_1f_gam[i] * dom_1f_gam[i] * r3_phi;
      }
    }
    do_dom = 1;
  } else {
    do_dom = 0;
  }

  // For FX DFs = DOMESTIC DFs
  if (total->num_fx_df > 0 && total->fx_idx != -1) {
    fx_dff = dvector(0, total->num_fx_df - 1);
    fx_1f_gam = dvector(0, total->num_fx_df - 1);
    fx_1f_gam2 = dvector(0, total->num_fx_df - 1);

    if (!fx_dff || !fx_1f_gam || !fx_1f_gam2) {
      err = "Memory allocation error (2.1) in grfn_payoff_4_3dfx_tree";
      goto FREE_RETURN;
    }

    if (lgm2F_is_dom_for == 0) {
      fx_2f_gam = dvector(0, total->num_fx_df - 1);
      fx_2f_gam2 = dvector(0, total->num_fx_df - 1);
      fx_cross_gam = dvector(0, total->num_fx_df - 1);
      if (!fx_2f_gam || !fx_2f_gam2 || !fx_cross_gam) {
        err = "Memory allocation error (2.2) in grfn_payoff_4_3fquanto_tree";
        goto FREE_RETURN;
      }
    }

    for (i = 0; i < total->num_fx_df; i++) {

      if (lgm2F_is_dom_for == 0) {
        // LGM2F is DOMESTIC
        fx_dff[i] = swp_f_df(evt_date, total->fx_df_dts[i], (char *)lgm2F_yc);
        fx_1f_gam[i] =
            (1.0 - exp(-r1_lambda * total->fx_df_tms[i])) / r1_lambda;
        fx_1f_gam2[i] = 0.5 * fx_1f_gam[i] * fx_1f_gam[i] * r1_phi;
        fx_2f_gam[i] =
            (1.0 - exp(-r2_lambda * total->fx_df_tms[i])) / r2_lambda;
        fx_2f_gam2[i] = 0.5 * fx_2f_gam[i] * fx_2f_gam[i] * r2_phi;
        fx_cross_gam[i] = rho * fx_1f_gam[i] * fx_2f_gam[i] * r1_r2_phi;
      } else {
        // LGM2F is FOREIGN
        fx_dff[i] = swp_f_df(evt_date, total->fx_df_dts[i], (char *)lgm1F_yc);
        fx_1f_gam[i] =
            (1.0 - exp(-r3_lambda * total->fx_df_tms[i])) / r3_lambda;
        fx_1f_gam2[i] = 0.5 * fx_1f_gam[i] * fx_1f_gam[i] * r3_phi;
      }
    }
    do_fx = 1;
  } else {
    do_fx = 0;
  }

  // For FOREIGN DFs
  if (total->num_for_df > 0 && total->for_idx != -1) {
    for_dff = dvector(0, total->num_for_df - 1);
    for_1f_gam = dvector(0, total->num_for_df - 1);
    for_1f_gam2 = dvector(0, total->num_for_df - 1);

    if (!for_dff || !for_1f_gam || !for_1f_gam2) {
      err = "Memory allocation error (3.1) in grfn_payoff_4_3fquanto_tree";
      goto FREE_RETURN;
    }

    if (lgm2F_is_dom_for == 1) {
      for_2f_gam = dvector(0, total->num_for_df - 1);
      for_2f_gam2 = dvector(0, total->num_for_df - 1);
      for_cross_gam = dvector(0, total->num_for_df - 1);
      if (!for_2f_gam || !for_2f_gam2 || !for_cross_gam) {
        err = "Memory allocation error (3.2) in grfn_payoff_4_3fquanto_tree";
        goto FREE_RETURN;
      }
    }

    for (i = 0; i < total->num_for_df; i++) {
      if (lgm2F_is_dom_for == 1) {
        // LGM2F is FOREIGN
        for_dff[i] = swp_f_df(evt_date, total->for_df_dts[i], (char *)lgm2F_yc);
        for_1f_gam[i] =
            (1.0 - exp(-r1_lambda * total->for_df_tms[i])) / r1_lambda;
        for_1f_gam2[i] = 0.5 * for_1f_gam[i] * for_1f_gam[i] * r1_phi;
        for_2f_gam[i] =
            (1.0 - exp(-r2_lambda * total->for_df_tms[i])) / r2_lambda;
        for_2f_gam2[i] = 0.5 * for_2f_gam[i] * for_2f_gam[i] * r2_phi;
        for_cross_gam[i] = rho * for_1f_gam[i] * for_2f_gam[i] * r1_r2_phi;
      } else {
        // LGM2F is DOMESTIC
        for_dff[i] = swp_f_df(evt_date, total->for_df_dts[i], (char *)lgm1F_yc);
        for_1f_gam[i] =
            (1.0 - exp(-r3_lambda * total->for_df_tms[i])) / r3_lambda;
        for_1f_gam2[i] = 0.5 * for_1f_gam[i] * for_1f_gam[i] * r3_phi;
      }
    }
    do_for = 1;
  } else {
    do_for = 0;
  }

  /* Eval payoff */
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) {
      for (k = 0; k < n3; k++) {
        // local->smp.und[total->fx_idx].sv[SPOT] = spot_fx * exp
        // (sv[i][j][k][2]) + 1.0e-12;

        if (do_dom) {
          for (l = 0; l < total->num_dom_df; l++) {
            if (lgm2F_is_dom_for == 0) {
              local->evt->df[total->dom_idx][l] =
                  dom_dff[l] *
                  exp(-dom_1f_gam[l] * sv[i][j][k][0] - dom_1f_gam2[l] -
                      dom_2f_gam[l] * sv[i][j][k][1] - dom_2f_gam2[l] -
                      dom_cross_gam[l]);
            } else {
              local->evt->df[total->dom_idx][l] =
                  dom_dff[l] *
                  exp(-dom_1f_gam[l] * sv[i][j][k][2] - dom_1f_gam2[l]);
            }
          }
        }

        if (do_fx) {
          for (l = 0; l < total->num_fx_df; l++) {
            if (lgm2F_is_dom_for == 0) {
              local->evt->df[total->fx_idx][l] =
                  fx_dff[l] *
                  exp(-fx_1f_gam[l] * sv[i][j][k][0] - fx_1f_gam2[l] -
                      fx_2f_gam[l] * sv[i][j][k][1] - fx_2f_gam2[l] -
                      fx_cross_gam[l]);
            } else {
              local->evt->df[total->fx_idx][l] =
                  fx_dff[l] *
                  exp(-fx_1f_gam[l] * sv[i][j][k][2] - fx_1f_gam2[l]);
            }
          }
        }

        if (do_for) {
          for (l = 0; l < total->num_for_df; l++) {
            if (lgm2F_is_dom_for == 0) {
              local->evt->df[total->for_idx][l] =
                  for_dff[l] *
                  exp(-for_1f_gam[l] * sv[i][j][k][2] - for_1f_gam2[l]);

            } else {
              local->evt->df[total->for_idx][l] =
                  for_dff[l] *
                  exp(-for_1f_gam[l] * sv[i][j][k][0] - for_1f_gam2[l] -
                      for_2f_gam[l] * sv[i][j][k][1] - for_2f_gam2[l] -
                      for_cross_gam[l]);
            }
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

FREE_RETURN:

  if (dom_dff)
    free_dvector(dom_dff, 0, total->num_dom_df - 1);
  if (dom_1f_gam)
    free_dvector(dom_1f_gam, 0, total->num_dom_df - 1);
  if (dom_1f_gam2)
    free_dvector(dom_1f_gam2, 0, total->num_dom_df - 1);
  if (dom_2f_gam)
    free_dvector(dom_2f_gam, 0, total->num_dom_df - 1);
  if (dom_2f_gam2)
    free_dvector(dom_2f_gam2, 0, total->num_dom_df - 1);
  if (dom_cross_gam)
    free_dvector(dom_cross_gam, 0, total->num_dom_df - 1);

  if (fx_dff)
    free_dvector(fx_dff, 0, total->num_fx_df - 1);
  if (fx_1f_gam)
    free_dvector(fx_1f_gam, 0, total->num_fx_df - 1);
  if (fx_1f_gam2)
    free_dvector(fx_1f_gam2, 0, total->num_fx_df - 1);
  if (fx_2f_gam)
    free_dvector(fx_2f_gam, 0, total->num_fx_df - 1);
  if (fx_2f_gam2)
    free_dvector(fx_2f_gam2, 0, total->num_fx_df - 1);
  if (fx_cross_gam)
    free_dvector(fx_cross_gam, 0, total->num_fx_df - 1);

  if (for_dff)
    free_dvector(for_dff, 0, total->num_for_df - 1);
  if (for_1f_gam)
    free_dvector(for_1f_gam, 0, total->num_for_df - 1);
  if (for_1f_gam2)
    free_dvector(for_1f_gam2, 0, total->num_for_df - 1);
  if (for_2f_gam)
    free_dvector(for_2f_gam, 0, total->num_for_df - 1);
  if (for_2f_gam2)
    free_dvector(for_2f_gam2, 0, total->num_for_df - 1);
  if (for_cross_gam)
    free_dvector(for_cross_gam, 0, total->num_for_df - 1);

  return err;
}
