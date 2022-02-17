/**********************************************************************
 *      Name: SrtGrfnMainFxBetaDLM.cxx                                  *
 *  Function: Entry point to GRFN with raw data                       *
 *            in the case of the 3Factor Beta DLM model * Copyright: (C) BNP
 *Paribas Capital Markets Ltd.                    *
 *--------------------------------------------------------------------*
 *    Author: Laurent Carlier	                                      *
 *      Date: 05/01/00                                                *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere (Excel or 2020)                  *
 *   Returns:                                                         *
 *   Globals: Expects mkt and request list structures to exist        *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 05/01/00 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/

#include "Fx3FBetaDLMCalculations.h"
#include "Fx3FBetaDLMCalibration.h"
#include "Fx3FBetaDLMGRFN.h"
#include "Fx3FBetaDLMMC.h"
#include "Fx3FBetaDLMTree.h"
#include "Fx3FBetaDLMUtil.h"
#include "SrtAccess.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

char *SrtGrfnFx3DFxBetaDLMTreeQTStar(
    char *und3dfx, int numeventdates, long *eventdates, long tableauRows,
    long tableauCols, char ***tableauStrings, int **tableauMask, long auxWidth,
    long *auxLen, double **aux, FxBetaDLM_GRFNNumerParams *Params,
    double *barrier, int *bar_col, long num_stp, int *num_prod,
    double **prod_val);

char *SrtGrfnFx3DFxBetaDLMTreeQBeta(char *und3dfx, int numeventdates,
                                    long *eventdates, long tableauRows,
                                    long tableauCols, char ***tableauStrings,
                                    int **tableauMask, long auxWidth,
                                    long *auxLen, double **aux,
                                    FxBetaDLM_GRFNNumerParams *Params,
                                    double *barrier, int *bar_col, long num_stp,
                                    int *num_prod, double **prod_val);

char *SrtGrfnFx3DFxBetaDLMTree(char *und3dfx, int numeventdates,
                               long *eventdates, long tableauRows,
                               long tableauCols, char ***tableauStrings,
                               int **tableauMask, long auxWidth, long *auxLen,
                               double **aux, FxBetaDLM_GRFNNumerParams *Params,
                               double *barrier, int *bar_col, long num_stp,
                               int *num_prod, double **prod_val) {
  Err err = NULL;

  if (Params->numeraire == 0) {
    err = SrtGrfnFx3DFxBetaDLMTreeQBeta(
        und3dfx, numeventdates, eventdates, tableauRows, tableauCols,
        tableauStrings, tableauMask, auxWidth, auxLen, aux, Params, barrier,
        bar_col, num_stp, num_prod, prod_val);
  } else {
    err = SrtGrfnFx3DFxBetaDLMTreeQTStar(
        und3dfx, numeventdates, eventdates, tableauRows, tableauCols,
        tableauStrings, tableauMask, auxWidth, auxLen, aux, Params, barrier,
        bar_col, num_stp, num_prod, prod_val);
  }

  return err;
}

Fx3DBetaDLM_DriftVolAndCorr(long nstp, double *time, FxBetaDLM_model *model,

                            double *dom_phi, double *dom_beta, double *for_phi,
                            double *for_beta, double *var_ffx,

                            double *dr_const_dom, double *dr_coef_dom,
                            double *dr_const_for, double *dr_coef1_for,
                            double *dr_coef2_for, double *dr_coef3_for,
                            double *dr_const_fx, double *dr_coef1_fx,
                            double *dr_coef2_fx, double *dr_coef3_fx,

                            int *vol_change, double *sig_dom, double *sig_for,
                            double *sig_fx, double *corr_dom_for,
                            double *corr_dom_fx, double *corr_for_fx) {
  int i, index;
  double t_mid;
  double dom_lam, for_lam;
  double vol_bond_dom, vol_bond_for;
  double small_bracket, big_bracket, coef_B, coef_C, coef_D;

  dom_lam = model->dLambdaDom;
  for_lam = model->dLambdaFor;

  vol_change[nstp - 1] = 1;
  index = Get_Index(time[nstp - 1], model->dPWTime, model->iNbPWTime);

  for (i = nstp - 1; i >= 0; i--) {
    /*
    if (i < nstp-1)
    {
            t_mid = 0.5 * (time[i] + time[i+1]);
    }
    else
    {
            t_mid = 0.5 * (time[i] + time[i-1]);
    }
    */

    t_mid = time[i];

    index = Get_Index(t_mid, model->dPWTime, model->iNbPWTime);

    sig_dom[i] = model->dSigmaDom[index];
    sig_for[i] = model->dSigmaFor[index];
    sig_fx[i] = model->dSigmaFx[index];

    corr_dom_for[i] = model->dCorrDomFor[index];
    corr_dom_fx[i] = model->dCorrDomFx[index];
    corr_for_fx[i] = model->dCorrForFx[index];

    if ((i == nstp - 1) ||
        (fabs(sig_dom[i] - sig_dom[i + 1]) + fabs(sig_for[i] - sig_for[i + 1]) +
             fabs(sig_fx[i] - sig_fx[i + 1]) +
             fabs(corr_dom_for[i] - corr_dom_for[i + 1]) +
             fabs(corr_dom_fx[i] - corr_dom_fx[i + 1]) +
             fabs(corr_for_fx[i] - corr_for_fx[i + 1]) >
         EPS)) {
      vol_change[i] = 1;
    } else {
      vol_change[i] = 0;
    }

    /* Get the drift */

    if (i < nstp - 1) {
      vol_bond_dom = sig_dom[i + 1] * dom_beta[i + 1];
      vol_bond_for = sig_for[i + 1] * for_beta[i + 1];

      dr_const_dom[i] = dom_phi[i];
      dr_coef_dom[i] = -dom_lam;

      coef_B = 1 / (1.0 + 2.0 * model->dC0 * var_ffx[i]);
      coef_C = model->dC0 * coef_B;
      coef_B *= model->dB0;
      coef_D = 0.5 * (var_ffx[i] + dom_beta[i] * dom_beta[i] * dom_phi[i] -
                      for_beta[i] * for_beta[i] * for_phi[i]);

      small_bracket = -vol_bond_for + corr_dom_for[i + 1] * vol_bond_dom;
      big_bracket = small_bracket + corr_for_fx[i + 1] * sig_fx[i + 1];

      small_bracket *= sig_for[i + 1];
      big_bracket *= sig_for[i + 1];

      dr_const_for[i] = for_phi[i] -
                        (coef_B + 2.0 * coef_C * coef_D) * big_bracket +
                        small_bracket;
      dr_coef1_for[i] = -2.0 * coef_C * dom_beta[i] * big_bracket;
      dr_coef2_for[i] = -for_lam + 2.0 * coef_C * for_beta[i] * big_bracket;
      dr_coef3_for[i] = -2.0 * coef_C * big_bracket;

      big_bracket *= -for_beta[i + 1];

      dr_const_fx[i] = (coef_B + 2.0 * coef_C * coef_D - 1.0) * big_bracket -
                       0.5 * sig_fx[i + 1] * sig_fx[i + 1];
      dr_coef1_fx[i] = 1.0 + 2.0 * coef_C * dom_beta[i] * big_bracket;
      dr_coef2_fx[i] = -1.0 - 2.0 * coef_C * for_beta[i] * big_bracket;
      dr_coef3_fx[i] = 2.0 * coef_C * big_bracket;
    }
  }
}

char *SrtGrfnFx3DFxBetaDLMTreeQBeta(char *und3dfx, int numeventdates,
                                    long *eventdates, long tableauRows,
                                    long tableauCols, char ***tableauStrings,
                                    int **tableauMask, long auxWidth,
                                    long *auxLen, double **aux,
                                    FxBetaDLM_GRFNNumerParams *Params,
                                    double *barrier, int *bar_col, long num_stp,
                                    int *num_prod, double **prod_val) {
  int free_str = 0;
  FIRSTAllMkts xStr;
  SrtGrfnParam defParm;
  int forback;
  long nstp;

  SrtUndPtr fx_und;

  double *time = NULL, *date = NULL;
  int *vol_change = NULL;

  double *sig_dom = NULL, *sig_for = NULL, *sig_fx = NULL, *corr_dom_for = NULL,
         *corr_dom_fx = NULL, *corr_for_fx = NULL;

  double *dr_const_dom = NULL, *dr_coef_dom = NULL, *dr_const_for = NULL,
         *dr_coef1_for = NULL, *dr_coef2_for = NULL, *dr_coef3_for = NULL,
         *dr_const_fx = NULL, *dr_coef1_fx = NULL, *dr_coef2_fx = NULL,
         *dr_coef3_fx = NULL;

  double *dom_fwd = NULL, *dom_var = NULL, *dom_phi = NULL, *dom_ifr = NULL,
         *dom_beta = NULL, *for_fwd = NULL, *for_var = NULL, *for_phi = NULL,
         *for_ifr = NULL, *for_beta = NULL, *fx_fwd = NULL, *fx_var = NULL,
         *var_ffx = NULL, *global_corr_dom_for = NULL,
         *global_corr_dom_fx = NULL, *global_corr_for_fx = NULL;

  void **void_prm = NULL;
  void **void_prm_down = NULL;
  GRFNPARMBETATREEDLM grfn_prm;
  int *is_event = NULL;

  long today;
  int i, j;
  double dom_lam, for_lam;
  double spot_fx;
  int fx_idx, dom_idx, for_idx;

  char *domname, *forname;

  int num_bar_times;
  double *bar_times = NULL;

  int *is_bar = NULL;
  double *bar_lvl = NULL;
  int *bar_cl = NULL;
  double next_bar;

  long next_d;

  /* For The Alpha Fudge */
  double *prod_val_temp = NULL;
  int alphaFudge = 0, UpOrDown;

  clock_t t1, t2;

  FxBetaDLM_model *model = NULL;

  Err err = NULL;

  t1 = clock();

  /*	Initialise the GRFN tableau */

  /*	First        , initialise the param struct */

  err = srt_f_set_default_GrfnParams(&defParm);
  defParm.min_nodes_per_path = num_stp;

  err = FIRSTInitMktStruct(numeventdates, eventdates, tableauRows, tableauCols,
                           tableauStrings, tableauMask, auxWidth, auxLen, aux,
                           und3dfx, &defParm, &forback, &xStr);

  if (err) {
    goto FREE_RETURN;
  }

  free_str = 1;

  /*	Now        , lookup underlyings involved and their term structures */
  model = calloc(1, sizeof(FxBetaDLM_model));

  if (!model) {
    err = "Memory allocation error in SrtGrfnFx3DBetaDLMTree";
    goto FREE_RETURN;
  }

  err = FxBetaDLM_Get_Model(und3dfx, model);

  if (err)
    goto FREE_RETURN;

  if (fabs(model->dAlpha) > 1e-6) {
    alphaFudge = 1;
  }

  today = model->lToday;
  dom_lam = model->dLambdaDom;
  for_lam = model->dLambdaFor;
  spot_fx = model->dSpotFx;

  fx_und = lookup_und(und3dfx);
  domname = get_domname_from_fxund(fx_und);
  forname = get_forname_from_fxund(fx_und);

  /*	Next        , get the time steps */

  /*	Copy event dates */
  nstp = xStr.num_evt;
  while (nstp >= 1 && xStr.evt[nstp - 1].evt == NULL) {
    nstp--;
  }
  if (nstp < 1) {
    err = "No event in Tableau";
    goto FREE_RETURN;
  }
  time = (double *)calloc(nstp, sizeof(double));
  if (!time) {
    err = "Memory allocation error (1) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }
  memcpy(time, xStr.tms, nstp * sizeof(double));

  /*	Get barrier dates */
  num_bar_times = 0;
  bar_times = (double *)calloc(nstp, sizeof(double));
  if (!bar_times) {
    err = "Memory allocation error (2) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }
  for (i = 0; i < xStr.num_evt; i++) {
    if (eventdates[i] > today && eventdates[i] <= xStr.dts[nstp - 1] &&
        barrier[i] > 1.0e-08) {
      bar_times[num_bar_times] = (eventdates[i] - today) * YEARS_IN_DAY;
      num_bar_times++;
    }
  }

  /*	Fill the time vector */

  err = fill_time_vector(&time, &nstp, num_bar_times, bar_times,
                         model->iNbPWTime, model->dPWTime, num_stp);

  if (err) {
    goto FREE_RETURN;
  }

  date = (double *)calloc(nstp, sizeof(double));
  if (!date) {
    err = "Memory allocation error (3) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }

  for (i = 0; i < nstp; i++) {
    date[i] = today + DAYS_IN_YEAR * time[i];

    if (i > 0 && date[i] - date[i - 1] >= 1) {
      date[i] = (long)(date[i] + 1.0e-08);
      time[i] = YEARS_IN_DAY * (date[i] - today);
    }
  }

  /*	Fill the model parameters as required by tree_main_3dfx */
  vol_change = (int *)calloc(nstp, sizeof(int));
  sig_dom = (double *)calloc(nstp, sizeof(double));
  sig_for = (double *)calloc(nstp, sizeof(double));
  sig_fx = (double *)calloc(nstp, sizeof(double));

  corr_dom_for = (double *)calloc(nstp, sizeof(double));
  corr_dom_fx = (double *)calloc(nstp, sizeof(double));
  corr_for_fx = (double *)calloc(nstp, sizeof(double));

  dr_const_dom = (double *)calloc(nstp, sizeof(double));
  dr_coef_dom = (double *)calloc(nstp, sizeof(double));
  dr_const_for = (double *)calloc(nstp, sizeof(double));
  dr_coef1_for = (double *)calloc(nstp, sizeof(double));
  dr_coef2_for = (double *)calloc(nstp, sizeof(double));
  dr_coef3_for = (double *)calloc(nstp, sizeof(double));
  dr_const_fx = (double *)calloc(nstp, sizeof(double));
  dr_coef1_fx = (double *)calloc(nstp, sizeof(double));
  dr_coef2_fx = (double *)calloc(nstp, sizeof(double));
  dr_coef3_fx = (double *)calloc(nstp, sizeof(double));

  /*	Get distributions */
  dom_ifr = (double *)calloc(nstp, sizeof(double));
  dom_fwd = (double *)calloc(nstp, sizeof(double));
  dom_var = (double *)calloc(nstp, sizeof(double));
  dom_phi = (double *)calloc(nstp, sizeof(double));
  dom_beta = (double *)calloc(nstp, sizeof(double));

  for_ifr = (double *)calloc(nstp, sizeof(double));
  for_fwd = (double *)calloc(nstp, sizeof(double));
  for_var = (double *)calloc(nstp, sizeof(double));
  for_phi = (double *)calloc(nstp, sizeof(double));
  for_beta = (double *)calloc(nstp, sizeof(double));

  fx_fwd = (double *)calloc(nstp, sizeof(double));
  fx_var = (double *)calloc(nstp, sizeof(double));
  var_ffx = (double *)calloc(nstp, sizeof(double));

  global_corr_dom_for = (double *)calloc(nstp, sizeof(double));
  global_corr_dom_fx = (double *)calloc(nstp, sizeof(double));
  global_corr_for_fx = (double *)calloc(nstp, sizeof(double));

  if (!vol_change || !sig_dom || !sig_for || !sig_fx || !corr_dom_for ||
      !corr_dom_fx || !corr_for_fx || !dr_const_dom || !dr_coef_dom ||
      !dr_const_for || !dr_coef1_for || !dr_coef2_for || !dr_coef3_for ||
      !dr_const_fx || !dr_coef1_fx || !dr_coef2_fx || !dr_coef3_fx ||
      !dom_ifr || !dom_fwd || !dom_var || !dom_phi || !dom_beta || !for_ifr ||
      !for_fwd || !for_var || !for_phi || !for_beta || !fx_fwd || !fx_var ||
      !var_ffx || !global_corr_dom_for || !global_corr_dom_fx ||
      !global_corr_for_fx) {
    err = "Memory allocation error (5) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }

  for (UpOrDown = 0; UpOrDown <= alphaFudge; UpOrDown++) {
    if (alphaFudge) {
      model->dAlpha = 0.0;

      if (!UpOrDown) {
        memcpy(model->dSigmaFx, model->dSigmaFxUp,
               model->iNbPWTime * sizeof(double));
      } else {
        memcpy(model->dSigmaFx, model->dSigmaFxDown,
               model->iNbPWTime * sizeof(double));
        memcpy(model->dCorrDomFx, model->dCorrDomFxDown,
               model->iNbPWTime * sizeof(double));
        memcpy(model->dCorrForFx, model->dCorrForFxDown,
               model->iNbPWTime * sizeof(double));
      }
    }

    /* Fill the global expectations and variances */
    Fx3DBetaDLM_PrecalcGRFNTreeQBeta(nstp, time, model, dom_fwd, dom_var,
                                     dom_phi, dom_beta, for_fwd, for_var,
                                     for_phi, for_beta, fx_fwd, fx_var, var_ffx,
                                     global_corr_dom_for, global_corr_dom_fx,
                                     global_corr_for_fx, Params->min_time);

    /* Fill the local drifts        , vols and correlations */
    Fx3DBetaDLM_DriftVolAndCorr(
        nstp, time, model, dom_phi, dom_beta, for_phi, for_beta, var_ffx,

        dr_const_dom, dr_coef_dom, dr_const_for, dr_coef1_for, dr_coef2_for,
        dr_coef3_for, dr_const_fx, dr_coef1_fx, dr_coef2_fx, dr_coef3_fx,

        vol_change, sig_dom, sig_for, sig_fx, corr_dom_for, corr_dom_fx,
        corr_for_fx);

    /*	Fill product structure */

    strupper(und3dfx);
    strip_white_space(und3dfx);
    strupper(domname);
    strip_white_space(domname);
    strupper(forname);
    strip_white_space(forname);

    for (i = 0; i < xStr.num_und; i++) {
      strupper(xStr.und_data[i].und_name);
      strip_white_space(xStr.und_data[i].und_name);
    }

    fx_idx = -1;
    for (i = 0; i < xStr.num_und; i++) {
      if (!strcmp(xStr.und_data[i].und_name, und3dfx)) {
        fx_idx = i;
      }
    }
    if (fx_idx == -1) {
      err = "The Fx underlying is not present in the mdlcomm structure";
      goto FREE_RETURN;
    }

    dom_idx = -1;
    for (i = 0; i < xStr.num_und; i++) {
      if (!strcmp(xStr.und_data[i].und_name, domname)) {
        dom_idx = i;
      }
    }
    if (dom_idx == -1) {
      err = "The domestic underlying is not present in the mdlcomm structure";
      goto FREE_RETURN;
    }

    for_idx = -1;
    for (i = 0; i < xStr.num_und; i++) {
      if (!strcmp(xStr.und_data[i].und_name, forname)) {
        for_idx = i;
      }
    }
    if (for_idx == -1) {
      err = "The foreign underlying is not present in the mdlcomm structure";
      goto FREE_RETURN;
    }

    if (!UpOrDown) {
      is_event = (int *)calloc(nstp, sizeof(int));
      void_prm = (void **)calloc(nstp, sizeof(void *));

      is_bar = (int *)calloc(nstp, sizeof(int));
      bar_lvl = (double *)calloc(nstp, sizeof(double));
      bar_cl = (int *)calloc(nstp, sizeof(int));

      if (!is_event || !void_prm || !is_bar || !bar_lvl || !bar_cl) {
        err = "Memory allocation error (6) in SrtGrfn3DFXTree";
        goto FREE_RETURN;
      }
    } else {
      void_prm_down = (void **)calloc(nstp, sizeof(void *));

      if (!void_prm_down) {
        err = "Memory allocation error (6) in SrtGrfn3DFXTree";
        goto FREE_RETURN;
      }
    }

    j = xStr.num_evt - 1;
    next_bar = 0.0;
    next_d = xStr.dts[j] + 1;

    for (i = nstp - 1; i >= 0; i--) {
      if (next_d > date[i]) {
        dom_ifr[i] = swp_f_zr(date[i], next_d, model->cYcDom);
      } else {
        dom_ifr[i] = swp_f_zr(date[i], date[i] + 1, model->cYcDom);
      }

      if (j >= 0 && fabs(date[i] - xStr.dts[j]) < 1.0e-04) {
        grfn_prm = malloc(sizeof(grfn_parm_beta_Tree_DLM));
        grfn_prm->global = &xStr;
        grfn_prm->local = xStr.evt + j;
        grfn_prm->fx_idx = fx_idx;
        grfn_prm->dom_idx = dom_idx;
        grfn_prm->for_idx = for_idx;

        grfn_prm->today = today;

        grfn_prm->lam_dom = dom_lam;
        grfn_prm->lam_for = for_lam;
        grfn_prm->beta_dom = dom_beta[i];
        grfn_prm->beta_for = for_beta[i];
        grfn_prm->phi_dom = dom_phi[i];
        grfn_prm->phi_for = for_phi[i];

        grfn_prm->B0 = model->dB0;
        grfn_prm->C0 = model->dC0;
        grfn_prm->varFFX = var_ffx[i];

        grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
        grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
        grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

        grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
        grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
        grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

        grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
        grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
        grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

        is_event[i] = 1;
        if (!UpOrDown)
          void_prm[i] = (void *)grfn_prm;
        else
          void_prm_down[i] = (void *)grfn_prm;

        if (barrier[j] > 1.0e-08) {
          next_bar = log(barrier[j] / spot_fx);
          if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols - 1) {
            is_bar[i] = 1;
            bar_cl[i] = bar_col[j];
          } else {
            is_bar[i] = 0;
          }
        } else {
          is_bar[i] = 0;
        }

        j--;
        while (j >= 0 && xStr.evt[j].evt == NULL) {
          j--;
        }
      } else if (j >= 0 && xStr.am[j]) {
        grfn_prm = malloc(sizeof(grfn_parm_tree));
        grfn_prm->global = &xStr;
        grfn_prm->local = xStr.evt + j;
        grfn_prm->fx_idx = fx_idx;
        grfn_prm->dom_idx = dom_idx;
        grfn_prm->for_idx = for_idx;

        grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
        grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
        grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

        grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
        grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
        grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

        grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
        grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
        grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

        if (barrier[j] > 1.0e-08) {
          next_bar = log(barrier[j] / spot_fx);
          if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols - 1) {
            is_bar[i] = 1;
            bar_cl[i] = bar_col[j];
          } else {
            is_bar[i] = 0;
          }
        } else {
          is_bar[i] = 0;
        }

        is_event[i] = 1;
        if (!UpOrDown)
          void_prm[i] = (void *)grfn_prm;
        else
          void_prm_down[i] = (void *)grfn_prm;
      } else {
        is_event[i] = 0;
        if (!UpOrDown)
          void_prm[i] = NULL;
        else
          void_prm_down[i] = NULL;
        is_bar[i] = 0;
      }

      if (fabs(next_bar) > 1.0e-08) {
        bar_lvl[i] = next_bar;
      } else {
        bar_lvl[i] = fx_fwd[i];
      }

      if ((i > 0) && (fabs(next_bar) > 1.0e-08)) {
        next_bar += (for_ifr[i - 1] + for_fwd[i - 1] - dom_ifr[i - 1] -
                     dom_fwd[i - 1]) *
                    (time[i] - time[i - 1]);
      }

      next_d = (long)(date[i] + 0.5);
    }

    is_bar[0] = 0;
    bar_lvl[0] = fx_fwd[0];
    bar_cl[0] = -1;

    is_bar[1] = 0;
    bar_lvl[1] = fx_fwd[1];
    bar_cl[1] = -1;

    /*	Eventually! call to function */

    *num_prod = xStr.num_cols;
    if (!UpOrDown) {
      *prod_val = (double *)calloc(*num_prod + 1, sizeof(double));
      prod_val_temp = (double *)calloc(*num_prod + 1, sizeof(double));

      if (!(*prod_val) || !prod_val_temp) {
        err = "Memory allocation error (7) in SrtGrfn3DFXTree";
        goto FREE_RETURN;
      }
    }

    t2 = clock();

    smessage("Phase 1 -preprocessing        , time in sec: %.2f",
             (double)(t2 - t1) / CLOCKS_PER_SEC);
    smessage("Number of times steps        , required: %d        , actual: %d",
             num_stp, nstp);

    err = tree_main_3dBetaDLM_QBeta(
        nstp, time, date,

        vol_change, dom_lam, for_lam, sig_dom, sig_for, sig_fx, corr_dom_for,
        corr_dom_fx, corr_for_fx, dr_const_dom, dr_coef_dom, dr_const_for,
        dr_coef1_for, dr_coef2_for, dr_coef3_for, dr_const_fx, dr_coef1_fx,
        dr_coef2_fx, dr_coef3_fx,

        dom_fwd, dom_var, for_fwd, for_var, fx_fwd, fx_var, corr_dom_for,
        corr_dom_fx, corr_for_fx,

        model->dCashFx, model->cYcDom, model->cYcFor, dom_ifr,

        UpOrDown ? void_prm_down : void_prm, is_event, bar_lvl, bar_cl, is_bar,
        grfn_payoff_4_3dfxBetaDLMQBeta_Tree, *num_prod, Params->do_discount,
        UpOrDown ? prod_val_temp : (*prod_val));
  }

  if (alphaFudge) {
    for (i = 0; i < (*num_prod + 1); i++) {
      (*prod_val)[i] = 0.5 * ((*prod_val)[i] + prod_val_temp[i]);
    }
  }
  (*prod_val)[*num_prod - 1] += xStr.gd->pv_of_past;

FREE_RETURN:

  if (time)
    free(time);
  if (date)
    free(date);

  if (vol_change)
    free(vol_change);
  if (sig_dom)
    free(sig_dom);
  if (sig_for)
    free(sig_for);
  if (sig_fx)
    free(sig_fx);
  if (corr_dom_for)
    free(corr_dom_for);
  if (corr_dom_fx)
    free(corr_dom_fx);
  if (corr_for_fx)
    free(corr_for_fx);

  if (dr_const_dom)
    free(dr_const_dom);
  if (dr_coef_dom)
    free(dr_coef_dom);
  if (dr_const_for)
    free(dr_const_for);
  if (dr_coef1_for)
    free(dr_coef1_for);
  if (dr_coef2_for)
    free(dr_coef2_for);
  if (dr_coef3_for)
    free(dr_coef3_for);
  if (dr_const_fx)
    free(dr_const_fx);
  if (dr_coef1_fx)
    free(dr_coef1_fx);
  if (dr_coef2_fx)
    free(dr_coef2_fx);
  if (dr_coef3_fx)
    free(dr_coef3_fx);

  if (dom_ifr)
    free(dom_ifr);
  if (dom_fwd)
    free(dom_fwd);
  if (dom_var)
    free(dom_var);
  if (dom_phi)
    free(dom_phi);
  if (dom_beta)
    free(dom_beta);

  if (for_ifr)
    free(for_ifr);
  if (for_fwd)
    free(for_fwd);
  if (for_var)
    free(for_var);
  if (for_phi)
    free(for_phi);
  if (for_beta)
    free(for_beta);

  if (fx_fwd)
    free(fx_fwd);
  if (fx_var)
    free(fx_var);
  if (var_ffx)
    free(var_ffx);

  if (global_corr_dom_for)
    free(global_corr_dom_for);
  if (global_corr_dom_fx)
    free(global_corr_dom_fx);
  if (global_corr_for_fx)
    free(global_corr_for_fx);

  if (bar_times)
    free(bar_times);

  if (is_event)
    free(is_event);

  if (bar_lvl)
    free(bar_lvl);
  if (is_bar)
    free(is_bar);
  if (bar_cl)
    free(bar_cl);

  if (void_prm) {
    for (i = 0; i < nstp; i++) {
      if (void_prm[i]) {
        grfn_prm = (GRFNPARMBETATREEDLM)void_prm[i];

        free(grfn_prm);
      }
    }

    free(void_prm);
  }

  if (void_prm_down) {
    for (i = 0; i < nstp; i++) {
      if (void_prm_down[i]) {
        grfn_prm = (GRFNPARMBETATREEDLM)void_prm_down[i];

        free(grfn_prm);
      }
    }

    free(void_prm_down);
  }

  if (free_str) {
    FIRSTFreeMktStruct(&xStr);
  }

  if (model) {
    free_FxBetaDLM_model(model);
    free(model);
  }

  if (prod_val_temp)
    free(prod_val_temp);

  return err;
}

char *SrtGrfnFx3DFxBetaDLMTreeQTStar(
    char *und3dfx, int numeventdates, long *eventdates, long tableauRows,
    long tableauCols, char ***tableauStrings, int **tableauMask, long auxWidth,
    long *auxLen, double **aux, FxBetaDLM_GRFNNumerParams *Params,
    double *barrier, int *bar_col, long num_stp, int *num_prod,
    double **prod_val) {
  int free_str = 0;
  FIRSTAllMkts xStr;
  SrtGrfnParam defParm;
  int forback;
  long nstp;

  SrtUndPtr fx_und;

  double *time = NULL, *date = NULL;
  int *vol_change = NULL;

  double *sig_dom = NULL, *mu_quanto_cons = NULL, *mu_quanto_lin = NULL,
         *sig_for = NULL, *sig_fx = NULL, *corr_dom_for = NULL,
         *corr_dom_fx = NULL, *corr_for_fx = NULL;

  double *dom_fwd = NULL, *dom_var = NULL, *dom_phi = NULL, *dom_beta = NULL,
         *for_fwd = NULL, *for_var = NULL, *for_phi = NULL, *for_beta = NULL,
         *fx_fwd = NULL, *fx_var = NULL, *ffx_var = NULL,
         *global_corr_dom_for = NULL, *global_corr_dom_fx = NULL,
         *global_corr_for_fx = NULL;

  void **void_prm = NULL;
  void **void_prm_down = NULL;
  GRFNPARMBETAMCDLM grfn_prm;
  int *is_event = NULL;

  long today;
  int i, j, k, index;
  double dom_lam, for_lam;
  double spot_fx;
  int fx_idx, dom_idx, for_idx;

  char *domname, *forname;

  int num_bar_times;
  double *bar_times = NULL;

  int *is_bar = NULL;
  double *bar_lvl = NULL;
  int *bar_cl = NULL;
  double next_bar;

  long pay_date;
  double df, var_ffx, temp, T, beta_temp;
  double t_mid;
  double temp_beta_dom, temp_beta_for, temp_vol_dom, temp_vol_for, temp_vol_fx;

  clock_t t1, t2;

  /* For The Alpha Fudge */
  double *prod_val_temp = NULL;
  int alphaFudge = 0, UpOrDown;

  FxBetaDLM_model *model = NULL;

  Err err = NULL;

  t1 = clock();

  /*	Initialise the GRFN tableau */

  /*	First        , initialise the param struct */

  err = srt_f_set_default_GrfnParams(&defParm);
  defParm.min_nodes_per_path = num_stp;

  err = FIRSTInitMktStruct(numeventdates, eventdates, tableauRows, tableauCols,
                           tableauStrings, tableauMask, auxWidth, auxLen, aux,
                           und3dfx, &defParm, &forback, &xStr);

  if (err) {
    goto FREE_RETURN;
  }

  free_str = 1;

  /*	Now        , lookup underlyings involved and their term structures */
  model = calloc(1, sizeof(FxBetaDLM_model));

  if (!model) {
    err = "Memory allocation error in SrtGrfnFx3DBetaDLMTree";
    goto FREE_RETURN;
  }

  err = FxBetaDLM_Get_Model(und3dfx, model);
  if (err)
    goto FREE_RETURN;
  if (fabs(model->dAlpha) > 1e-6) {
    alphaFudge = 1;
  }

  today = model->lToday;
  dom_lam = model->dLambdaDom;
  for_lam = model->dLambdaFor;
  spot_fx = model->dSpotFx;

  fx_und = lookup_und(und3dfx);
  domname = get_domname_from_fxund(fx_und);
  forname = get_forname_from_fxund(fx_und);

  /*	Next        , get the time steps */

  /*	Copy event dates */
  nstp = xStr.num_evt;
  while (nstp >= 1 && xStr.evt[nstp - 1].evt == NULL) {
    nstp--;
  }
  if (nstp < 1) {
    err = "No event in Tableau";
    goto FREE_RETURN;
  }
  time = (double *)calloc(nstp, sizeof(double));
  if (!time) {
    err = "Memory allocation error (1) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }
  memcpy(time, xStr.tms, nstp * sizeof(double));

  /*	Get barrier dates */
  num_bar_times = 0;
  bar_times = (double *)calloc(nstp, sizeof(double));
  if (!bar_times) {
    err = "Memory allocation error (2) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }
  for (i = 0; i < xStr.num_evt; i++) {
    if (eventdates[i] > today && eventdates[i] <= xStr.dts[nstp - 1] &&
        barrier[i] > 1.0e-08) {
      bar_times[num_bar_times] = (eventdates[i] - today) * YEARS_IN_DAY;
      num_bar_times++;
    }
  }

  /*	Fill the time vector */

  err = fill_time_vector(&time, &nstp, num_bar_times, bar_times,
                         model->iNbPWTime, model->dPWTime, num_stp);

  if (err) {
    goto FREE_RETURN;
  }

  date = (double *)calloc(nstp, sizeof(double));
  if (!date) {
    err = "Memory allocation error (3) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }

  for (i = 0; i < nstp; i++) {
    date[i] = today + DAYS_IN_YEAR * time[i];

    if (i > 0 && date[i] - date[i - 1] >= 1) {
      date[i] = (long)(date[i] + 1.0e-08);
      time[i] = YEARS_IN_DAY * (date[i] - today);
    }
  }

  /*	Fill the model parameters as required by tree_main_3dfx */
  vol_change = (int *)calloc(nstp, sizeof(int));
  sig_dom = (double *)calloc(nstp, sizeof(double));
  sig_for = (double *)calloc(nstp, sizeof(double));
  mu_quanto_cons = (double *)calloc(nstp, sizeof(double));
  mu_quanto_lin = (double *)calloc(nstp, sizeof(double));
  sig_fx = (double *)calloc(nstp, sizeof(double));

  corr_dom_for = (double *)calloc(nstp, sizeof(double));
  corr_dom_fx = (double *)calloc(nstp, sizeof(double));
  corr_for_fx = (double *)calloc(nstp, sizeof(double));

  /*	Get distributions */
  dom_fwd = (double *)calloc(nstp, sizeof(double));
  dom_var = (double *)calloc(nstp, sizeof(double));
  dom_phi = (double *)calloc(nstp, sizeof(double));
  dom_beta = (double *)calloc(nstp, sizeof(double));

  for_fwd = (double *)calloc(nstp, sizeof(double));
  for_var = (double *)calloc(nstp, sizeof(double));
  for_phi = (double *)calloc(nstp, sizeof(double));
  for_beta = (double *)calloc(nstp, sizeof(double));

  fx_fwd = (double *)calloc(nstp, sizeof(double));
  fx_var = (double *)calloc(nstp, sizeof(double));
  ffx_var = (double *)calloc(nstp, sizeof(double));

  global_corr_dom_for = (double *)calloc(nstp, sizeof(double));
  global_corr_dom_fx = (double *)calloc(nstp, sizeof(double));
  global_corr_for_fx = (double *)calloc(nstp, sizeof(double));

  if (!vol_change || !sig_dom || !mu_quanto_cons || !mu_quanto_lin ||
      !sig_for || !sig_fx || !corr_dom_for || !corr_dom_fx || !corr_for_fx ||
      !dom_fwd || !dom_var || !dom_phi || !dom_beta || !for_fwd || !for_var ||
      !for_phi || !for_beta || !fx_fwd || !fx_var || !ffx_var ||
      !global_corr_dom_for || !global_corr_dom_fx || !global_corr_for_fx) {
    err = "Memory allocation error (5) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }

  for (UpOrDown = 0; UpOrDown <= alphaFudge; UpOrDown++) {

    if (alphaFudge) {
      model->dAlpha = 0.0;

      if (!UpOrDown) {
        memcpy(model->dSigmaFx, model->dSigmaFxUp,
               model->iNbPWTime * sizeof(double));
      } else {
        memcpy(model->dSigmaFx, model->dSigmaFxDown,
               model->iNbPWTime * sizeof(double));
        memcpy(model->dCorrDomFx, model->dCorrDomFxDown,
               model->iNbPWTime * sizeof(double));
        memcpy(model->dCorrForFx, model->dCorrForFxDown,
               model->iNbPWTime * sizeof(double));
      }
    }

    /* Fill the global expectations and variances */
    Fx3DBetaDLM_ExpectAndVarGrfn(nstp, time, model, dom_fwd, dom_var, dom_phi,
                                 dom_beta, for_fwd, NULL, for_var, for_phi,
                                 for_beta, fx_fwd, fx_var, ffx_var, 1,
                                 global_corr_dom_for, global_corr_dom_fx,
                                 global_corr_for_fx, Params->min_time);

    /* Fill the local vols */
    vol_change[nstp - 1] = 1;
    index = Get_Index(time[nstp - 1], model->dPWTime, model->iNbPWTime);

    for (i = nstp - 1; i >= 0; i--) {
      if (i < nstp - 1) {
        t_mid = 0.5 * (time[i] + time[i + 1]);
      } else {
        t_mid = 0.5 * (time[i] + time[i - 1]);
      }

      index = Get_Index(t_mid, model->dPWTime, model->iNbPWTime);

      temp_beta_dom = exp(-dom_lam * (model->dTStar - t_mid));
      temp_beta_for = exp(-for_lam * (model->dTStar - t_mid));
      temp_vol_dom = model->dSigmaDom[index] * (1.0 - temp_beta_dom) / dom_lam;
      temp_vol_for = model->dSigmaFor[index] * (1.0 - temp_beta_for) / for_lam;
      temp_vol_fx = model->dSigmaFx[index];

      sig_dom[i] = model->dSigmaDom[index] * temp_beta_dom;
      sig_for[i] = model->dSigmaFor[index] * temp_beta_for;

      sig_fx[i] =
          pow(temp_vol_fx, 2) + pow(temp_vol_dom, 2) + pow(temp_vol_for, 2) -
          2.0 * model->dCorrDomFor[index] * temp_vol_dom * temp_vol_for +
          2.0 * model->dCorrDomFx[index] * temp_vol_dom * temp_vol_fx -
          2.0 * model->dCorrForFx[index] * temp_vol_for * temp_vol_fx;

      if (sig_fx[i] < 0.0) {
        err = "Correlation matrix is not positive definite in "
              "SrtgrfnMainFx3fBetaDLM";
        goto FREE_RETURN;
      }

      sig_fx[i] = sqrt(sig_fx[i]);

      corr_dom_for[i] = model->dCorrDomFor[index];

      corr_dom_fx[i] =
          (temp_vol_dom - model->dCorrDomFor[index] * temp_vol_for +
           model->dCorrDomFx[index] * temp_vol_fx) /
          sig_fx[i];

      corr_for_fx[i] = (model->dCorrDomFor[index] * temp_vol_dom -
                        temp_vol_for + model->dCorrForFx[index] * temp_vol_fx) /
                       sig_fx[i];

      if (i < nstp - 1) {
        mu_quanto_cons[i] =
            -corr_for_fx[i] * sig_for[i] * sig_fx[i] /
            (1.0 + 2.0 * model->dC0 * 0.5 * (fx_var[i] + fx_var[i + 1]));
      } else {
        mu_quanto_cons[i] =
            -corr_for_fx[i] * sig_for[i] * sig_fx[i] /
            (1.0 + 2.0 * model->dC0 * 0.5 * (fx_var[i] + fx_var[i - 1]));
      }

      mu_quanto_lin[i] = 2.0 * model->dC0 * mu_quanto_cons[i];

      /* calculation of the real volatility and correlations */

      if ((i == nstp - 1) || (fabs(sig_dom[i] - sig_dom[i + 1]) +
                                  fabs(sig_for[i] - sig_for[i + 1]) +
                                  fabs(sig_fx[i] - sig_fx[i + 1]) +
                                  fabs(corr_dom_for[i] - corr_dom_for[i + 1]) +
                                  fabs(corr_dom_fx[i] - corr_dom_fx[i + 1]) +
                                  fabs(corr_for_fx[i] - corr_for_fx[i + 1]) >
                              EPS)) {
        vol_change[i] = 1;
      } else {
        vol_change[i] = 0;
      }
    }

    /*	Fill product structure */

    strupper(und3dfx);
    strip_white_space(und3dfx);
    strupper(domname);
    strip_white_space(domname);
    strupper(forname);
    strip_white_space(forname);

    for (i = 0; i < xStr.num_und; i++) {
      strupper(xStr.und_data[i].und_name);
      strip_white_space(xStr.und_data[i].und_name);
    }

    fx_idx = -1;
    for (i = 0; i < xStr.num_und; i++) {
      if (!strcmp(xStr.und_data[i].und_name, und3dfx)) {
        fx_idx = i;
      }
    }
    if (fx_idx == -1) {
      err = "The Fx underlying is not present in the mdlcomm structure";
      goto FREE_RETURN;
    }

    dom_idx = -1;
    for (i = 0; i < xStr.num_und; i++) {
      if (!strcmp(xStr.und_data[i].und_name, domname)) {
        dom_idx = i;
      }
    }
    if (dom_idx == -1) {
      err = "The domestic underlying is not present in the mdlcomm structure";
      goto FREE_RETURN;
    }

    for_idx = -1;
    for (i = 0; i < xStr.num_und; i++) {
      if (!strcmp(xStr.und_data[i].und_name, forname)) {
        for_idx = i;
      }
    }
    if (for_idx == -1) {
      err = "The foreign underlying is not present in the mdlcomm structure";
      goto FREE_RETURN;
    }

    if (!UpOrDown) {
      is_event = (int *)calloc(nstp, sizeof(int));
      void_prm = (void **)calloc(nstp, sizeof(void *));

      is_bar = (int *)calloc(nstp, sizeof(int));
      bar_lvl = (double *)calloc(nstp, sizeof(double));
      bar_cl = (int *)calloc(nstp, sizeof(int));

      if (!is_event || !void_prm || !is_bar || !bar_lvl || !bar_cl) {
        err = "Memory allocation error (6) in SrtGrfn3DFXTree";
        goto FREE_RETURN;
      }
    } else {
      void_prm_down = (void **)calloc(nstp, sizeof(void *));

      if (!void_prm_down) {
        err = "Memory allocation error (6) in SrtGrfn3DFXTree";
        goto FREE_RETURN;
      }
    }

    pay_date = model->lTStarDate;
    df = swp_f_df(today, pay_date, model->cYcDom);

    j = xStr.num_evt - 1;
    next_bar = 0.0;
    for (i = nstp - 1; i >= 0; i--) {
      if (j >= 0 && fabs(date[i] - xStr.dts[j]) < 1.0e-04) {
        grfn_prm = malloc(sizeof(grfn_parm_beta_MC_DLM));

        if (!grfn_prm) {
          err = "Memory allocation error (7) in SrtGrfn3DFXMc";
          goto FREE_RETURN;
        }

        grfn_prm->global = &xStr;
        grfn_prm->local = xStr.evt + j;
        grfn_prm->fx_idx = fx_idx;
        grfn_prm->dom_idx = dom_idx;
        grfn_prm->for_idx = for_idx;

        /* Pay Bond reconstruction */
        grfn_prm->bond_pay_const =
            log(swp_f_df(xStr.dts[j], pay_date, model->cYcDom));
        grfn_prm->bond_pay_const +=
            0.5 * dom_beta[i] * dom_beta[i] * dom_phi[i];
        grfn_prm->bond_pay_lin = dom_beta[i];

        /* Spot Fx reconstruction */
        var_ffx = fx_var[i];
        temp = 1.0 + 2.0 * model->dC0 * var_ffx;

        grfn_prm->fx_coef_const =
            log(model->dCashFx * swp_f_df(today, xStr.dts[j], model->cYcFor) /
                swp_f_df(today, xStr.dts[j], model->cYcDom));

        grfn_prm->fx_coef_const +=
            -0.5 * (model->dB0 * model->dB0 * var_ffx / temp + log(temp));

        grfn_prm->fx_coef_const +=
            0.5 * (dom_beta[i] * dom_beta[i] * dom_phi[i] -
                   for_beta[i] * for_beta[i] * for_phi[i]);

        grfn_prm->fx_coef_lin = model->dB0 / temp;
        grfn_prm->fx_coef_quad = model->dC0 / temp;
        grfn_prm->fx_coef_dom = dom_beta[i];
        grfn_prm->fx_coef_for = -for_beta[i];

        /* Fx discount factors */
        grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
        grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
        grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

        if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1) {
          grfn_prm->df_fx_coef_const = dvector(0, grfn_prm->num_fx_df - 1);
          grfn_prm->df_fx_coef_lin = dvector(0, grfn_prm->num_fx_df - 1);

          if (!grfn_prm->df_fx_coef_const || !grfn_prm->df_fx_coef_lin) {
            err = "Memory allocation error (7) in SrtGrfn3DFXMc";
            goto FREE_RETURN;
          }

          for (k = 0; k < grfn_prm->num_fx_df; k++) {
            T = xStr.tms[j] + grfn_prm->fx_df_tms[k];
            beta_temp = (1.0 - exp(-model->dLambdaDom * (T - model->dTStar))) /
                        model->dLambdaDom;
            grfn_prm->df_fx_coef_const[k] = log(
                swp_f_df(xStr.dts[j], grfn_prm->fx_df_dts[k], model->cYcDom));
            grfn_prm->df_fx_coef_const[k] -=
                0.5 * (beta_temp * beta_temp - dom_beta[i] * dom_beta[i]) *
                dom_phi[i];
            grfn_prm->df_fx_coef_lin[k] = -(beta_temp - dom_beta[i]);
          }

          grfn_prm->do_fx = 1;
        } else {
          grfn_prm->do_fx = 0;
        }

        /* Domestic discount factors */
        grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
        grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
        grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

        if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1) {
          grfn_prm->df_dom_coef_const = dvector(0, grfn_prm->num_dom_df - 1);
          grfn_prm->df_dom_coef_lin = dvector(0, grfn_prm->num_dom_df - 1);

          if (!grfn_prm->df_dom_coef_const || !grfn_prm->df_dom_coef_lin) {
            err = "Memory allocation error (7) in SrtGrfn3DFXMc";
            goto FREE_RETURN;
          }

          for (k = 0; k < grfn_prm->num_dom_df; k++) {
            T = xStr.tms[j] + grfn_prm->dom_df_tms[k];
            beta_temp = (1.0 - exp(-model->dLambdaDom * (T - model->dTStar))) /
                        model->dLambdaDom;
            grfn_prm->df_dom_coef_const[k] = log(
                swp_f_df(xStr.dts[j], grfn_prm->dom_df_dts[k], model->cYcDom));
            grfn_prm->df_dom_coef_const[k] -=
                0.5 * (beta_temp * beta_temp - dom_beta[i] * dom_beta[i]) *
                dom_phi[i];
            grfn_prm->df_dom_coef_lin[k] = -(beta_temp - dom_beta[i]);
          }

          grfn_prm->do_dom = 1;
        } else {
          grfn_prm->do_dom = 0;
        }

        /* Foreign discount factors */
        grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
        grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
        grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

        if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1) {
          grfn_prm->df_for_coef_const = dvector(0, grfn_prm->num_for_df - 1);
          grfn_prm->df_for_coef_lin = dvector(0, grfn_prm->num_for_df - 1);

          if (!grfn_prm->df_for_coef_const || !grfn_prm->df_for_coef_lin) {
            err = "Memory allocation error (7) in SrtGrfn3DFXMc";
            goto FREE_RETURN;
          }

          for (k = 0; k < grfn_prm->num_for_df; k++) {
            T = xStr.tms[j] + grfn_prm->for_df_tms[k];
            beta_temp = (1.0 - exp(-model->dLambdaFor * (T - model->dTStar))) /
                        model->dLambdaFor;
            grfn_prm->df_for_coef_const[k] = log(
                swp_f_df(xStr.dts[j], grfn_prm->for_df_dts[k], model->cYcFor));
            grfn_prm->df_for_coef_const[k] -=
                0.5 * (beta_temp * beta_temp - for_beta[i] * for_beta[i]) *
                for_phi[i];
            grfn_prm->df_for_coef_lin[k] = -(beta_temp - for_beta[i]);
          }

          grfn_prm->do_for = 1;
        } else {
          grfn_prm->do_for = 0;
        }

        is_event[i] = 1;
        void_prm[i] = (void *)grfn_prm;

        if (barrier[j] > 1.0e-08) {
          next_bar = log(barrier[j] / spot_fx);
          if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols - 1) {
            is_bar[i] = 1;
            bar_cl[i] = bar_col[j];
          } else {
            is_bar[i] = 0;
          }
        } else {
          is_bar[i] = 0;
        }

        j--;
        while (j >= 0 && xStr.evt[j].evt == NULL) {
          j--;
        }
      } else if (j >= 0 && xStr.am[j]) {
        grfn_prm = malloc(sizeof(grfn_parm_beta_MC_DLM));

        if (!grfn_prm) {
          err = "Memory allocation error (7) in SrtGrfn3DFXMc";
          goto FREE_RETURN;
        }

        grfn_prm->global = &xStr;
        grfn_prm->local = xStr.evt + j;
        grfn_prm->fx_idx = fx_idx;
        grfn_prm->dom_idx = dom_idx;
        grfn_prm->for_idx = for_idx;

        /* Pay Bond reconstruction */
        grfn_prm->bond_pay_const =
            log(swp_f_df(xStr.dts[j], pay_date, model->cYcDom));
        grfn_prm->bond_pay_const +=
            0.5 * dom_beta[i] * dom_beta[i] * dom_phi[i];
        grfn_prm->bond_pay_lin = dom_beta[i];

        /* Spot Fx reconstruction */
        var_ffx = fx_var[i] * fx_var[i];
        temp = 1.0 + 2.0 * model->dC0 * var_ffx;

        grfn_prm->fx_coef_const =
            log(model->dCashFx * swp_f_df(today, xStr.dts[j], model->cYcFor) /
                swp_f_df(today, xStr.dts[j], model->cYcDom));

        grfn_prm->fx_coef_const +=
            -0.5 * (model->dB0 * model->dB0 * var_ffx / temp + log(temp));

        grfn_prm->fx_coef_const +=
            0.5 * (dom_beta[i] * dom_beta[i] * dom_phi[i] -
                   for_beta[i] * for_beta[i] * for_phi[i]);

        grfn_prm->fx_coef_lin = model->dB0 / temp;
        grfn_prm->fx_coef_quad = model->dC0 / temp;
        grfn_prm->fx_coef_dom = dom_beta[i];
        grfn_prm->fx_coef_for = -for_beta[i];

        /* Fx discount factors */
        grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
        grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
        grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

        if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1) {
          grfn_prm->df_fx_coef_const = dvector(0, grfn_prm->num_fx_df - 1);
          grfn_prm->df_fx_coef_lin = dvector(0, grfn_prm->num_fx_df - 1);

          if (!grfn_prm->df_fx_coef_const || !grfn_prm->df_fx_coef_lin) {
            err = "Memory allocation error (7) in SrtGrfn3DFXMc";
            goto FREE_RETURN;
          }

          for (k = 0; k < grfn_prm->num_fx_df; k++) {
            T = xStr.tms[j] + grfn_prm->fx_df_tms[k];
            beta_temp = (1.0 - exp(-model->dLambdaDom * (T - model->dTStar))) /
                        model->dLambdaDom;
            grfn_prm->df_fx_coef_const[k] = log(
                swp_f_df(xStr.dts[j], grfn_prm->fx_df_dts[k], model->cYcDom));
            grfn_prm->df_fx_coef_const[k] -=
                0.5 * (beta_temp * beta_temp - dom_beta[i] * dom_beta[i]) *
                dom_phi[i];
            grfn_prm->df_fx_coef_lin[k] = -(beta_temp - dom_beta[i]);
          }

          grfn_prm->do_fx = 1;
        } else {
          grfn_prm->do_fx = 0;
        }

        /* Domestic discount factors */
        grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
        grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
        grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

        if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1) {
          grfn_prm->df_dom_coef_const = dvector(0, grfn_prm->num_dom_df - 1);
          grfn_prm->df_dom_coef_lin = dvector(0, grfn_prm->num_dom_df - 1);

          if (!grfn_prm->df_dom_coef_const || !grfn_prm->df_dom_coef_lin) {
            err = "Memory allocation error (7) in SrtGrfn3DFXMc";
            goto FREE_RETURN;
          }

          for (k = 0; k < grfn_prm->num_dom_df; k++) {
            T = xStr.tms[j] + grfn_prm->dom_df_tms[k];
            beta_temp = (1.0 - exp(-model->dLambdaDom * (T - model->dTStar))) /
                        model->dLambdaDom;
            grfn_prm->df_dom_coef_const[k] = log(
                swp_f_df(xStr.dts[j], grfn_prm->dom_df_dts[k], model->cYcDom));
            grfn_prm->df_dom_coef_const[k] -=
                0.5 * (beta_temp * beta_temp - dom_beta[i] * dom_beta[i]) *
                dom_phi[i];
            grfn_prm->df_dom_coef_lin[k] = -(beta_temp - dom_beta[i]);
          }

          grfn_prm->do_dom = 1;
        } else {
          grfn_prm->do_dom = 0;
        }

        /* Foreign discount factors */
        grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
        grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
        grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

        if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1) {
          grfn_prm->df_for_coef_const = dvector(0, grfn_prm->num_for_df - 1);
          grfn_prm->df_for_coef_lin = dvector(0, grfn_prm->num_for_df - 1);

          if (!grfn_prm->df_for_coef_const || !grfn_prm->df_for_coef_lin) {
            err = "Memory allocation error (7) in SrtGrfn3DFXMc";
            goto FREE_RETURN;
          }

          for (k = 0; k < grfn_prm->num_for_df; k++) {
            T = xStr.tms[j] + grfn_prm->for_df_tms[k];
            beta_temp = (1.0 - exp(-model->dLambdaFor * (T - model->dTStar))) /
                        model->dLambdaFor;
            grfn_prm->df_for_coef_const[k] = log(
                swp_f_df(xStr.dts[j], grfn_prm->for_df_dts[k], model->cYcFor));
            grfn_prm->df_for_coef_const[k] -=
                0.5 * (beta_temp * beta_temp - for_beta[i] * for_beta[i]) *
                for_phi[i];
            grfn_prm->df_for_coef_lin[k] = -(beta_temp - for_beta[i]);
          }

          grfn_prm->do_for = 1;
        } else {
          grfn_prm->do_for = 0;
        }

        if (barrier[j] > 1.0e-08) {
          next_bar = log(barrier[j] / spot_fx);
          if (bar_col[j] >= 0 && bar_col[j] <= xStr.num_cols - 1) {
            is_bar[i] = 1;
            bar_cl[i] = bar_col[j];
          } else {
            is_bar[i] = 0;
          }
        } else {
          is_bar[i] = 0;
        }

        is_event[i] = 1;
        void_prm[i] = (void *)grfn_prm;
      } else {
        is_event[i] = 0;
        void_prm[i] = NULL;
        is_bar[i] = 0;
      }

      if (fabs(next_bar) > 1.0e-08) {
        bar_lvl[i] = next_bar;
      } else {
        bar_lvl[i] = fx_fwd[i];
      }

      if ((i > 0) && (fabs(next_bar) > 1.0e-08)) {
        next_bar +=
            (for_fwd[i - 1] - -dom_fwd[i - 1]) * (time[i] - time[i - 1]);
      }
    }

    is_bar[0] = 0;
    bar_lvl[0] = fx_fwd[0];
    bar_cl[0] = -1;

    is_bar[1] = 0;
    bar_lvl[1] = fx_fwd[1];
    bar_cl[1] = -1;

    /*	Eventually! call to function */

    *num_prod = xStr.num_cols;
    if (!UpOrDown) {
      *prod_val = (double *)calloc(*num_prod + 1, sizeof(double));
      if (alphaFudge) {
        prod_val_temp = (double *)calloc(*num_prod + 1, sizeof(double));
      }
    }
    if (!(*prod_val) || !prod_val_temp) {
      err = "Memory allocation error (7) in SrtGrfn3DFXTree";
      goto FREE_RETURN;
    }

    t2 = clock();

    smessage("Phase 1 -preprocessing        , time in sec: %.2f",
             (double)(t2 - t1) / CLOCKS_PER_SEC);
    smessage("Number of times steps        , required: %d        , actual: %d",
             num_stp, nstp);

    err = tree_main_3dBetaDLM_QTStar(
        nstp, time, date, vol_change, sig_dom, mu_quanto_cons, mu_quanto_lin,
        sig_for, sig_fx, dom_fwd, dom_var, for_fwd, for_var, fx_fwd, fx_var,
        dom_lam, for_lam, corr_dom_for, corr_dom_fx, corr_for_fx,
        global_corr_dom_for, global_corr_dom_fx, global_corr_for_fx,
        UpOrDown ? void_prm_down : void_prm, is_event, bar_lvl, bar_cl, is_bar,
        grfn_payoff_4_3dfxBetaDLMQTStar_Tree, *num_prod,
        UpOrDown ? prod_val_temp : (*prod_val));
  }
  /*	Add PV of Past */
  for (i = 0; i < *num_prod; i++) {
    if (alphaFudge) {
      (*prod_val)[i] = 0.5 * ((*prod_val)[i] + prod_val_temp[i]) * df;
    } else {
      (*prod_val)[i] *= df;
    }
  }

  (*prod_val)[*num_prod - 1] += xStr.gd->pv_of_past;

FREE_RETURN:

  if (time)
    free(time);
  if (date)
    free(date);

  if (vol_change)
    free(vol_change);
  if (sig_dom)
    free(sig_dom);
  if (sig_for)
    free(sig_for);
  if (mu_quanto_cons)
    free(mu_quanto_cons);
  if (mu_quanto_lin)
    free(mu_quanto_lin);
  if (sig_fx)
    free(sig_fx);
  if (corr_dom_for)
    free(corr_dom_for);
  if (corr_dom_fx)
    free(corr_dom_fx);
  if (corr_for_fx)
    free(corr_for_fx);

  if (dom_fwd)
    free(dom_fwd);
  if (dom_var)
    free(dom_var);
  if (dom_phi)
    free(dom_phi);
  if (dom_beta)
    free(dom_beta);

  if (for_fwd)
    free(for_fwd);
  if (for_var)
    free(for_var);
  if (for_phi)
    free(for_phi);
  if (for_beta)
    free(for_beta);

  if (fx_fwd)
    free(fx_fwd);
  if (fx_var)
    free(fx_var);
  if (ffx_var)
    free(ffx_var);

  if (global_corr_dom_for)
    free(global_corr_dom_for);
  if (global_corr_dom_fx)
    free(global_corr_dom_fx);
  if (global_corr_for_fx)
    free(global_corr_for_fx);

  if (bar_times)
    free(bar_times);

  if (is_event)
    free(is_event);

  if (bar_lvl)
    free(bar_lvl);
  if (is_bar)
    free(is_bar);
  if (bar_cl)
    free(bar_cl);

  if (void_prm) {
    for (i = 0; i < nstp; i++) {
      if (void_prm[i]) {
        grfn_prm = (GRFNPARMBETAMCDLM)void_prm[i];

        if (grfn_prm->do_fx && grfn_prm->num_fx_df > 0 &&
            grfn_prm->fx_idx != -1) {
          if (grfn_prm->df_fx_coef_const)
            free_dvector(grfn_prm->df_fx_coef_const, 0,
                         grfn_prm->num_fx_df - 1);
          if (grfn_prm->df_fx_coef_lin)
            free_dvector(grfn_prm->df_fx_coef_lin, 0, grfn_prm->num_fx_df - 1);
        }

        if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0 &&
            grfn_prm->dom_idx != -1) {
          if (grfn_prm->df_dom_coef_const)
            free_dvector(grfn_prm->df_dom_coef_const, 0,
                         grfn_prm->num_dom_df - 1);
          if (grfn_prm->df_dom_coef_lin)
            free_dvector(grfn_prm->df_dom_coef_lin, 0,
                         grfn_prm->num_dom_df - 1);
        }

        if (grfn_prm->do_for && grfn_prm->num_for_df > 0 &&
            grfn_prm->for_idx != -1) {
          if (grfn_prm->df_for_coef_const)
            free_dvector(grfn_prm->df_for_coef_const, 0,
                         grfn_prm->num_for_df - 1);
          if (grfn_prm->df_for_coef_lin)
            free_dvector(grfn_prm->df_for_coef_lin, 0,
                         grfn_prm->num_for_df - 1);
        }
      }
    }

    free(void_prm);
  }

  if (void_prm_down) {
    for (i = 0; i < nstp; i++) {
      if (void_prm_down[i]) {
        grfn_prm = (GRFNPARMBETAMCDLM)void_prm_down[i];

        if (grfn_prm->do_fx && grfn_prm->num_fx_df > 0 &&
            grfn_prm->fx_idx != -1) {
          if (grfn_prm->df_fx_coef_const)
            free_dvector(grfn_prm->df_fx_coef_const, 0,
                         grfn_prm->num_fx_df - 1);
          if (grfn_prm->df_fx_coef_lin)
            free_dvector(grfn_prm->df_fx_coef_lin, 0, grfn_prm->num_fx_df - 1);
        }

        if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0 &&
            grfn_prm->dom_idx != -1) {
          if (grfn_prm->df_dom_coef_const)
            free_dvector(grfn_prm->df_dom_coef_const, 0,
                         grfn_prm->num_dom_df - 1);
          if (grfn_prm->df_dom_coef_lin)
            free_dvector(grfn_prm->df_dom_coef_lin, 0,
                         grfn_prm->num_dom_df - 1);
        }

        if (grfn_prm->do_for && grfn_prm->num_for_df > 0 &&
            grfn_prm->for_idx != -1) {
          if (grfn_prm->df_for_coef_const)
            free_dvector(grfn_prm->df_for_coef_const, 0,
                         grfn_prm->num_for_df - 1);
          if (grfn_prm->df_for_coef_lin)
            free_dvector(grfn_prm->df_for_coef_lin, 0,
                         grfn_prm->num_for_df - 1);
        }
      }
    }

    free(void_prm_down);
  }

  if (free_str) {
    FIRSTFreeMktStruct(&xStr);
  }

  if (model) {
    free_FxBetaDLM_model(model);
    free(model);
  }

  if (prod_val_temp)
    free(prod_val_temp);

  return err;
}

char *SrtGrfn3DFxBetaDLMMc(char *und3dfx, int numeventdates, long *eventdates,
                           long tableauRows, long *tableauCols,
                           char ***tableauStrings, int **tableauMask,
                           long auxWidth, long *auxLen, double **aux,

                           // Params
                           FxBetaDLM_GRFNNumerParams *Params,

                           // for Optimisation of exercise boundary
                           int do_optimisation, int *optimise,
                           MCEBPARAMS params, long *resRows,

                           long num_paths, int *nb_prod, double ***prod_val) {
  int free_str = 0;
  FIRSTAllMkts xStr;
  SrtGrfnParam defParm;
  int forback;
  int flag = 0;
  long nstp;

  double *time = NULL, *date = NULL, *dom_fwd = NULL, *dom_std = NULL,
         *dom_phi = NULL, *dom_beta = NULL, *bond_pay_const = NULL,
         *bond_pay_lin = NULL, *for_ifr = NULL, *for_fwd_const = NULL,
         *for_fwd_lin = NULL, *for_std = NULL, *for_phi = NULL,
         *for_beta = NULL, *fx_fwd = NULL, *fx_std = NULL, *ffx_var = NULL,
         *dom_for_cov = NULL, *dom_fx_cov = NULL, *for_fx_cov = NULL,
         **prod_val_temp = NULL;

  void **void_prm = NULL;
  void **void_prm_down = NULL;

  GRFNPARMBETAMCDLM grfn_prm;

  long today;

  int i, j, k, num_col;
  SrtUndPtr fx_und;
  char *domname, *forname;

  int fx_idx, dom_idx, for_idx;

  double pay_time, df;
  long pay_date;
  double T;

  double var_ffx, temp, beta_temp;

  FxBetaDLM_model *model = NULL;

  clock_t t1, t2;

  int alphaFudge = 0, UpOrDown; /* For The Alpha Fudge */

  Err err = NULL;

  t1 = clock();

  /*	Initialise the GRFN tableau */

  /*	First        , initialise the param struct */

  err = srt_f_set_default_GrfnParams(&defParm);
  defParm.num_MCarlo_paths = num_paths;
  defParm.max_time_per_slice = 1000;
  defParm.min_nodes_per_path = 1;
  defParm.force_mc = 1;
  defParm.jumping = 1;

  err = FIRSTInitMktStruct(numeventdates, eventdates, tableauRows, *tableauCols,
                           tableauStrings, tableauMask, auxWidth, auxLen, aux,
                           und3dfx, &defParm, &forback, &xStr);

  if (err) {
    goto FREE_RETURN;
  }

  free_str = 1;

  /*	Now        , lookup underlyings involved and their term structures */

  fx_und = lookup_und(und3dfx);

  if (!fx_und) {
    err = serror("Couldn't find underlying named %s", und3dfx);
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(fx_und);

  model = calloc(1, sizeof(FxBetaDLM_model));

  if (!model) {
    err = "Memory allocation error in SrtGrfnFx3DBetaDLMTree";
    goto FREE_RETURN;
  }

  err = FxBetaDLM_Get_Model(und3dfx, model);

  if (err)
    goto FREE_RETURN;

  if (fabs(model->dAlpha) > 1e-6) {
    alphaFudge = 1;
  }

  domname = get_domname_from_fxund(fx_und);
  forname = get_forname_from_fxund(fx_und);

  num_col = xStr.num_cols;
  /*	Next        , get the time steps */

  /*	Copy event dates */
  nstp = xStr.num_evt;
  while (nstp >= 1 && xStr.evt[nstp - 1].evt == NULL) {
    nstp--;
  }
  if (nstp < 1) {
    err = "No event in Tableau";
    goto FREE_RETURN;
  }
  time = (double *)calloc(nstp, sizeof(double));
  if (!time) {
    err = "Memory allocation error (1) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }
  memcpy(time, xStr.tms, nstp * sizeof(double));

  date = (double *)calloc(nstp, sizeof(double));
  if (!date) {
    err = "Memory allocation error (3) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }

  for (i = 0; i < nstp; i++) {
    date[i] = today + DAYS_IN_YEAR * time[i];

    if (i > 0 && date[i] - date[i - 1] >= 1) {
      date[i] = (long)(date[i] + 1.0e-08);
      time[i] = YEARS_IN_DAY * (date[i] - today);
    }
  }

  if (time[0] > 0) {
    /* add the zero time */
    num_f_add_number(&nstp, &time, 0);
    num_f_sort_vector(nstp, time);
    nstp -= 1;
    num_f_add_number(&nstp, &date, today);
    num_f_sort_vector(nstp, date);
    flag = 1;
  }

  pay_date = model->lTStarDate;
  pay_time = model->dTStar;

  /*	Get distributions */
  dom_fwd = (double *)calloc(nstp, sizeof(double));
  dom_std = (double *)calloc(nstp, sizeof(double));
  dom_phi = (double *)calloc(nstp, sizeof(double));
  dom_beta = (double *)calloc(nstp, sizeof(double));

  bond_pay_const = (double *)calloc(nstp, sizeof(double));
  bond_pay_lin = (double *)calloc(nstp, sizeof(double));

  for_fwd_const = (double *)calloc(nstp, sizeof(double));
  for_fwd_lin = (double *)calloc(nstp, sizeof(double));
  for_std = (double *)calloc(nstp, sizeof(double));
  for_phi = (double *)calloc(nstp, sizeof(double));
  for_beta = (double *)calloc(nstp, sizeof(double));

  fx_fwd = (double *)calloc(nstp, sizeof(double));
  fx_std = (double *)calloc(nstp, sizeof(double));
  ffx_var = (double *)calloc(nstp, sizeof(double));

  dom_for_cov = (double *)calloc(nstp, sizeof(double));
  dom_fx_cov = (double *)calloc(nstp, sizeof(double));
  for_fx_cov = (double *)calloc(nstp, sizeof(double));

  if (!dom_fwd || !dom_std || !dom_phi || !dom_beta || !for_fwd_const ||
      !for_fwd_lin || !for_std || !for_phi || !for_beta || !fx_fwd || !fx_std ||
      !ffx_var || !bond_pay_const || !bond_pay_lin || !dom_for_cov ||
      !dom_fx_cov || !for_fx_cov) {
    err = "Memory allocation error (5) in SrtGrfn3DFXMc";
    goto FREE_RETURN;
  }

  for (UpOrDown = 0; UpOrDown <= alphaFudge; UpOrDown++) {

    if (alphaFudge) {
      model->dAlpha = 0.0;

      if (!UpOrDown) {
        memcpy(model->dSigmaFx, model->dSigmaFxUp,
               model->iNbPWTime * sizeof(double));
      } else {
        memcpy(model->dSigmaFx, model->dSigmaFxDown,
               model->iNbPWTime * sizeof(double));
        memcpy(model->dCorrDomFx, model->dCorrDomFxDown,
               model->iNbPWTime * sizeof(double));
        memcpy(model->dCorrForFx, model->dCorrForFxDown,
               model->iNbPWTime * sizeof(double));
      }
    }

    Fx3DBetaDLM_ExpectAndVarGrfn(
        nstp, time, model, dom_fwd, dom_std, dom_phi, dom_beta, for_fwd_const,
        for_fwd_lin, for_std, for_phi, for_beta, fx_fwd, fx_std, ffx_var, 0,
        dom_for_cov, dom_fx_cov, for_fx_cov, Params->min_time);

    /*	Fill product structure */

    strupper(und3dfx);
    strip_white_space(und3dfx);
    strupper(domname);
    strip_white_space(domname);
    strupper(forname);
    strip_white_space(forname);
    for (i = 0; i < xStr.num_und; i++) {
      strupper(xStr.und_data[i].und_name);
      strip_white_space(xStr.und_data[i].und_name);
    }

    fx_idx = -1;
    for (i = 0; i < xStr.num_und; i++) {
      if (!strcmp(xStr.und_data[i].und_name, und3dfx)) {
        fx_idx = i;
      }
    }
    if (fx_idx == -1) {
      err = "The Fx underlying is not present in the mdlcomm structure";
      goto FREE_RETURN;
    }

    dom_idx = -1;
    for (i = 0; i < xStr.num_und; i++) {
      if (!strcmp(xStr.und_data[i].und_name, domname)) {
        dom_idx = i;
      }
    }
    if (dom_idx == -1) {
      err = "The domestic underlying is not present in the mdlcomm structure";
      goto FREE_RETURN;
    }

    for_idx = -1;
    for (i = 0; i < xStr.num_und; i++) {
      if (!strcmp(xStr.und_data[i].und_name, forname)) {
        for_idx = i;
      }
    }
    if (for_idx == -1) {
      err = "The foreign underlying is not present in the mdlcomm structure";
      goto FREE_RETURN;
    }

    if (!UpOrDown) {
      void_prm = (void **)calloc(nstp, sizeof(void *));

      if (!void_prm) {
        err = "Memory allocation error (6) in SrtGrfn3DFXMc";
        goto FREE_RETURN;
      }
    } else {
      void_prm_down = (void **)calloc(nstp, sizeof(void *));

      if (!void_prm_down) {
        err = "Memory allocation error (6) in SrtGrfn3DFXMc";
        goto FREE_RETURN;
      }
    }

    df = swp_f_df(today, pay_date, model->cYcDom);

    for (i = xStr.num_evt - 1; i >= 0; i--) {
      if (xStr.evt[i].evt) {
        grfn_prm = malloc(sizeof(grfn_parm_beta_MC_DLM));

        if (!grfn_prm) {
          err = "Memory allocation error (7) in SrtGrfn3DFXMc";
          goto FREE_RETURN;
        }

        grfn_prm->global = &xStr;
        grfn_prm->local = xStr.evt + i;
        grfn_prm->fx_idx = fx_idx;
        grfn_prm->dom_idx = dom_idx;
        grfn_prm->for_idx = for_idx;

        /* Pay Bond reconstruction */
        bond_pay_const[i + flag] =
            log(swp_f_df(xStr.dts[i], pay_date, model->cYcDom));
        bond_pay_const[i + flag] +=
            0.5 * dom_beta[i + flag] * dom_beta[i + flag] * dom_phi[i + flag];
        bond_pay_lin[i + flag] = dom_beta[i + flag];

        /* Spot Fx reconstruction */
        var_ffx = ffx_var[i + flag];
        temp = 1.0 + 2.0 * model->dC0 * var_ffx;

        grfn_prm->fx_coef_const =
            log(model->dCashFx * swp_f_df(today, xStr.dts[i], model->cYcFor) /
                swp_f_df(today, xStr.dts[i], model->cYcDom));

        grfn_prm->fx_coef_const +=
            -0.5 * (model->dB0 * model->dB0 * var_ffx / temp + log(temp));

        grfn_prm->fx_coef_const +=
            0.5 * (dom_beta[i + flag] * dom_beta[i + flag] * dom_phi[i + flag] -
                   for_beta[i + flag] * for_beta[i + flag] * for_phi[i + flag]);

        grfn_prm->fx_coef_lin = model->dB0 / temp;
        grfn_prm->fx_coef_quad = model->dC0 / temp;
        grfn_prm->fx_coef_dom = dom_beta[i + flag];
        grfn_prm->fx_coef_for = -for_beta[i + flag];

        /* Fx discount factors */
        grfn_prm->num_fx_df = xStr.evt[i].evt->dflen[fx_idx];
        grfn_prm->fx_df_tms = xStr.evt[i].evt->dft[fx_idx];
        grfn_prm->fx_df_dts = xStr.evt[i].evt->dfd[fx_idx];

        if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1) {
          grfn_prm->df_fx_coef_const = dvector(0, grfn_prm->num_fx_df - 1);
          grfn_prm->df_fx_coef_lin = dvector(0, grfn_prm->num_fx_df - 1);

          if (!grfn_prm->df_fx_coef_const || !grfn_prm->df_fx_coef_lin) {
            err = "Memory allocation error (7) in SrtGrfn3DFXMc";
            goto FREE_RETURN;
          }

          for (j = 0; j < grfn_prm->num_fx_df; j++) {
            T = xStr.tms[i] + grfn_prm->fx_df_tms[j];
            beta_temp = (1.0 - exp(-model->dLambdaDom * (T - model->dTStar))) /
                        model->dLambdaDom;
            grfn_prm->df_fx_coef_const[j] = log(
                swp_f_df(xStr.dts[i], grfn_prm->fx_df_dts[j], model->cYcDom));
            grfn_prm->df_fx_coef_const[j] -=
                0.5 *
                (beta_temp * beta_temp -
                 dom_beta[i + flag] * dom_beta[i + flag]) *
                dom_phi[i + flag];
            grfn_prm->df_fx_coef_lin[j] = -(beta_temp - dom_beta[i + flag]);
          }

          grfn_prm->do_fx = 1;
        } else {
          grfn_prm->do_fx = 0;
        }

        /* Domestic discount factors */
        grfn_prm->num_dom_df = xStr.evt[i].evt->dflen[dom_idx];
        grfn_prm->dom_df_tms = xStr.evt[i].evt->dft[dom_idx];
        grfn_prm->dom_df_dts = xStr.evt[i].evt->dfd[dom_idx];

        if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1) {
          grfn_prm->df_dom_coef_const = dvector(0, grfn_prm->num_dom_df - 1);
          grfn_prm->df_dom_coef_lin = dvector(0, grfn_prm->num_dom_df - 1);

          if (!grfn_prm->df_dom_coef_const || !grfn_prm->df_dom_coef_lin) {
            err = "Memory allocation error (7) in SrtGrfn3DFXMc";
            goto FREE_RETURN;
          }

          for (j = 0; j < grfn_prm->num_dom_df; j++) {
            T = xStr.tms[i] + grfn_prm->dom_df_tms[j];
            beta_temp = (1.0 - exp(-model->dLambdaDom * (T - model->dTStar))) /
                        model->dLambdaDom;
            grfn_prm->df_dom_coef_const[j] = log(
                swp_f_df(xStr.dts[i], grfn_prm->dom_df_dts[j], model->cYcDom));
            grfn_prm->df_dom_coef_const[j] -=
                0.5 *
                (beta_temp * beta_temp -
                 dom_beta[i + flag] * dom_beta[i + flag]) *
                dom_phi[i + flag];
            grfn_prm->df_dom_coef_lin[j] = -(beta_temp - dom_beta[i + flag]);
          }

          grfn_prm->do_dom = 1;
        } else {
          grfn_prm->do_dom = 0;
        }

        /* Foreign discount factors */
        grfn_prm->num_for_df = xStr.evt[i].evt->dflen[for_idx];
        grfn_prm->for_df_tms = xStr.evt[i].evt->dft[for_idx];
        grfn_prm->for_df_dts = xStr.evt[i].evt->dfd[for_idx];

        if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1) {
          grfn_prm->df_for_coef_const = dvector(0, grfn_prm->num_for_df - 1);
          grfn_prm->df_for_coef_lin = dvector(0, grfn_prm->num_for_df - 1);

          if (!grfn_prm->df_for_coef_const || !grfn_prm->df_for_coef_lin) {
            err = "Memory allocation error (7) in SrtGrfn3DFXMc";
            goto FREE_RETURN;
          }

          for (j = 0; j < grfn_prm->num_for_df; j++) {
            T = xStr.tms[i] + grfn_prm->for_df_tms[j];
            beta_temp = (1.0 - exp(-model->dLambdaFor * (T - model->dTStar))) /
                        model->dLambdaFor;
            grfn_prm->df_for_coef_const[j] = log(
                swp_f_df(xStr.dts[i], grfn_prm->for_df_dts[j], model->cYcFor));
            grfn_prm->df_for_coef_const[j] -=
                0.5 *
                (beta_temp * beta_temp -
                 for_beta[i + flag] * for_beta[i + flag]) *
                for_phi[i + flag];
            grfn_prm->df_for_coef_lin[j] = -(beta_temp - for_beta[i + flag]);
          }

          grfn_prm->do_for = 1;
        } else {
          grfn_prm->do_for = 0;
        }

        if (!UpOrDown) {
          void_prm[i + flag] = (void *)grfn_prm;
        } else {
          void_prm_down[i + flag] = (void *)grfn_prm;
        }
      } else {
        if (!UpOrDown) {
          void_prm[i + flag] = NULL;
        } else {
          void_prm_down[i + flag] = NULL;
        }
      }
    }

    /*	Eventually! call to function */

    *nb_prod = num_col;
    *tableauCols = num_col;
    *resRows = max(num_col + 1, xStr.num_evt);

    if (!UpOrDown) {
      if (do_optimisation) {
        if (params->iMultiIndex) {
          params->iNbIndex = params->iMultiIndex;
        } else {
          params->iNbIndex = 1;
        }

        err = mceb_allocate_params(params, nstp);

        if (err)
          goto FREE_RETURN;

        *prod_val = dmatrix(0, *resRows - 1, 0, 2 + params->iNbIndex);
        prod_val_temp = dmatrix(0, *resRows - 1, 0, 2 + params->iNbIndex);
      } else {
        *prod_val = dmatrix(0, *nb_prod - 1, 0, 2);
        prod_val_temp = dmatrix(0, *nb_prod - 1, 0, 2);
      }

      if (!*prod_val || !prod_val_temp) {
        err = "Memory allocation failure in SrtGrfnMainLGMSVMC";
        goto FREE_RETURN;
      }
    }

    t2 = clock();

    smessage("Phase 1 -preprocessing        , time in sec: %.2f",
             (double)(t2 - t1) / CLOCKS_PER_SEC);

    err = mc_main_3dBetaDLMfx(
        /*	Time data */
        num_paths, num_col, time, date, nstp, dom_fwd, dom_std, for_fwd_const,
        for_fwd_lin, for_std, fx_fwd, fx_std, bond_pay_const, bond_pay_lin,
        dom_for_cov, dom_fx_cov, for_fx_cov,
        UpOrDown ? void_prm_down : void_prm, Params->do_pecs, do_optimisation,
        optimise, params, NULL,
        /*	Payoff function */
        grfn_payoff_4_3dfxBetaDLM_mc, /*	Result */
        *prod_val);

    if (do_optimisation) {
      /* Recopy Barrier / CoefLin for the moment */
      for (i = 0; i < nstp; i++) {
        (*prod_val)[i][2] = params->dBarrier[i];

        for (j = 0; j < params->iNbIndex; j++) {
          (*prod_val)[i][3 + j] = params->dCoefLin[i][j + 1];
        }
      }
    }

    if (!UpOrDown && alphaFudge) {
      for (i = 0; i < num_col; i++) {
        (prod_val_temp)[i][0] = (*prod_val)[i][0] * df;
        (prod_val_temp)[i][1] = (*prod_val)[i][1] * df;
      }
    } else {
      for (i = 0; i < num_col; i++) {
        (*prod_val)[i][0] *= df;
        (*prod_val)[i][1] *= df;
      }
    }

    if (do_optimisation) {
      if (!UpOrDown && alphaFudge) {
        (prod_val_temp)[num_col][0] = (*prod_val)[i][0] * df;
        (prod_val_temp)[num_col][1] = (*prod_val)[i][1] * df;
      } else {
        (*prod_val)[num_col][0] *= df;
        (*prod_val)[num_col][1] *= df;
      }
    }

    if (flag) {
      if (!UpOrDown && alphaFudge) {
        for (i = 0; i < xStr.num_evt - 1; i++) {
          if (do_optimisation) {
            (prod_val_temp)[i][2] = (prod_val_temp)[i + 1][2];

            for (k = 0; k < params->iNbIndex; k++) {
              (prod_val_temp)[i][3 + k] = (prod_val_temp)[i + 1][3 + k];
            }
          }
        }
      } else {
        for (i = 0; i < xStr.num_evt - 1; i++) {
          if (do_optimisation) {
            (*prod_val)[i][2] = (*prod_val)[i + 1][2];

            for (k = 0; k < params->iNbIndex; k++) {
              (*prod_val)[i][3 + k] = (*prod_val)[i + 1][3 + k];
            }
          }
        }
      }
    }
  }

  /* Average PV */
  if (alphaFudge) {
    for (i = 0; i < num_col; i++) {
      (*prod_val)[i][0] = 0.5 * ((*prod_val)[i][0] + (prod_val_temp)[i][0]);
      (*prod_val)[i][1] = 0.5 * ((*prod_val)[i][1] + (prod_val_temp)[i][1]);
    }

    if (do_optimisation) {
      (*prod_val)[num_col][0] =
          0.5 * ((*prod_val)[num_col][0] + (prod_val_temp)[num_col][0]);
      (*prod_val)[num_col][1] =
          0.5 * ((*prod_val)[num_col][1] + (prod_val_temp)[num_col][1]);
    }

    if (flag) {
      for (i = 0; i < xStr.num_evt - 1; i++) {
        if (do_optimisation) {
          (*prod_val)[i][2] = 0.5 * ((*prod_val)[i][2] + (prod_val_temp)[i][2]);

          for (k = 0; k < params->iNbIndex; k++) {
            (*prod_val)[i][3 + k] =
                0.5 * ((*prod_val)[i][3 + k] + (prod_val_temp)[i][3 + k]);
          }
        }
      }
    }
  }

  /*	Add PV of Past */
  (*prod_val)[num_col - 1][0] += xStr.gd->pv_of_past;

FREE_RETURN:

  if (time)
    free(time);
  if (date)
    free(date);

  if (dom_fwd)
    free(dom_fwd);
  if (dom_std)
    free(dom_std);
  if (dom_phi)
    free(dom_phi);
  if (dom_beta)
    free(dom_beta);

  if (bond_pay_const)
    free(bond_pay_const);
  if (bond_pay_lin)
    free(bond_pay_lin);

  if (for_fwd_const)
    free(for_fwd_const);
  if (for_fwd_lin)
    free(for_fwd_lin);
  if (for_std)
    free(for_std);
  if (for_phi)
    free(for_phi);
  if (for_beta)
    free(for_beta);

  if (fx_fwd)
    free(fx_fwd);
  if (fx_std)
    free(fx_std);
  if (ffx_var)
    free(ffx_var);

  if (dom_for_cov)
    free(dom_for_cov);
  if (dom_fx_cov)
    free(dom_fx_cov);
  if (for_fx_cov)
    free(for_fx_cov);

  if (void_prm) {
    for (i = 0; i < nstp; i++) {
      if (void_prm[i]) {
        grfn_prm = (GRFNPARMBETAMCDLM)void_prm[i];

        if (grfn_prm->do_fx && grfn_prm->num_fx_df > 0 &&
            grfn_prm->fx_idx != -1) {
          if (grfn_prm->df_fx_coef_const)
            free_dvector(grfn_prm->df_fx_coef_const, 0,
                         grfn_prm->num_fx_df - 1);
          if (grfn_prm->df_fx_coef_lin)
            free_dvector(grfn_prm->df_fx_coef_lin, 0, grfn_prm->num_fx_df - 1);
        }

        if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0 &&
            grfn_prm->dom_idx != -1) {
          if (grfn_prm->df_dom_coef_const)
            free_dvector(grfn_prm->df_dom_coef_const, 0,
                         grfn_prm->num_dom_df - 1);
          if (grfn_prm->df_dom_coef_lin)
            free_dvector(grfn_prm->df_dom_coef_lin, 0,
                         grfn_prm->num_dom_df - 1);
        }

        if (grfn_prm->do_for && grfn_prm->num_for_df > 0 &&
            grfn_prm->for_idx != -1) {
          if (grfn_prm->df_for_coef_const)
            free_dvector(grfn_prm->df_for_coef_const, 0,
                         grfn_prm->num_for_df - 1);
          if (grfn_prm->df_for_coef_lin)
            free_dvector(grfn_prm->df_for_coef_lin, 0,
                         grfn_prm->num_for_df - 1);
        }

        free(grfn_prm);
      }
    }

    free(void_prm);
  }

  if (void_prm_down) {
    for (i = 0; i < nstp; i++) {
      if (void_prm_down[i]) {
        grfn_prm = (GRFNPARMBETAMCDLM)void_prm_down[i];

        if (grfn_prm->do_fx && grfn_prm->num_fx_df > 0 &&
            grfn_prm->fx_idx != -1) {
          if (grfn_prm->df_fx_coef_const)
            free_dvector(grfn_prm->df_fx_coef_const, 0,
                         grfn_prm->num_fx_df - 1);
          if (grfn_prm->df_fx_coef_lin)
            free_dvector(grfn_prm->df_fx_coef_lin, 0, grfn_prm->num_fx_df - 1);
        }

        if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0 &&
            grfn_prm->dom_idx != -1) {
          if (grfn_prm->df_dom_coef_const)
            free_dvector(grfn_prm->df_dom_coef_const, 0,
                         grfn_prm->num_dom_df - 1);
          if (grfn_prm->df_dom_coef_lin)
            free_dvector(grfn_prm->df_dom_coef_lin, 0,
                         grfn_prm->num_dom_df - 1);
        }

        if (grfn_prm->do_for && grfn_prm->num_for_df > 0 &&
            grfn_prm->for_idx != -1) {
          if (grfn_prm->df_for_coef_const)
            free_dvector(grfn_prm->df_for_coef_const, 0,
                         grfn_prm->num_for_df - 1);
          if (grfn_prm->df_for_coef_lin)
            free_dvector(grfn_prm->df_for_coef_lin, 0,
                         grfn_prm->num_for_df - 1);
        }

        free(grfn_prm);
      }
    }

    free(void_prm_down);
  }

  if (free_str) {
    FIRSTFreeMktStruct(&xStr);
  }

  if (model) {
    free_FxBetaDLM_model(model);
    free(model);
  }

  if (prod_val_temp) {
    if (do_optimisation) {
      free_dmatrix(prod_val_temp, 0, *resRows - 1, 0, 2 + params->iNbIndex);
    } else {
      free_dmatrix(prod_val_temp, 0, *nb_prod - 1, 0, 2);
    }
  }
  return err;
}