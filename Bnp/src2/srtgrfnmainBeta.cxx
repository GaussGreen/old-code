/**********************************************************************
 *      Name: SrtGrfnMain.cxx                                           *
 *  Function: Entry point to GRFN with raw data                       *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 18/10/95                                                *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere (Excel or 2020)                  *
 *   Returns:                                                         *
 *   Globals: Expects mkt and request list structures to exist        *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 18/10/95 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/
#include "BGMEval.h"
#include "SrtAccess.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

#define MAX_TIME 0.02
#define MAX_STP 3000

char *SrtGrfn3DFXBetaTree(char *und3dfx, double alpha, double beta,
                          int numeventdates, long *eventdates, long tableauRows,
                          long tableauCols, char ***tableauStrings,
                          int **tableauMask, long auxWidth, long *auxLen,
                          double **aux, int is_end_of_day_fixing,
                          int is_end_of_day_payment, double *barrier,
                          int *bar_col, long num_stp, int *num_prod,
                          int discount, double **prod_val) {
  int free_str = 0;
  FIRSTAllMkts xStr;
  SrtGrfnParam defParm;
  int forback;
  long nstp;

  double *time = NULL, *date = NULL;
  int *vol_change = NULL;
  double *dom_vol = NULL, *for_vol = NULL, *fx_vol = NULL,
         *sig_fx_approx = NULL, *fx_fwd_approx = NULL;

  double *dom_ifr = NULL, *dom_fwd = NULL, *dom_var = NULL, *for_ifr = NULL,
         *for_fwd = NULL, *for_var = NULL, *fx_fwd = NULL, *fx_var = NULL,
         *corr_dom_for = NULL, *corr_dom_fx = NULL, *corr_for_fx = NULL;

  void **void_prm = NULL;
  GRFNPARMTREE grfn_prm;
  int *is_event = NULL;

  long today, spot_date;
  int i, j;
  SrtUndPtr fx_und, dom_und, for_und;
  TermStruct *fx_ts, *dom_ts, *for_ts;
  SrtCorrLstPtr sCorrlist;
  char *domname, *forname;
  double dom_lam, for_lam;
  double spot_fx;
  char *dom_yc, *for_yc;
  int fx_idx, dom_idx, for_idx;

  int num_vol_times;
  double *vol_times = NULL;

  int num_bar_times;
  double *bar_times = NULL;
  double *beta_tab = NULL;

  int *is_bar = NULL;
  double *bar_lvl = NULL;
  int *bar_cl = NULL;
  double next_bar;

  double beta2, logSpot;

  clock_t t1, t2;

  Err err = NULL;

  t1 = clock();

  /*	Check if it is lognormal */
  if (fabs(beta - 1.0) < 1.0E-04) {
    return SrtGrfn3DFXTree(und3dfx, numeventdates, eventdates, tableauRows,
                           tableauCols, tableauStrings, tableauMask, auxWidth,
                           auxLen, aux, is_end_of_day_fixing,
                           is_end_of_day_payment, barrier, bar_col, num_stp,
                           num_prod, discount, prod_val);
  }

  /*	Initialise the GRFN tableau */

  /*	First        , initialise the param struct */

  beta2 = 1.0 - beta;

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

  fx_und = lookup_und(und3dfx);

  if (!fx_und) {
    err = serror("Couldn't find underlying named %s", und3dfx);
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(fx_und);

  if (get_underlying_type(fx_und) != FOREX_UND) {
    err = serror("Underlying %s is not of type FX", und3dfx);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_fxund(fx_und) != FX_STOCH_RATES) {
    err = serror("Underlying %s is not of type FX Stoch Rates", und3dfx);
    goto FREE_RETURN;
  }

  fx_ts = get_ts_from_fxund(fx_und);

  domname = get_domname_from_fxund(fx_und);
  dom_und = lookup_und(domname);
  if (!dom_und) {
    err = serror("Couldn't find underlying named %s", domname);
    goto FREE_RETURN;
  }
  dom_ts = get_ts_from_irund(dom_und);

  forname = get_forname_from_fxund(fx_und);
  for_und = lookup_und(forname);
  if (!for_und) {
    err = serror("Couldn't find underlying named %s", forname);
    goto FREE_RETURN;
  }
  for_ts = get_ts_from_irund(for_und);

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

  /*	Get vol dates */
  err = compute_vol_times(und3dfx, &num_vol_times, &vol_times, time[nstp - 1]);
  if (err) {
    goto FREE_RETURN;
  }

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

  err = fill_time_vector(&time, &nstp, num_bar_times, bar_times, num_vol_times,
                         vol_times, num_stp);

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

  sCorrlist = srt_f_GetTheCorrelationList();
  if (!sCorrlist->head->element) {
    err = "correlation list improperly initialised";
    goto FREE_RETURN;
  }

  vol_change = (int *)calloc(nstp, sizeof(int));
  dom_vol = (double *)calloc(nstp, sizeof(double));
  for_vol = (double *)calloc(nstp, sizeof(double));
  fx_vol = (double *)calloc(nstp, sizeof(double));

  beta_tab = (double *)calloc(nstp, sizeof(double));

  corr_dom_for = (double *)calloc(nstp, sizeof(double));
  corr_dom_fx = (double *)calloc(nstp, sizeof(double));
  corr_for_fx = (double *)calloc(nstp, sizeof(double));

  if (!vol_change || !dom_vol || !for_vol || !fx_vol || !corr_dom_for ||
      !corr_dom_fx || !corr_for_fx || !beta_tab) {
    err = "Memory allocation error (4) in SrtGrfn3DFXBetaTree";
    goto FREE_RETURN;
  }

  vol_change[nstp - 1] = 1;
  dom_vol[nstp - 1] = find_sig(time[nstp - 1], dom_ts);
  for_vol[nstp - 1] = find_sig(time[nstp - 1], for_ts);
  fx_vol[nstp - 1] =
      find_fx_sig(time[nstp - 1], fx_ts) *
      exp(alpha * sqrt(time[nstp - 1]) - 0.5 * alpha * alpha * time[nstp - 1]);
  beta_tab[nstp - 1] = beta;

  err = srt_f_get_corr_from_CorrList(sCorrlist, domname, forname,
                                     time[nstp - 1], &(corr_dom_for[nstp - 1]));
  if (err) {
    goto FREE_RETURN;
  }
  err = srt_f_get_corr_from_CorrList(sCorrlist, domname, und3dfx,
                                     time[nstp - 1], &(corr_dom_fx[nstp - 1]));
  if (err) {
    goto FREE_RETURN;
  }
  err = srt_f_get_corr_from_CorrList(sCorrlist, und3dfx, forname,
                                     time[nstp - 1], &(corr_for_fx[nstp - 1]));
  if (err) {
    goto FREE_RETURN;
  }

  for (i = nstp - 2; i >= 0; i--) {
    dom_vol[i] = find_sig(time[i], dom_ts);
    for_vol[i] = find_sig(time[i], for_ts);
    fx_vol[i] = find_fx_sig(time[i], fx_ts) *
                exp(alpha * sqrt(time[i]) - 0.5 * alpha * alpha * time[i]);

    err = srt_f_get_corr_from_CorrList(sCorrlist, domname, forname, time[i],
                                       &(corr_dom_for[i]));
    if (err)
      goto FREE_RETURN;
    err = srt_f_get_corr_from_CorrList(sCorrlist, domname, und3dfx, time[i],
                                       &(corr_dom_fx[i]));
    if (err)
      goto FREE_RETURN;
    err = srt_f_get_corr_from_CorrList(sCorrlist, und3dfx, forname, time[i],
                                       &(corr_for_fx[i]));
    if (err)
      goto FREE_RETURN;

    beta_tab[i] = beta;

    if (fabs(dom_vol[i] - dom_vol[i + 1]) + fabs(for_vol[i] - for_vol[i + 1]) +
            fabs(fx_vol[i] - fx_vol[i + 1]) +
            fabs(corr_dom_for[i] - corr_dom_for[i + 1]) +
            fabs(corr_dom_fx[i] - corr_dom_fx[i + 1]) +
            fabs(corr_for_fx[i] - corr_for_fx[i + 1]) >
        EPS) {
      vol_change[i] = 1;
    } else {
      vol_change[i] = 0;
    }
  }

  /*	Get lambdas and correls */

  err = get_lambda_from_ir_ts(dom_ts, &dom_lam);
  err = get_lambda_from_ir_ts(for_ts, &for_lam);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Get Fx spot and yield curves */

  dom_yc = get_ycname_from_irund(dom_und);
  for_yc = get_ycname_from_irund(for_und);

  spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  spot_fx = get_spot_from_fxund(fx_und) * swp_f_df(today, spot_date, dom_yc) /
            swp_f_df(today, spot_date, for_yc);

  /*	Get distributions */

  dom_ifr = (double *)calloc(nstp, sizeof(double));
  dom_fwd = (double *)calloc(nstp, sizeof(double));
  dom_var = (double *)calloc(nstp, sizeof(double));
  for_ifr = (double *)calloc(nstp, sizeof(double));
  for_fwd = (double *)calloc(nstp, sizeof(double));
  for_var = (double *)calloc(nstp, sizeof(double));
  fx_fwd = (double *)calloc(nstp, sizeof(double));
  fx_var = (double *)calloc(nstp, sizeof(double));

  sig_fx_approx = (double *)calloc(nstp, sizeof(double));
  fx_fwd_approx = (double *)calloc(nstp, sizeof(double));

  if (!dom_ifr || !dom_fwd || !dom_var || !for_ifr || !for_fwd || !for_var ||
      !fx_fwd || !fx_var || !sig_fx_approx || !fx_fwd_approx) {
    err = "Memory allocation error (5) in SrtGrfn3DFXBetaTree";
    goto FREE_RETURN;
  }

  logSpot = log(spot_fx);

  /* first get the coresponding lognormal volatilities */
  err = Fxbeta_log_approx_corr(today, time, nstp, dom_vol, for_vol, time, nstp,
                               fx_vol, beta_tab, dom_lam, for_lam, time,
                               corr_dom_for, corr_dom_fx, corr_for_fx, nstp,
                               spot_fx, dom_yc, for_yc, time, nstp,
                               fx_fwd_approx, sig_fx_approx, MAX_TIME);

  fill_fwd_var_corr(nstp, time, date, dom_vol, for_vol, sig_fx_approx, dom_lam,
                    for_lam, corr_dom_for, corr_dom_fx, corr_for_fx, dom_yc,
                    for_yc, dom_ifr, dom_fwd, dom_var, for_ifr, for_fwd,
                    for_var, fx_fwd, fx_var);

  /*	Overwrite beta-forward fx by a more accurate calculation
  for (i=0; i<nstp; i++)
  {
          err = Fx3DFwdBeta (und3dfx        , time[i]        , time[i]        ,
  0        ,
  &(fx_fwd[i])); if (err)
          {
                  goto FREE_RETURN;
          }
  }
  */

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

  is_event = (int *)calloc(nstp, sizeof(int));
  void_prm = (void **)calloc(nstp, sizeof(void *));

  is_bar = (int *)calloc(nstp, sizeof(int));
  bar_lvl = (double *)calloc(nstp, sizeof(double));
  bar_cl = (int *)calloc(nstp, sizeof(int));

  if (!is_event || !void_prm || !is_bar || !bar_lvl || !bar_cl) {
    err = "Memory allocation error (6) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }

  j = xStr.num_evt - 1;
  next_bar = 0.0;
  for (i = nstp - 1; i >= 0; i--) {
    if (j >= 0 && fabs(date[i] - xStr.dts[j]) < 1.0e-04) {
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

      is_event[i] = 1;
      void_prm[i] = (void *)grfn_prm;

      if (barrier[j] > 1.0e-08) {
        next_bar = (exp(beta2 * log(barrier[j])) - 1.0) / beta2;
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
        next_bar = (exp(beta2 * log(barrier[j])) - 1.0) / beta2;
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
      bar_lvl[i] = (exp(beta2 * (fx_fwd_approx[i])) - 1.0) / beta2;
    }

    if ((i > 0) && (fabs(next_bar) > 1.0e-08)) {
      next_bar = (next_bar + (for_ifr[i - 1] + for_fwd[i - 1] - dom_ifr[i - 1] -
                              dom_fwd[i - 1]) *
                                 (time[i] - time[i - 1])) /
                 (1.0 - beta2 *
                            (for_ifr[i - 1] + for_fwd[i - 1] - dom_ifr[i - 1] -
                             dom_fwd[i - 1]) *
                            (time[i] - time[i - 1]));
    }
  }

  is_bar[0] = 0;
  bar_lvl[0] = (exp(beta2 * (fx_fwd_approx[0])) - 1.0) / beta2;
  bar_cl[0] = -1;

  is_bar[1] = 0;
  bar_lvl[1] = (exp(beta2 * (fx_fwd_approx[1])) - 1.0) / beta2;
  bar_cl[1] = -1;

  /*	Eventually! call to function */

  *num_prod = xStr.num_cols;
  *prod_val = (double *)calloc(*num_prod + 1, sizeof(double));

  if (!(*prod_val)) {
    err = "Memory allocation error (7) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }

  t2 = clock();

  smessage("Phase 1 -preprocessing        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);
  smessage("Number of times steps        , required: %d        , actual: %d",
           num_stp, nstp);

  err = treeBeta_main_3dfx(
      nstp, time, date, vol_change, dom_vol, for_vol, fx_vol, beta_tab, dom_ifr,
      dom_fwd, dom_var, for_ifr, for_fwd, for_var, fx_fwd, fx_var, void_prm,
      is_event, bar_lvl, bar_cl, is_bar, dom_lam, for_lam, corr_dom_for,
      corr_dom_fx, corr_for_fx, spot_fx, dom_yc, for_yc,
      grfn_payoff_4_3dfxBeta_tree, *num_prod, discount, *prod_val);

  /*	Add PV of Past */
  (*prod_val)[*num_prod - 1] += xStr.gd->pv_of_past;

FREE_RETURN:

  if (time)
    free(time);
  if (date)
    free(date);
  if (vol_change)
    free(vol_change);
  if (dom_vol)
    free(dom_vol);
  if (for_vol)
    free(for_vol);
  if (fx_vol)
    free(fx_vol);
  if (dom_ifr)
    free(dom_ifr);
  if (dom_fwd)
    free(dom_fwd);
  if (dom_var)
    free(dom_var);
  if (for_ifr)
    free(for_ifr);
  if (for_fwd)
    free(for_fwd);
  if (for_var)
    free(for_var);
  if (fx_fwd)
    free(fx_fwd);
  if (fx_var)
    free(fx_var);
  if (corr_dom_for)
    free(corr_dom_for);
  if (corr_dom_fx)
    free(corr_dom_fx);
  if (corr_for_fx)
    free(corr_for_fx);

  if (fx_fwd_approx)
    free(fx_fwd_approx);
  if (sig_fx_approx)
    free(sig_fx_approx);
  if (beta_tab)
    free(beta_tab);

  if (vol_times)
    free(vol_times);
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
        grfn_prm = (GRFNPARMTREE)void_prm[i];
        free(grfn_prm);
      }
    }

    free(void_prm);
  }

  if (free_str) {
    FIRSTFreeMktStruct(&xStr);
  }

  return err;
}

char *SrtGrfn3DFXAlphaBetaTree(char *und3dfx, double alpha, double beta,
                               int numeventdates, long *eventdates,
                               long tableauRows, long tableauCols,
                               char ***tableauStrings, int **tableauMask,
                               long auxWidth, long *auxLen, double **aux,
                               int is_end_of_day_fixing,
                               int is_end_of_day_payment, double *barrier,
                               int *bar_col, long num_stp, int *num_prod,
                               int discount, double **prod_val) {
  char *cerr = NULL;
  double *pv1 = NULL, *pv2 = NULL;
  int i, npv;

  *prod_val = NULL;

  if (fabs(alpha) < 1.0e-04) {
    alpha = 0.0;
  }

  cerr = SrtGrfn3DFXBetaTree(
      und3dfx, alpha, beta, numeventdates, eventdates, tableauRows, tableauCols,
      tableauStrings, tableauMask, auxWidth, auxLen, aux, is_end_of_day_fixing,
      is_end_of_day_payment, barrier, bar_col, num_stp, &npv, discount, &pv1);

  if (cerr) {
    goto FREE_RETURN;
  }

  if (fabs(alpha) >= 1.0e-04) {
    cerr = SrtGrfn3DFXBetaTree(und3dfx, -alpha, beta, numeventdates, eventdates,
                               tableauRows, tableauCols, tableauStrings,
                               tableauMask, auxWidth, auxLen, aux,
                               is_end_of_day_fixing, is_end_of_day_payment,
                               barrier, bar_col, num_stp, &npv, discount, &pv2);

    if (cerr) {
      goto FREE_RETURN;
    }
  } else {
    pv2 = (double *)calloc(npv, sizeof(double));
    if (!pv2) {
      cerr = "Memory allocation error in SrtGrfn3DFXAlphaBetaTree";
      goto FREE_RETURN;
    }
    memcpy(pv2, pv1, npv * sizeof(double));
  }

  *prod_val = (double *)calloc(npv, sizeof(double));

  if (!(*prod_val)) {
    cerr = "Memory allocation error in SrtGrfn3DFXAlphaBetaTree";
    goto FREE_RETURN;
  }

  for (i = 0; i < npv; i++) {
    (*prod_val)[i] = 0.5 * (pv1[i] + pv2[i]);
  }

  *num_prod = npv;

FREE_RETURN:

  if (pv1) {
    free(pv1);
  }

  if (pv2) {
    free(pv2);
  }

  if (cerr) {
    free(*prod_val);
  }

  return cerr;
}

char *SrtGrfn3DFXAlphaBetaMcTree(char *und3dfx, double alpha, double lambda,
                                 double beta, long nbPaths, int numeventdates,
                                 long *eventdates, long tableauRows,
                                 long tableauCols, char ***tableauStrings,
                                 int **tableauMask, long auxWidth, long *auxLen,
                                 double **aux, double *barrier, int *bar_col,
                                 long num_stp, int *num_prod, int discount,
                                 double ***prod_val) {
  int free_str = 0;
  FIRSTAllMkts xStr;
  SrtGrfnParam defParm;
  int forback;
  long nstp, path;

  double *time = NULL, *date = NULL;
  int *vol_change = NULL;
  double *dom_vol = NULL, *for_vol = NULL, *fx_vol = NULL,
         *sig_fx_approx = NULL, *fx_fwd_approx = NULL;

  double *dom_ifr = NULL, *dom_fwd = NULL, *dom_var = NULL, *for_ifr = NULL,
         *for_fwd = NULL, *for_var = NULL, *fx_fwd = NULL, *fx_var = NULL,
         *corr_dom_for = NULL, *corr_dom_fx = NULL, *corr_for_fx = NULL;

  void **void_prm = NULL;
  GRFNPARMTREE grfn_prm;
  int *is_event = NULL;

  long today, spot_date;
  int i, j;
  SrtUndPtr fx_und, dom_und, for_und;
  TermStruct *fx_ts, *dom_ts, *for_ts;
  SrtCorrLstPtr sCorrlist;
  char *domname, *forname;
  double dom_lam, for_lam;
  double spot_fx;
  char *dom_yc, *for_yc;
  int fx_idx, dom_idx, for_idx;

  int num_vol_times;
  double *vol_times = NULL;

  int num_bar_times;
  double *bar_times = NULL;
  double *beta_tab = NULL;

  int *is_bar = NULL;
  double *bar_lvl = NULL;
  int *bar_cl = NULL;
  double next_bar;

  double beta2, logSpot;

  long start_nb = -123456789;
  double *prod_val_aux;

  double U;

  clock_t t1, t2;

  Err err = NULL;

  double def, t, prev_t, drift, dt;

  t1 = clock();

  /*	Initialise the GRFN tableau */

  /*	First        , initialise the param struct */

  beta2 = 1.0 - beta;

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

  fx_und = lookup_und(und3dfx);

  if (!fx_und) {
    err = serror("Couldn't find underlying named %s", und3dfx);
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(fx_und);

  if (get_underlying_type(fx_und) != FOREX_UND) {
    err = serror("Underlying %s is not of type FX", und3dfx);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_fxund(fx_und) != FX_STOCH_RATES) {
    err = serror("Underlying %s is not of type FX Stoch Rates", und3dfx);
    goto FREE_RETURN;
  }

  fx_ts = get_ts_from_fxund(fx_und);

  domname = get_domname_from_fxund(fx_und);
  dom_und = lookup_und(domname);
  if (!dom_und) {
    err = serror("Couldn't find underlying named %s", domname);
    goto FREE_RETURN;
  }
  dom_ts = get_ts_from_irund(dom_und);

  forname = get_forname_from_fxund(fx_und);
  for_und = lookup_und(forname);
  if (!for_und) {
    err = serror("Couldn't find underlying named %s", forname);
    goto FREE_RETURN;
  }
  for_ts = get_ts_from_irund(for_und);

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

  /*	Get vol dates */
  err = compute_vol_times(und3dfx, &num_vol_times, &vol_times, time[nstp - 1]);
  if (err) {
    goto FREE_RETURN;
  }

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

  err = fill_time_vector(&time, &nstp, num_bar_times, bar_times, num_vol_times,
                         vol_times, num_stp);

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

  sCorrlist = srt_f_GetTheCorrelationList();
  if (!sCorrlist->head->element) {
    err = "correlation list improperly initialised";
    goto FREE_RETURN;
  }

  vol_change = (int *)calloc(nstp, sizeof(int));
  dom_vol = (double *)calloc(nstp, sizeof(double));
  for_vol = (double *)calloc(nstp, sizeof(double));
  fx_vol = (double *)calloc(nstp, sizeof(double));

  beta_tab = (double *)calloc(nstp, sizeof(double));

  corr_dom_for = (double *)calloc(nstp, sizeof(double));
  corr_dom_fx = (double *)calloc(nstp, sizeof(double));
  corr_for_fx = (double *)calloc(nstp, sizeof(double));

  if (!vol_change || !dom_vol || !for_vol || !fx_vol || !corr_dom_for ||
      !corr_dom_fx || !corr_for_fx || !beta_tab) {
    err = "Memory allocation error (4) in SrtGrfn3DFXAlphaBetaTree";
    goto FREE_RETURN;
  }

  /*	Get lambdas and correls */
  err = get_lambda_from_ir_ts(dom_ts, &dom_lam);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_lambda_from_ir_ts(for_ts, &for_lam);
  if (err) {
    goto FREE_RETURN;
  }

  /*	Get Fx spot and yield curves */

  dom_yc = get_ycname_from_irund(dom_und);
  for_yc = get_ycname_from_irund(for_und);

  spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  spot_fx = get_spot_from_fxund(fx_und) * swp_f_df(today, spot_date, dom_yc) /
            swp_f_df(today, spot_date, for_yc);

  /*	Get distributions */

  dom_ifr = (double *)calloc(nstp, sizeof(double));
  dom_fwd = (double *)calloc(nstp, sizeof(double));
  dom_var = (double *)calloc(nstp, sizeof(double));
  for_ifr = (double *)calloc(nstp, sizeof(double));
  for_fwd = (double *)calloc(nstp, sizeof(double));
  for_var = (double *)calloc(nstp, sizeof(double));
  fx_fwd = (double *)calloc(nstp, sizeof(double));
  fx_var = (double *)calloc(nstp, sizeof(double));

  sig_fx_approx = (double *)calloc(nstp, sizeof(double));
  fx_fwd_approx = (double *)calloc(nstp, sizeof(double));

  if (!dom_ifr || !dom_fwd || !dom_var || !for_ifr || !for_fwd || !for_var ||
      !fx_fwd || !fx_var || !sig_fx_approx || !fx_fwd_approx) {
    err = "Memory allocation error (5) in SrtGrfn3DFXBetaTree";
    goto FREE_RETURN;
  }

  logSpot = log(spot_fx);

  is_event = (int *)calloc(nstp, sizeof(int));
  void_prm = (void **)calloc(nstp, sizeof(void *));

  is_bar = (int *)calloc(nstp, sizeof(int));
  bar_lvl = (double *)calloc(nstp, sizeof(double));
  bar_cl = (int *)calloc(nstp, sizeof(int));

  if (!is_event || !void_prm || !is_bar || !bar_lvl || !bar_cl) {
    err = "Memory allocation error (6) in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }

  *num_prod = xStr.num_cols;
  *prod_val = dmatrix(0, *num_prod + 1 - 1, 0, 1);
  prod_val_aux = dvector(0, *num_prod + 1 - 1);

  if (!(*prod_val) || !prod_val_aux) {
    err = "Memory allocation error (7) in SrtGrfn3DAlphaBetaFXTree";
    goto FREE_RETURN;
  }

  drift = -(lambda + 0.5 * alpha * alpha);

  dom_vol[0] = find_sig(0, dom_ts);
  for_vol[0] = find_sig(0, for_ts);
  fx_vol[0] = find_fx_sig(0, fx_ts);
  beta_tab[0] = beta;
  err = srt_f_get_corr_from_CorrList(sCorrlist, domname, forname, 0.0,
                                     &(corr_dom_for[0]));
  if (err)
    goto FREE_RETURN;
  err = srt_f_get_corr_from_CorrList(sCorrlist, domname, und3dfx, 0.0,
                                     &(corr_dom_fx[0]));
  if (err)
    goto FREE_RETURN;
  err = srt_f_get_corr_from_CorrList(sCorrlist, und3dfx, forname, 0.0,
                                     &(corr_for_fx[0]));
  if (err)
    goto FREE_RETURN;

  for (i = 1; i < nstp; i++) {
    dom_vol[i] = find_sig(time[i], dom_ts);
    for_vol[i] = find_sig(time[i], for_ts);

    err = srt_f_get_corr_from_CorrList(sCorrlist, domname, forname, time[i],
                                       &(corr_dom_for[i]));
    if (err)
      goto FREE_RETURN;
    err = srt_f_get_corr_from_CorrList(sCorrlist, domname, und3dfx, time[i],
                                       &(corr_dom_fx[i]));
    if (err)
      goto FREE_RETURN;
    err = srt_f_get_corr_from_CorrList(sCorrlist, und3dfx, forname, time[i],
                                       &(corr_for_fx[i]));
    if (err)
      goto FREE_RETURN;

    beta_tab[i] = beta;
  }

  for (path = 0; path < nbPaths; path++) {
    smessage("Starting path %d\n", path);

    def = 0;
    t = 0;

    vol_change[nstp - 1] = 1;

    for (i = 1; i < nstp; i++) {
      prev_t = t;
      t = time[i];
      dt = (t - prev_t);

      U = inv_cumnorm_fast(uniform(&start_nb));
      def += drift * dt + alpha * sqrt(dt) * U;

      fx_vol[i] = find_fx_sig(time[i], fx_ts) * exp(def);

      if (fabs(dom_vol[i] - dom_vol[i + 1]) +
              fabs(for_vol[i] - for_vol[i + 1]) +
              fabs(fx_vol[i] - fx_vol[i + 1]) +
              fabs(corr_dom_for[i] - corr_dom_for[i + 1]) +
              fabs(corr_dom_fx[i] - corr_dom_fx[i + 1]) +
              fabs(corr_for_fx[i] - corr_for_fx[i + 1]) >
          EPS) {
        vol_change[i - 1] = 1;
      } else {
        vol_change[i - 1] = 0;
      }
    }

    /* first get the coresponding lognormal volatilities */
    err = Fxbeta_log_approx_corr(today, time, nstp, dom_vol, for_vol, time,
                                 nstp, fx_vol, beta_tab, dom_lam, for_lam, time,
                                 corr_dom_for, corr_dom_fx, corr_for_fx, nstp,
                                 spot_fx, dom_yc, for_yc, time, nstp,
                                 fx_fwd_approx, sig_fx_approx, MAX_TIME);

    fill_fwd_var_corr(nstp, time, date, dom_vol, for_vol, sig_fx_approx,
                      dom_lam, for_lam, corr_dom_for, corr_dom_fx, corr_for_fx,
                      dom_yc, for_yc, dom_ifr, dom_fwd, dom_var, for_ifr,
                      for_fwd, for_var, fx_fwd, fx_var);

    /*	Overwrite beta-forward fx by a more accurate calculation
    for (i=0; i<nstp; i++)
    {
            err = Fx3DFwdBeta (und3dfx        , time[i]        , time[i] , 0 ,
    &(fx_fwd[i])); if (err)
            {
                    goto FREE_RETURN;
            }
    }
    */

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

    j = xStr.num_evt - 1;
    next_bar = 0.0;
    for (i = nstp - 1; i >= 0; i--) {
      if (j >= 0 && fabs(date[i] - xStr.dts[j]) < 1.0e-04) {
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

        is_event[i] = 1;
        void_prm[i] = (void *)grfn_prm;

        if (barrier[j] > 1.0e-08) {
          next_bar = (exp(beta2 * log(barrier[j])) - 1.0) / beta2;
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
      } else if (j > 0 && xStr.am[j - 1]) {
        grfn_prm = malloc(sizeof(grfn_parm_tree));
        grfn_prm->global = &xStr;
        grfn_prm->local = xStr.evt + j - 1;
        grfn_prm->fx_idx = fx_idx;
        grfn_prm->dom_idx = dom_idx;
        grfn_prm->for_idx = for_idx;

        grfn_prm->num_fx_df = xStr.evt[j - 1].evt->dflen[fx_idx];
        grfn_prm->fx_df_tms = xStr.evt[j - 1].evt->dft[fx_idx];
        grfn_prm->fx_df_dts = xStr.evt[j - 1].evt->dfd[fx_idx];

        grfn_prm->num_dom_df = xStr.evt[j - 1].evt->dflen[dom_idx];
        grfn_prm->dom_df_tms = xStr.evt[j - 1].evt->dft[dom_idx];
        grfn_prm->dom_df_dts = xStr.evt[j - 1].evt->dfd[dom_idx];

        grfn_prm->num_for_df = xStr.evt[j - 1].evt->dflen[for_idx];
        grfn_prm->for_df_tms = xStr.evt[j - 1].evt->dft[for_idx];
        grfn_prm->for_df_dts = xStr.evt[j - 1].evt->dfd[for_idx];

        if (barrier[j - 1] > 1.0e-08) {
          next_bar = (exp(beta2 * log(barrier[j - 1])) - 1.0) / beta2;
          if (bar_col[j - 1] >= 0 && bar_col[j - 1] <= xStr.num_cols - 1) {
            is_bar[i] = 1;
            bar_cl[i] = bar_col[j - 1];
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
        bar_lvl[i] = (exp(beta2 * (fx_fwd_approx[i])) - 1.0) / beta2;
      }

      if (i > 0) {
        next_bar = (next_bar + (for_ifr[i - 1] + for_fwd[i - 1] -
                                dom_ifr[i - 1] - dom_fwd[i - 1]) *
                                   (time[i] - time[i - 1])) /
                   (1.0 - beta2 *
                              (for_ifr[i - 1] + for_fwd[i - 1] -
                               dom_ifr[i - 1] - dom_fwd[i - 1]) *
                              (time[i] - time[i - 1]));
      }
    }

    is_bar[0] = 0;
    bar_lvl[0] = (exp(beta2 * (fx_fwd_approx[0])) - 1.0) / beta2;
    bar_cl[0] = -1;

    is_bar[1] = 0;
    bar_lvl[1] = (exp(beta2 * (fx_fwd_approx[1])) - 1.0) / beta2;
    bar_cl[1] = -1;

    /*	Eventually! call to function */

    t2 = clock();

    smessage("Phase 1 -preprocessing        , time in sec: %.2f",
             (double)(t2 - t1) / CLOCKS_PER_SEC);
    smessage("Number of times steps        , required: %d        , actual: %d",
             num_stp, nstp);

    err = treeBeta_main_3dfx(
        nstp, time, date, vol_change, dom_vol, for_vol, fx_vol, beta_tab,
        dom_ifr, dom_fwd, dom_var, for_ifr, for_fwd, for_var, fx_fwd, fx_var,
        void_prm, is_event, bar_lvl, bar_cl, is_bar, dom_lam, for_lam,
        corr_dom_for, corr_dom_fx, corr_for_fx, spot_fx, dom_yc, for_yc,
        grfn_payoff_4_3dfxBeta_tree, *num_prod, discount, prod_val_aux);

    /*	Add PV of Past */
    prod_val_aux[*num_prod - 1] += xStr.gd->pv_of_past;

    for (i = 0; i < *num_prod; i++) {
      (*prod_val)[i][0] += prod_val_aux[i] / nbPaths;
      (*prod_val)[i][1] += prod_val_aux[i] * prod_val_aux[i] / nbPaths;
    }

    if (void_prm) {
      for (i = 0; i < nstp; i++) {
        if (void_prm[i]) {
          grfn_prm = (GRFNPARMTREE)void_prm[i];
          free(grfn_prm);
          void_prm[i] = NULL;
        }
      }
    }
  }

  for (i = 0; i < *num_prod; i++) {
    (*prod_val)[i][1] = sqrt(
        ((*prod_val)[i][1] - (*prod_val)[i][0] * (*prod_val)[i][0]) / nbPaths);
  }

FREE_RETURN:

  if (time)
    free(time);
  if (date)
    free(date);
  if (vol_change)
    free(vol_change);
  if (dom_vol)
    free(dom_vol);
  if (for_vol)
    free(for_vol);
  if (fx_vol)
    free(fx_vol);
  if (dom_ifr)
    free(dom_ifr);
  if (dom_fwd)
    free(dom_fwd);
  if (dom_var)
    free(dom_var);
  if (for_ifr)
    free(for_ifr);
  if (for_fwd)
    free(for_fwd);
  if (for_var)
    free(for_var);
  if (fx_fwd)
    free(fx_fwd);
  if (fx_var)
    free(fx_var);
  if (corr_dom_for)
    free(corr_dom_for);
  if (corr_dom_fx)
    free(corr_dom_fx);
  if (corr_for_fx)
    free(corr_for_fx);

  if (fx_fwd_approx)
    free(fx_fwd_approx);
  if (sig_fx_approx)
    free(sig_fx_approx);
  if (beta_tab)
    free(beta_tab);

  if (vol_times)
    free(vol_times);
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
    free(void_prm);
  }

  if (free_str) {
    FIRSTFreeMktStruct(&xStr);
  }

  if (prod_val_aux) {
    free_dvector(prod_val_aux, 0, *num_prod);
  }

  return err;
}

char *SrtGrfn3DFXBetaMc(
    char *und3dfx,
    double (*vol_ln_func_Beta)(double t, double S, double F, double *param),
    double *param, int numeventdates, long *eventdates, long tableauRows,
    long *tableauCols, char ***tableauStrings, int **tableauMask, long auxWidth,
    long *auxLen, double **aux, long num_paths, double max_time, int do_pecs,
    double ***prod_val) {
  int free_str = 0;
  FIRSTAllMkts xStr;
  SrtGrfnParam defParm;
  int forback;
  int flag = 0;
  long nstp;

  double *time = NULL, *date = NULL;

  double *dom_ifr = NULL, *dom_vol = NULL, *dom_phi = NULL,

         *for_ifr = NULL, *for_vol = NULL, *for_phi = NULL,

         *fx_vol = NULL,

         *corr_dom_for = NULL, *corr_dom_fx = NULL, *corr_for_fx = NULL;

  void **void_prm = NULL;
  GRFNPARMMC grfn_prm;

  long today, spot_date;
  int i, j, k, num_col;
  SrtUndPtr fx_und, dom_und, for_und;
  TermStruct *fx_ts, *dom_ts, *for_ts;
  char *domname, *forname;
  double dom_lam, for_lam;
  double spot_fx;
  char *dom_yc, *for_yc;
  int fx_idx, dom_idx, for_idx;

  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *sigma_date_fx = NULL,
         *sigma_fx = NULL, *dom_for_cov = NULL, *dom_fx_cov = NULL,
         *for_fx_cov = NULL, *correl_mat = NULL, *correl_dom_for = NULL,
         *correl_dom_fx = NULL, *correl_for_fx = NULL;

  double *merge_dates = NULL, *sig_dom = NULL, *sig_for = NULL, *sig_fx = NULL;

  long nb_merge_dates;

  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx, nb_correl;

  double *beta_tab;
  int *has_evt = NULL;

  long num_stp;

  clock_t t1, t2;

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

  if (get_underlying_type(fx_und) != FOREX_UND) {
    err = serror("Underlying %s is not of type FX", und3dfx);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_fxund(fx_und) != FX_STOCH_RATES) {
    err = serror("Underlying %s is not of type FX Stoch Rates", und3dfx);
    goto FREE_RETURN;
  }

  fx_ts = get_ts_from_fxund(fx_und);

  domname = get_domname_from_fxund(fx_und);
  dom_und = lookup_und(domname);
  if (!dom_und) {
    err = serror("Couldn't find underlying named %s", domname);
    goto FREE_RETURN;
  }
  dom_ts = get_ts_from_irund(dom_und);

  forname = get_forname_from_fxund(fx_und);
  for_und = lookup_und(forname);
  if (!for_und) {
    err = serror("Couldn't find underlying named %s", forname);
    goto FREE_RETURN;
  }
  for_ts = get_ts_from_irund(for_und);

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

  num_stp = (long)(time[nstp - 1] / max_time + 1.0E-08);

  /* fill time vector */
  err = fill_time_vector(&time, &nstp, 0, NULL, 0, NULL, num_stp);

  date = (double *)calloc(nstp, sizeof(double));
  has_evt = (int *)calloc(nstp, sizeof(int));
  beta_tab = (double *)calloc(nstp, sizeof(double));

  if (!date || !has_evt || !beta_tab) {
    err = "Memory allocation error (3) in SrtGrfn3DFXBetaTree";
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

  /* Get all the term structures */
  err = Get_FX_StochRate_TermStructures_corr(
      und3dfx, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_mat, &correl_dom_for, &correl_dom_fx, &correl_for_fx,
      &nb_correl);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &dom_lam);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &for_lam);
  if (err) {
    goto FREE_RETURN;
  }

  /*	Get Fx spot and yield curves */
  dom_yc = get_ycname_from_irund(dom_und);
  for_yc = get_ycname_from_irund(for_und);

  spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  spot_fx = get_spot_from_fxund(fx_und) * swp_f_df(today, spot_date, dom_yc) /
            swp_f_df(today, spot_date, for_yc);

  /* now merge all the term structure */
  merge_dates = (double *)calloc(sigma_n_dom, sizeof(double));

  if (!merge_dates) {
    err = "Memory allocation error (3) in SrtGrfn3DFXBetaMc";
    goto FREE_RETURN;
  }

  memcpy(merge_dates, sigma_date_dom, sigma_n_dom * sizeof(double));
  nb_merge_dates = sigma_n_dom;
  num_f_concat_vector(&nb_merge_dates, &merge_dates, sigma_n_for,
                      sigma_date_for);
  num_f_concat_vector(&nb_merge_dates, &merge_dates, sigma_n_fx, sigma_date_fx);
  num_f_concat_vector(&nb_merge_dates, &merge_dates, nb_correl, correl_mat);
  num_f_sort_vector(nb_merge_dates, merge_dates);
  num_f_unique_vector(&nb_merge_dates, merge_dates);

  /*	Fill the new term structures */

  sig_dom = (double *)calloc(nb_merge_dates, sizeof(double));
  sig_for = (double *)calloc(nb_merge_dates, sizeof(double));
  sig_fx = (double *)calloc(nb_merge_dates, sizeof(double));

  dom_vol = (double *)calloc(nstp, sizeof(double));
  for_vol = (double *)calloc(nstp, sizeof(double));
  fx_vol = (double *)calloc(nstp, sizeof(double));

  dom_ifr = (double *)calloc(nstp, sizeof(double));
  dom_phi = (double *)calloc(nstp, sizeof(double));
  for_ifr = (double *)calloc(nstp, sizeof(double));
  for_phi = (double *)calloc(nstp, sizeof(double));

  corr_dom_for = (double *)calloc(nstp, sizeof(double));
  corr_dom_fx = (double *)calloc(nstp, sizeof(double));
  corr_for_fx = (double *)calloc(nstp, sizeof(double));

  if (!sig_dom || !sig_for || !sig_fx || !dom_vol || !for_vol || !fx_vol ||
      !corr_dom_for || !corr_dom_fx || !corr_for_fx || !dom_ifr || !dom_phi ||
      !for_ifr || !for_phi) {
    err = "Memory allocation error (4) in SrtGrfn3DFXBetaMc";
    goto FREE_RETURN;
  }

  for (i = nb_merge_dates - 1; i >= 0; i--) {
    sig_dom[i] = find_sig(merge_dates[i], dom_ts);
    sig_for[i] = find_sig(merge_dates[i], for_ts);
    sig_fx[i] = find_fx_sig(merge_dates[i], fx_ts);
  }

  /*	Get phi        , std... */
  err = fill_Betamc_init(date, time, nstp, merge_dates, nb_merge_dates, sig_dom,
                         dom_lam, sig_for, for_lam, sig_fx, dom_yc, for_yc,
                         dom_ifr, dom_vol, dom_phi, for_ifr, for_vol, for_phi,
                         fx_vol);

  if (err) {
    goto FREE_RETURN;
  }

  for (i = nstp - 1; i >= 0; i--) {
    corr_dom_for[i] = correl_dom_for[Get_Index(time[i], correl_mat, nb_correl)];
    corr_dom_fx[i] = correl_dom_fx[Get_Index(time[i], correl_mat, nb_correl)];
    corr_for_fx[i] = correl_for_fx[Get_Index(time[i], correl_mat, nb_correl)];
  }

  /* Check if there is an event */
  j = 0;
  for (i = 0; (i < nstp && j < xStr.num_evt); i++) {
    if (time[i] == xStr.tms[j]) {
      has_evt[i] = 1;
      j += 1;
    } else if (j > 0 && xStr.am[j - 1]) {
      has_evt[i] = 1;
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

  /* get the index of each underlying */
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

  void_prm = (void **)calloc(nstp, sizeof(void *));

  if (!void_prm) {
    err = "Memory allocation error (6) in SrtGrfn3DFXMc";
    goto FREE_RETURN;
  }

  j = xStr.num_evt - 1;

  for (i = nstp - 1; i >= 0; i--) {
    if (j >= 0 && fabs(date[i] - xStr.dts[j]) < 1.0e-04) {
      grfn_prm = malloc(sizeof(grfn_parm_mc));
      grfn_prm->global = &xStr;
      grfn_prm->local = xStr.evt + j;
      grfn_prm->fx_idx = fx_idx;
      grfn_prm->dom_idx = dom_idx;
      grfn_prm->for_idx = for_idx;

      grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
      grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
      grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

      if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1) {
        grfn_prm->fx_dff = dvector(0, grfn_prm->num_fx_df - 1);
        grfn_prm->fx_gam = dvector(0, grfn_prm->num_fx_df - 1);
        grfn_prm->fx_gam2 = dvector(0, grfn_prm->num_fx_df - 1);

        if (!grfn_prm->fx_dff || !grfn_prm->fx_gam || !grfn_prm->fx_gam2) {
          err = "Memory allocation error (7) in SrtGrfn3DFXBetaMc";
          goto FREE_RETURN;
        }

        for (k = 0; k < grfn_prm->num_fx_df; k++) {
          grfn_prm->fx_dff[k] =
              swp_f_df(xStr.dts[j], grfn_prm->fx_df_dts[k], (char *)dom_yc);
          grfn_prm->fx_gam[k] =
              (1.0 - exp(-dom_lam * grfn_prm->fx_df_tms[k])) / dom_lam;
          grfn_prm->fx_gam2[k] =
              0.5 * grfn_prm->fx_gam[k] * grfn_prm->fx_gam[k] * dom_phi[i];
        }

        grfn_prm->do_fx = 1;
      } else {
        grfn_prm->do_fx = 0;
      }

      grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
      grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
      grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

      if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1) {
        grfn_prm->dom_dff = dvector(0, grfn_prm->num_dom_df - 1);
        grfn_prm->dom_gam = dvector(0, grfn_prm->num_dom_df - 1);
        grfn_prm->dom_gam2 = dvector(0, grfn_prm->num_dom_df - 1);

        if (!grfn_prm->dom_dff || !grfn_prm->dom_gam || !grfn_prm->dom_gam2) {
          err = "Memory allocation error (8) in SrtGrfn3DFXMc";
          goto FREE_RETURN;
        }

        for (k = 0; k < grfn_prm->num_dom_df; k++) {
          grfn_prm->dom_dff[k] =
              swp_f_df(xStr.dts[j], grfn_prm->dom_df_dts[k], (char *)dom_yc);
          grfn_prm->dom_gam[k] =
              (1.0 - exp(-dom_lam * grfn_prm->dom_df_tms[k])) / dom_lam;
          grfn_prm->dom_gam2[k] =
              0.5 * grfn_prm->dom_gam[k] * grfn_prm->dom_gam[k] * dom_phi[i];
        }

        grfn_prm->do_dom = 1;
      } else {
        grfn_prm->do_dom = 0;
      }

      grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
      grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
      grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

      if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1) {
        grfn_prm->for_dff = dvector(0, grfn_prm->num_for_df - 1);
        grfn_prm->for_gam = dvector(0, grfn_prm->num_for_df - 1);
        grfn_prm->for_gam2 = dvector(0, grfn_prm->num_for_df - 1);

        if (!grfn_prm->for_dff || !grfn_prm->for_gam || !grfn_prm->for_gam2) {
          err = "Memory allocation error (9) in SrtGrfn3DFXMc";
          goto FREE_RETURN;
        }

        for (k = 0; k < grfn_prm->num_for_df; k++) {
          grfn_prm->for_dff[k] =
              swp_f_df(xStr.dts[j], grfn_prm->for_df_dts[k], (char *)for_yc);
          grfn_prm->for_gam[k] =
              (1.0 - exp(-for_lam * grfn_prm->for_df_tms[k])) / for_lam;
          grfn_prm->for_gam2[k] =
              0.5 * grfn_prm->for_gam[k] * grfn_prm->for_gam[k] * for_phi[i];
        }

        grfn_prm->do_for = 1;
      } else {
        grfn_prm->do_for = 0;
      }

      j--;
      while (j >= 0 && xStr.evt[j].evt == NULL) {
        j--;
      }

      void_prm[i + flag] = (void *)grfn_prm;
      has_evt[i + flag] = 1;
    } else if (j >= 0 && xStr.am[j]) {
      grfn_prm = malloc(sizeof(grfn_parm_mc));
      grfn_prm->global = &xStr;
      grfn_prm->local = xStr.evt + j;
      grfn_prm->fx_idx = fx_idx;
      grfn_prm->dom_idx = dom_idx;
      grfn_prm->for_idx = for_idx;

      grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
      grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
      grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

      if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1) {
        grfn_prm->fx_dff = dvector(0, grfn_prm->num_fx_df - 1);
        grfn_prm->fx_gam = dvector(0, grfn_prm->num_fx_df - 1);
        grfn_prm->fx_gam2 = dvector(0, grfn_prm->num_fx_df - 1);

        if (!grfn_prm->fx_dff || !grfn_prm->fx_gam || !grfn_prm->fx_gam2) {
          err = "Memory allocation error (7) in SrtGrfn3DFXBetaMc";
          goto FREE_RETURN;
        }

        for (k = 0; k < grfn_prm->num_fx_df; k++) {
          grfn_prm->fx_dff[k] =
              swp_f_df(date[i], grfn_prm->fx_df_dts[k], (char *)dom_yc);
          grfn_prm->fx_gam[k] =
              (1.0 - exp(-dom_lam * (grfn_prm->fx_df_dts[k] - date[i]) *
                         YEARS_IN_DAY)) /
              dom_lam;
          grfn_prm->fx_gam2[k] =
              0.5 * grfn_prm->fx_gam[k] * grfn_prm->fx_gam[k] * dom_phi[i];
        }

        grfn_prm->do_fx = 1;
      } else {
        grfn_prm->do_fx = 0;
      }

      grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
      grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
      grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

      if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1) {
        grfn_prm->dom_dff = dvector(0, grfn_prm->num_dom_df - 1);
        grfn_prm->dom_gam = dvector(0, grfn_prm->num_dom_df - 1);
        grfn_prm->dom_gam2 = dvector(0, grfn_prm->num_dom_df - 1);

        if (!grfn_prm->dom_dff || !grfn_prm->dom_gam || !grfn_prm->dom_gam2) {
          err = "Memory allocation error (8) in SrtGrfn3DFXMc";
          goto FREE_RETURN;
        }

        for (k = 0; k < grfn_prm->num_dom_df; k++) {
          grfn_prm->dom_dff[k] =
              swp_f_df(date[i], grfn_prm->dom_df_dts[k], (char *)dom_yc);
          grfn_prm->dom_gam[k] =
              (1.0 - exp(-dom_lam * (grfn_prm->dom_df_dts[k] - date[i]) *
                         YEARS_IN_DAY)) /
              dom_lam;
          grfn_prm->dom_gam2[k] =
              0.5 * grfn_prm->dom_gam[k] * grfn_prm->dom_gam[k] * dom_phi[i];
        }

        grfn_prm->do_dom = 1;
      } else {
        grfn_prm->do_dom = 0;
      }

      grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
      grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
      grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

      if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1) {
        grfn_prm->for_dff = dvector(0, grfn_prm->num_for_df - 1);
        grfn_prm->for_gam = dvector(0, grfn_prm->num_for_df - 1);
        grfn_prm->for_gam2 = dvector(0, grfn_prm->num_for_df - 1);

        if (!grfn_prm->for_dff || !grfn_prm->for_gam || !grfn_prm->for_gam2) {
          err = "Memory allocation error (9) in SrtGrfn3DFXMc";
          goto FREE_RETURN;
        }

        for (k = 0; k < grfn_prm->num_for_df; k++) {
          grfn_prm->for_dff[k] =
              swp_f_df(date[i], grfn_prm->for_df_dts[k], (char *)for_yc);
          grfn_prm->for_gam[k] =
              (1.0 - exp(-for_lam * (grfn_prm->for_df_dts[k] - date[i]) *
                         YEARS_IN_DAY)) /
              for_lam;
          grfn_prm->for_gam2[k] =
              0.5 * grfn_prm->for_gam[k] * grfn_prm->for_gam[k] * for_phi[i];
        }

        grfn_prm->do_for = 1;
      } else {
        grfn_prm->do_for = 0;
      }

      void_prm[i + flag] = (void *)grfn_prm;
      has_evt[i + flag] = 1;
    } else {
      void_prm[i + flag] = NULL;
      has_evt[i + flag] = 0;
    }
  }

  /*	Eventually! call to function */

  *prod_val = dmatrix(0, num_col - 1, 0, 1);

  t2 = clock();

  smessage("Phase 1 -preprocessing        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

  err = mcLocal_main_3dfx2(
      /*	Time data */
      num_paths, nstp, num_col, time, date, dom_ifr, /*	Distributions */
      dom_vol, dom_phi, for_ifr, for_vol, for_phi, fx_vol, fx_vol,
      vol_ln_func_Beta, param, void_prm, has_evt,
      /*	Model data */
      dom_lam, for_lam, corr_dom_for, corr_dom_fx, corr_for_fx,
      /*	Market data */
      spot_fx, dom_yc, for_yc,
      /* Do PECS adjustment */
      do_pecs,
      /*	Payoff function */
      grfn_payoff_4_3dfx_Betamc, /*	Result */
      *prod_val);

  *tableauCols = num_col;

  /*	Add PV of Past */
  (*prod_val)[num_col - 1][0] += xStr.gd->pv_of_past;

FREE_RETURN:

  if (time)
    free(time);
  if (date)
    free(date);

  if (dom_ifr)
    free(dom_ifr);
  if (dom_vol)
    free(dom_vol);
  if (dom_phi)
    free(dom_phi);

  if (for_ifr)
    free(for_ifr);
  if (for_vol)
    free(for_vol);
  if (for_phi)
    free(for_phi);

  if (fx_vol)
    free(fx_vol);

  if (sigma_date_dom)
    free(sigma_date_dom);
  if (sigma_dom)
    free(sigma_dom);
  if (tau_date_dom)
    free(tau_date_dom);
  if (tau_dom)
    free(tau_dom);

  if (sigma_date_for)
    free(sigma_date_for);
  if (sigma_for)
    free(sigma_for);
  if (tau_date_for)
    free(tau_date_for);
  if (tau_for)
    free(tau_for);

  if (sigma_date_fx)
    free(sigma_date_fx);
  if (sigma_fx)
    free(sigma_fx);

  if (correl_mat)
    free(correl_mat);
  if (correl_dom_for)
    free(correl_dom_for);
  if (correl_dom_fx)
    free(correl_dom_fx);
  if (correl_for_fx)
    free(correl_for_fx);

  if (corr_dom_for)
    free(corr_dom_for);
  if (corr_dom_fx)
    free(corr_dom_fx);
  if (corr_for_fx)
    free(corr_for_fx);

  if (dom_for_cov)
    free(dom_for_cov);
  if (dom_fx_cov)
    free(dom_fx_cov);
  if (for_fx_cov)
    free(for_fx_cov);

  if (beta_tab)
    free(beta_tab);

  if (void_prm) {
    for (i = 0; i < nstp; i++) {
      if (void_prm[i]) {
        grfn_prm = (GRFNPARMMC)void_prm[i];

        if (grfn_prm->do_fx && grfn_prm->num_fx_df > 0 &&
            grfn_prm->fx_idx != -1) {
          if (grfn_prm->fx_dff)
            free_dvector(grfn_prm->fx_dff, 0, grfn_prm->num_fx_df - 1);
          if (grfn_prm->fx_gam)
            free_dvector(grfn_prm->fx_gam, 0, grfn_prm->num_fx_df - 1);
          if (grfn_prm->fx_gam2)
            free_dvector(grfn_prm->fx_gam2, 0, grfn_prm->num_fx_df - 1);
        }

        if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0 &&
            grfn_prm->dom_idx != -1) {
          if (grfn_prm->dom_dff)
            free_dvector(grfn_prm->dom_dff, 0, grfn_prm->num_dom_df - 1);
          if (grfn_prm->dom_gam)
            free_dvector(grfn_prm->dom_gam, 0, grfn_prm->num_dom_df - 1);
          if (grfn_prm->dom_gam2)
            free_dvector(grfn_prm->dom_gam2, 0, grfn_prm->num_dom_df - 1);
        }

        if (grfn_prm->do_for && grfn_prm->num_for_df > 0 &&
            grfn_prm->for_idx != -1) {
          if (grfn_prm->for_dff)
            free_dvector(grfn_prm->for_dff, 0, grfn_prm->num_for_df - 1);
          if (grfn_prm->for_gam)
            free_dvector(grfn_prm->for_gam, 0, grfn_prm->num_for_df - 1);
          if (grfn_prm->for_gam2)
            free_dvector(grfn_prm->for_gam2, 0, grfn_prm->num_for_df - 1);
        }

        free(grfn_prm);
      }
    }

    free(void_prm);
  }

  if (free_str) {
    FIRSTFreeMktStruct(&xStr);
  }

  return err;
}