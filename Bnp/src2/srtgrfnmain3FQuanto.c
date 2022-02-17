/**********************************************************************
 *      Name: SrtGrfnMain3FQuanto.c                                   *
 *  Function: Entry point to GRFN with raw data                       *
 *            in the case of the 3Factor QUANTO model * Copyright: (C) BNP
 *Paribas Capital Markets Ltd.                    *
 *--------------------------------------------------------------------*
 *    Author: Jean-Michel LY	                                      *
 *      Date: 22/10/2003                                              *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere (Excel or 2020)                  *
 *   Returns:                                                         *
 *   Globals: Expects mkt and request list structures to exist        *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 **********************************************************************/

#include "SrtAccess.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_all3FQuanto.h"

char *SrtGrfn3FQuantoTree(char *und3Fquanto, double *correl_dates,
                          int n_correl_dates, double ***correl_matrix_4x4_ts,
                          int numeventdates, long *eventdates, long tableauRows,
                          long tableauCols, char ***tableauStrings,
                          int **tableauMask, long auxWidth, long *auxLen,
                          double **aux, long num_stp, int *num_prod,
                          int discount, double **prod_val) {
  int free_str = 0;
  FIRSTAllMkts xStr;
  SrtGrfnParam defParm;
  int forback;
  long nstp;

  double *time = NULL, *date = NULL;

  int *vol_change = NULL;

  double *r1_sigma_dates = NULL, *r1_sigma = NULL, *r3_sigma_dates = NULL,
         *r3_sigma = NULL, *merge_dates = NULL, *r1_merge_sigma = NULL,
         *r2_merge_sigma = NULL, *r3_merge_sigma = NULL, *fx_merge_sigma = NULL,
         *corr_merge_r1_r2 = NULL, *corr_merge_r1_r3 = NULL,
         *corr_merge_r1_fx = NULL, *corr_merge_r2_r3 = NULL,
         *corr_merge_r2_fx = NULL, *corr_merge_r3_fx = NULL;

  double *lgm2F_ifr = NULL, *r1_fwd = NULL, *r1_var = NULL, *r2_fwd = NULL,
         *r2_var = NULL, *r1_r2_covar = NULL, *lgm1F_ifr = NULL, *r3_fwd = NULL,
         *r3_var = NULL;

  void **void_prm = NULL;
  GRFNPARMTREE grfn_prm;
  int *is_event = NULL;

  long today, spot_date;
  int i, j;
  SrtUndPtr fx_und, dom_und, for_und;

  char *dom_name = NULL, *for_name = NULL;

  // LGM2F Parameters
  double alpha, gamma, rho;

  char *lgm2F_yc, *lgm1F_yc;

  int fx_idx, dom_idx, for_idx;

  int num_vol_times;
  double *vol_times = NULL;

  clock_t t1, t2;

  // To get term structure details of the 3F quanto underlying
  double *dom_sigma_dates = NULL, *dom_sigma_values = NULL,
         *for_sigma_dates = NULL, *for_sigma_values = NULL,
         *fx_sigma_dates = NULL, *fx_sigma_values = NULL,
         *dom_sigma_merge_values = NULL, *for_sigma_merge_values = NULL;

  int lgm2F_is_dom_for, correl_index;
  long n_dom_sigma_dates, n_for_sigma_dates, n_fx_sigma_dates;
  long n_r1_sigma_dates, n_r3_sigma_dates;
  double dom_fixed_lambda, dom_fixed_alpha, dom_fixed_gamma, dom_fixed_rho;
  double for_fixed_lambda, for_fixed_alpha, for_fixed_gamma, for_fixed_rho;
  double r1_lambda, r2_lambda, r3_lambda;

  Err err = NULL;

  t1 = clock();

  // Initialise the GRFN tableau

  // Initialise the param struct

  err = srt_f_set_default_GrfnParams(&defParm);
  defParm.min_nodes_per_path = num_stp;

  err = FIRSTInitMktStruct(numeventdates, eventdates, tableauRows, tableauCols,
                           tableauStrings, tableauMask, auxWidth, auxLen, aux,
                           und3Fquanto, &defParm, &forback, &xStr);

  if (err) {
    goto FREE_RETURN;
  }

  if (xStr.num_und > 3) {
    err = "Fx3FQuantoGrfn can only accept 3 underlyings  , please check the "
          "GRFN tableau or use the MultiGRFN functions";
    goto FREE_RETURN;
  }

  free_str = 1;

  // Now  , lookup underlyings involved and their term structures

  fx_und = lookup_und(und3Fquanto);

  if (!fx_und) {
    err = serror("Couldn't find underlying named %s", und3Fquanto);
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(fx_und);

  if (get_underlying_type(fx_und) != FOREX_UND) {
    err = serror("Underlying %s is not of type FX", und3Fquanto);
    goto FREE_RETURN;
  }

  if (get_mdltype_from_fxund(fx_und) != FX_STOCH_RATES) {
    err = serror("Underlying %s is not of type FX Stoch Rates", und3Fquanto);
    goto FREE_RETURN;
  }

  // Get the underlying parameters of the 3F Quanto model
  dom_name = (char *)calloc(256, sizeof(char));
  for_name = (char *)calloc(256, sizeof(char));

  err = Get_FX_StochRate_TermStructures3FQuanto(
      und3Fquanto, &dom_name, &dom_sigma_dates, &n_dom_sigma_dates,
      &dom_sigma_values, &dom_fixed_lambda, &dom_fixed_alpha, &dom_fixed_gamma,
      &dom_fixed_rho, &for_name, &for_sigma_dates, &n_for_sigma_dates,
      &for_sigma_values, &for_fixed_lambda, &for_fixed_alpha, &for_fixed_gamma,
      &for_fixed_rho, &lgm2F_is_dom_for, &fx_sigma_dates, &n_fx_sigma_dates,
      &fx_sigma_values);
  if (err)
    goto FREE_RETURN;

  // Locate the underlyings
  dom_und = lookup_und(dom_name);
  if (!dom_und) {
    err = serror("Couldn't find underlying named %s", dom_name);
    goto FREE_RETURN;
  }

  for_und = lookup_und(for_name);
  if (!for_und) {
    err = serror("Couldn't find underlying named %s", for_name);
    goto FREE_RETURN;
  }

  // Differentiate between r1 and r3
  if (lgm2F_is_dom_for == 0) {
    // LGM2F is DOMESTIC
    alpha = dom_fixed_alpha;
    gamma = dom_fixed_gamma;
    rho = dom_fixed_rho;
    r1_lambda = dom_fixed_lambda;
    r2_lambda = r1_lambda + gamma;
    r1_sigma_dates = dom_sigma_dates;
    r1_sigma = dom_sigma_values;
    n_r1_sigma_dates = n_dom_sigma_dates;
    lgm2F_yc = get_ycname_from_irund(dom_und);

    // LGM1F is FOREIGN
    r3_lambda = for_fixed_lambda;
    r3_sigma_dates = for_sigma_dates;
    r3_sigma = for_sigma_values;
    n_r3_sigma_dates = n_for_sigma_dates;
    lgm1F_yc = get_ycname_from_irund(for_und);

  } else {
    // LGM1F is DOMESTIC
    r3_lambda = dom_fixed_lambda;
    r3_sigma_dates = dom_sigma_dates;
    r3_sigma = dom_sigma_values;
    n_r3_sigma_dates = n_dom_sigma_dates;
    lgm1F_yc = get_ycname_from_irund(dom_und);

    // LGM2F is FOREIGN
    alpha = for_fixed_alpha;
    gamma = for_fixed_gamma;
    rho = for_fixed_rho;
    r1_lambda = for_fixed_lambda;
    r2_lambda = r1_lambda + gamma;
    r1_sigma_dates = for_sigma_dates;
    r1_sigma = for_sigma_values;
    n_r1_sigma_dates = n_for_sigma_dates;
    lgm2F_yc = get_ycname_from_irund(for_und);
  }

  // dom_yc = get_ycname_from_irund (dom_und);
  // for_yc = get_ycname_from_irund (for_und);

  spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  // Copy event dates in array "time"
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
    err = "Memory allocation error (1) in SrtGrfn3FQuantoTree";
    goto FREE_RETURN;
  }
  memcpy(time, xStr.tms, nstp * sizeof(double));

  // Create an array containing all dates of r1_sigma  , r3_sigma and fx_sigma
  // which are <= the last event date
  err = compute_vol_times_3FQuanto(
      r1_sigma_dates, r1_sigma, n_r1_sigma_dates, r3_sigma_dates, r3_sigma,
      n_r3_sigma_dates, fx_sigma_dates, fx_sigma_values, n_fx_sigma_dates,
      &num_vol_times, &vol_times, time[nstp - 1]);
  if (err)
    goto FREE_RETURN;

  // Takes the time vector containing the event dates and slices it into a
  // bigger time vector using the number of steps required in the tree. Don't
  // know why this function needs bar_times and vol_times since they are not
  // used....
  err = fill_time_vector(&time, &nstp, 0, (double *)NULL, num_vol_times,
                         vol_times, num_stp);
  if (err)
    goto FREE_RETURN;

  date = (double *)calloc(nstp, sizeof(double));
  if (!date) {
    err = "Memory allocation error (3) in SrtGrfn3FQuantoTree";
    goto FREE_RETURN;
  }

  for (i = 0; i < nstp; i++) {
    date[i] = today + DAYS_IN_YEAR * time[i];

    if (i > 0 && date[i] - date[i - 1] >= 1) {
      date[i] = (long)(date[i] + 1.0e-08);
      time[i] = YEARS_IN_DAY * (date[i] - today);
    }
  }

  // Allocates memory for vols and all correlations TS
  // The KEY point is that we will project all the term structures (vols  ,
  // correls) on the single "time" schedule which contains the GRFN event dates
  // and intermediary dates which are only calculated from the number of steps
  // we want. (no vols or corr dates) We MAY therefore be WRONG on the
  // vols/corrs since the vols/corr dates will not always match the "time"
  // schedule but we do not care since increasing the number of steps will lead
  // to convergence.

  vol_change = (int *)calloc(nstp, sizeof(int));
  r1_merge_sigma = (double *)calloc(nstp, sizeof(double));
  r2_merge_sigma = (double *)calloc(nstp, sizeof(double));
  r3_merge_sigma = (double *)calloc(nstp, sizeof(double));
  fx_merge_sigma = (double *)calloc(nstp, sizeof(double));

  corr_merge_r1_r2 = (double *)calloc(nstp, sizeof(double));
  corr_merge_r1_r3 = (double *)calloc(nstp, sizeof(double));
  corr_merge_r1_fx = (double *)calloc(nstp, sizeof(double));
  corr_merge_r2_r3 = (double *)calloc(nstp, sizeof(double));
  corr_merge_r2_fx = (double *)calloc(nstp, sizeof(double));
  corr_merge_r3_fx = (double *)calloc(nstp, sizeof(double));

  if (!vol_change || !r1_merge_sigma || !r2_merge_sigma || !r3_merge_sigma ||
      !fx_merge_sigma || !corr_merge_r1_r2 || !corr_merge_r1_r3 ||
      !corr_merge_r1_fx || !corr_merge_r2_r3 || !corr_merge_r2_fx ||
      !corr_merge_r3_fx) {
    err = "Memory allocation error (4) in SrtGrfn3FQuantoTree";
    goto FREE_RETURN;
  }

  // Filling for the last event date time[nstp-1]
  vol_change[nstp - 1] = 1;
  r1_merge_sigma[nstp - 1] =
      r1_sigma[Get_Index(time[nstp - 1], r1_sigma_dates, n_r1_sigma_dates)];
  r2_merge_sigma[nstp - 1] = alpha * r1_merge_sigma[nstp - 1];
  r3_merge_sigma[nstp - 1] =
      r3_sigma[Get_Index(time[nstp - 1], r3_sigma_dates, n_r3_sigma_dates)];
  fx_merge_sigma[nstp - 1] = fx_sigma_values[Get_Index(
      time[nstp - 1], fx_sigma_dates, n_fx_sigma_dates)];

  correl_index = Get_Index(time[nstp - 1], correl_dates, n_correl_dates);

  corr_merge_r1_r2[nstp - 1] = rho;
  corr_merge_r1_r3[nstp - 1] =
      correl_matrix_4x4_ts[correl_index][LGM2F_W1][LGM1F_W3];
  corr_merge_r2_r3[nstp - 1] =
      correl_matrix_4x4_ts[correl_index][LGM2F_W2][LGM1F_W3];

  if (lgm2F_is_dom_for == 0) {
    corr_merge_r1_fx[nstp - 1] = 0.0;
    corr_merge_r2_fx[nstp - 1] = 0.0;
    corr_merge_r3_fx[nstp - 1] =
        correl_matrix_4x4_ts[correl_index][LGM1F_W3][FX_3F_QUANTO];
  } else {
    corr_merge_r1_fx[nstp - 1] =
        correl_matrix_4x4_ts[correl_index][LGM2F_W1][FX_3F_QUANTO];
    corr_merge_r2_fx[nstp - 1] =
        correl_matrix_4x4_ts[correl_index][LGM2F_W2][FX_3F_QUANTO];
    corr_merge_r3_fx[nstp - 1] = 0.0;
  }

  // Filling for all previous dates in the master schedule time[i]
  for (i = nstp - 2; i >= 0; i--) {
    r1_merge_sigma[i] =
        r1_sigma[Get_Index(time[i], r1_sigma_dates, n_r1_sigma_dates)];
    r2_merge_sigma[i] = alpha * r1_merge_sigma[i];
    r3_merge_sigma[i] =
        r3_sigma[Get_Index(time[i], r3_sigma_dates, n_r3_sigma_dates)];
    fx_merge_sigma[i] =
        fx_sigma_values[Get_Index(time[i], fx_sigma_dates, n_fx_sigma_dates)];

    correl_index = Get_Index(time[i], correl_dates, n_correl_dates);

    corr_merge_r1_r2[i] = rho;
    corr_merge_r1_r3[i] =
        correl_matrix_4x4_ts[correl_index][LGM2F_W1][LGM1F_W3];
    corr_merge_r2_r3[i] =
        correl_matrix_4x4_ts[correl_index][LGM2F_W2][LGM1F_W3];

    if (lgm2F_is_dom_for == 0) {
      corr_merge_r1_fx[i] = 0.0;
      corr_merge_r2_fx[i] = 0.0;
      corr_merge_r3_fx[i] =
          correl_matrix_4x4_ts[correl_index][LGM1F_W3][FX_3F_QUANTO];
    } else {
      corr_merge_r1_fx[i] =
          correl_matrix_4x4_ts[correl_index][LGM2F_W1][FX_3F_QUANTO];
      corr_merge_r2_fx[i] =
          correl_matrix_4x4_ts[correl_index][LGM2F_W2][FX_3F_QUANTO];
      corr_merge_r3_fx[i] = 0.0;
    }

    if (fabs(r1_merge_sigma[i] - r1_merge_sigma[i + 1]) +
            fabs(r2_merge_sigma[i] - r2_merge_sigma[i + 1]) +
            fabs(r3_merge_sigma[i] - r3_merge_sigma[i + 1]) +
            fabs(fx_merge_sigma[i] - fx_merge_sigma[i + 1]) +
            fabs(corr_merge_r1_r2[i] - corr_merge_r1_r2[i + 1]) +
            fabs(corr_merge_r1_r3[i] - corr_merge_r1_r3[i + 1]) +
            fabs(corr_merge_r1_fx[i] - corr_merge_r1_fx[i + 1]) +
            fabs(corr_merge_r2_r3[i] - corr_merge_r2_r3[i + 1]) +
            fabs(corr_merge_r2_fx[i] - corr_merge_r2_fx[i + 1]) +
            fabs(corr_merge_r3_fx[i] - corr_merge_r3_fx[i + 1]) >
        EPS) {
      vol_change[i] = 1;
    } else {
      vol_change[i] = 0;
    }
  }

  /*	Get distributions */

  lgm2F_ifr = (double *)calloc(nstp, sizeof(double));
  r1_fwd = (double *)calloc(nstp, sizeof(double));
  r1_var = (double *)calloc(nstp, sizeof(double));
  r2_fwd = (double *)calloc(nstp, sizeof(double));
  r2_var = (double *)calloc(nstp, sizeof(double));
  r1_r2_covar = (double *)calloc(nstp, sizeof(double));
  lgm1F_ifr = (double *)calloc(nstp, sizeof(double));
  r3_fwd = (double *)calloc(nstp, sizeof(double));
  r3_var = (double *)calloc(nstp, sizeof(double));

  if (!lgm2F_ifr || !r1_fwd || !r1_var || !r2_fwd || !r2_var || !r1_r2_covar ||
      !lgm1F_ifr || !r3_fwd || !r3_var) {
    err = "Memory allocation error (5) in SrtGrfn3FQuantoTree";
    goto FREE_RETURN;
  }

  // Now fills the forwards  , variances and ifr
  fill_fwd_var_corr_3FQuanto(
      nstp, time, date, r1_merge_sigma, r2_merge_sigma, r3_merge_sigma,
      fx_merge_sigma, r1_lambda, r2_lambda, r3_lambda, corr_merge_r1_r2,
      corr_merge_r1_r3, corr_merge_r1_fx, corr_merge_r2_r3, corr_merge_r2_fx,
      corr_merge_r3_fx, lgm2F_yc, lgm1F_yc, lgm2F_is_dom_for, r1_fwd, r1_var,
      r2_fwd, r2_var, r1_r2_covar, r3_fwd, r3_var, lgm2F_ifr, lgm1F_ifr);

  /*	Fill product structure */

  strupper(und3Fquanto);
  strip_white_space(und3Fquanto);
  strupper(dom_name);
  strip_white_space(dom_name);
  strupper(for_name);
  strip_white_space(for_name);
  for (i = 0; i < xStr.num_und; i++) {
    strupper(xStr.und_data[i].und_name);
    strip_white_space(xStr.und_data[i].und_name);
  }

  fx_idx = -1;
  for (i = 0; i < xStr.num_und; i++) {
    if (!strcmp(xStr.und_data[i].und_name, und3Fquanto)) {
      fx_idx = i;
    }
  }
  if (fx_idx == -1) {
    err = "The Fx underlying is not present in the mdlcomm structure";
    goto FREE_RETURN;
  }

  dom_idx = -1;
  for (i = 0; i < xStr.num_und; i++) {
    if (!strcmp(xStr.und_data[i].und_name, dom_name)) {
      dom_idx = i;
    }
  }
  if (dom_idx == -1) {
    err = "The domestic underlying is not present in the mdlcomm structure";
    goto FREE_RETURN;
  }

  for_idx = -1;
  for (i = 0; i < xStr.num_und; i++) {
    if (!strcmp(xStr.und_data[i].und_name, for_name)) {
      for_idx = i;
    }
  }
  if (for_idx == -1) {
    err = "The foreign underlying is not present in the mdlcomm structure";
    goto FREE_RETURN;
  }

  is_event = (int *)calloc(nstp, sizeof(int));
  void_prm = (void **)calloc(nstp, sizeof(void *));

  if (!is_event || !void_prm) {
    err = "Memory allocation error (6) in SrtGrfn3FQuantoTree";
    goto FREE_RETURN;
  }

  j = xStr.num_evt - 1;
  for (i = nstp - 1; i >= 0; i--) {

    // Different fillings depending on the nature of date[i]
    if (j >= 0 && fabs(date[i] - xStr.dts[j]) < 1.0e-04) {
      // Filling if date[i] is a GRFN event date (tracked by index j)
      grfn_prm = malloc(sizeof(grfn_parm_tree));
      grfn_prm->global = &xStr;
      grfn_prm->local = xStr.evt + j;
      grfn_prm->fx_idx = fx_idx;
      grfn_prm->dom_idx = dom_idx;
      grfn_prm->for_idx = for_idx;

      // Detects if any DF without underlying is requested at the GRFN date (FX
      // DF = DOM DF)
      grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
      grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
      grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

      // Detects if any DOM DF is requested at the GRFN date
      grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
      grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
      grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

      // Detects if any FOR DF is requested at the GRFN date
      grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
      grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
      grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

      is_event[i] = 1;
      void_prm[i] = (void *)grfn_prm;

      j--;
      while (j >= 0 && xStr.evt[j].evt == NULL) {
        j--;
      }
    } else if (j >= 0 && xStr.am[j]) {
      // Filling if AM is used

      grfn_prm = malloc(sizeof(grfn_parm_tree));
      grfn_prm->global = &xStr;
      grfn_prm->local = xStr.evt + j;
      grfn_prm->fx_idx = fx_idx;
      grfn_prm->dom_idx = dom_idx;
      grfn_prm->for_idx = for_idx;

      // Detects if any DF without underlying is requested at the GRFN date (FX
      // DF = DOM DF)
      grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
      grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
      grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

      // Detects if any DOM DF is requested at the GRFN date
      grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
      grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
      grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

      // Detects if any FOR DF is requested at the GRFN date
      grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
      grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
      grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

      is_event[i] = 1;
      void_prm[i] = (void *)grfn_prm;
    } else {
      // Filling when i is a discretization date which is not a GRFN date (most
      // cases)
      is_event[i] = 0;
      void_prm[i] = NULL;
    }
  }

  // Call to pricing function
  *num_prod = xStr.num_cols;
  *prod_val = (double *)calloc(*num_prod + 1, sizeof(double));

  if (!(*prod_val)) {
    err = "Memory allocation error (7) in SrtGrfn3FQuantoTree";
    goto FREE_RETURN;
  }

  t2 = clock();

  smessage("Phase 1 -preprocessing  , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);
  smessage("Number of times steps  , required: %d  , actual: %d", num_stp,
           nstp);

  err = tree_main_3fquanto(
      nstp, time, date, vol_change,

      // LGM2F state variables r1 and r2
      r1_lambda, r1_merge_sigma, r1_fwd, r1_var, r2_lambda, r2_merge_sigma,
      r2_fwd, r2_var, r1_r2_covar, lgm2F_ifr, lgm2F_yc, lgm2F_is_dom_for,

      // LGM1F state variable r3
      r3_lambda, r3_merge_sigma, r3_fwd, r3_var, lgm1F_ifr, lgm1F_yc,

      // FX vols
      fx_merge_sigma,

      // GRFN stuff
      void_prm, is_event,

      // Correlations
      corr_merge_r1_r2, corr_merge_r1_r3, corr_merge_r1_fx, corr_merge_r2_r3,
      corr_merge_r2_fx, corr_merge_r3_fx,

      grfn_payoff_4_3fquanto_tree,

      *num_prod, discount, *prod_val);

  /*	Add PV of Past */
  (*prod_val)[*num_prod - 1] += xStr.gd->pv_of_past;

FREE_RETURN:

  if (time)
    free(time);
  if (date)
    free(date);

  if (dom_name)
    free(dom_name);
  if (for_name)
    free(for_name);

  // Free raw term structures returned from the underlyings
  if (dom_sigma_dates)
    free(dom_sigma_dates);
  if (dom_sigma_values)
    free(dom_sigma_values);
  if (for_sigma_dates)
    free(for_sigma_dates);
  if (for_sigma_values)
    free(for_sigma_values);
  if (fx_sigma_dates)
    free(fx_sigma_dates);
  if (fx_sigma_values)
    free(fx_sigma_values);
  if (dom_sigma_merge_values)
    free(dom_sigma_merge_values);
  if (for_sigma_merge_values)
    free(for_sigma_merge_values);

  // Free the merged term structures
  if (vol_change)
    free(vol_change);
  if (r1_merge_sigma)
    free(r1_merge_sigma);
  if (r2_merge_sigma)
    free(r2_merge_sigma);
  if (r3_merge_sigma)
    free(r3_merge_sigma);
  if (fx_merge_sigma)
    free(fx_merge_sigma);
  if (corr_merge_r1_r2)
    free(corr_merge_r1_r2);
  if (corr_merge_r1_r3)
    free(corr_merge_r1_r3);
  if (corr_merge_r1_fx)
    free(corr_merge_r1_fx);
  if (corr_merge_r2_r3)
    free(corr_merge_r2_r3);
  if (corr_merge_r2_fx)
    free(corr_merge_r2_fx);
  if (corr_merge_r3_fx)
    free(corr_merge_r3_fx);

  // Free the ifr  , expectations and variance
  if (lgm2F_ifr)
    free(lgm2F_ifr);
  if (r1_fwd)
    free(r1_fwd);
  if (r1_var)
    free(r1_var);
  if (r2_fwd)
    free(r2_fwd);
  if (r2_var)
    free(r2_var);
  if (r1_r2_covar)
    free(r1_r2_covar);
  if (lgm1F_ifr)
    free(lgm1F_ifr);
  if (r3_fwd)
    free(r3_fwd);
  if (r3_var)
    free(r3_var);

  // Free the other times (barriers  , events)
  if (vol_times)
    free(vol_times);
  if (is_event)
    free(is_event);

  // Free the GRFN stuff
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
