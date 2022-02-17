/* ==========================================================================
   FILE_NAME:	Fx3FCalib.cxx

   PURPOSE:		Modelling of the spot FX vol by taking in consideration
   a LGM 1 factor on the domestic and foreign market and a lognormal model on
   the Spot Fx

   DATE:		05/25/00

   AUTHOR:		L.C.
   ========================================================================== */

#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

#define MAX_TIME 0.02
#define SHIFT_VOL 0.005
#define MAX_ERROR 0.0001
#define MAX_ERROR2 0.0005
#define SHIFT_VOL2 0.05
#define MAX_ERROR_VOL 0.0002

#define MAX_STP 3000
/*	Fill the time vector */
static Err fill_time_vector(double **time, int *nstp, int num_bar_times,
                            double *bar_times, int num_vol_times,
                            double *vol_times, int target_nstp) {
  Err err = NULL;

  /*	Add today if required */
  if ((*time)[0] < -EPS) {
    err = "Past event date in SrtGrfn3DFXTree";
    goto FREE_RETURN;
  }
  if ((*time)[0] > EPS) {
    num_f_add_number(nstp, time, 0.0);
    num_f_sort_vector(*nstp, *time);
    num_f_unique_vector(nstp, *time);
  }

  /*	If only one event today        , add empty event */
  if (*nstp == 1) {
    num_f_add_number(nstp, time, 1.0);
  }

  /*	Fill the vector */
  /*	New algorithm */
  num_f_fill_vector_newalgo(nstp, time, target_nstp);

  /*	Old algorithm */
  /*	num_f_fill_vector (*nstp        , time        , target_nstp        , 0);
          *nstp = (*nstp > target_nstp? *nstp: target_nstp);
          num_f_fill_vector_maxtime (nstp        , time        , 1.1 *
     (*time)[*nstp-1]/(*nstp-1)); */

FREE_RETURN:

  return err;
}

Err Fxbeta_log_approx(long today, double *maturity, long nb_mat,
                      double *sig_dom, double *sig_for, double *mat_fx,
                      long nb_mat_fx, double *sig_fx, double *beta,
                      double lam_dom, double lam_for, double corr_dom_for,
                      double corr_dom_fx, double corr_for_fx, double spot,
                      char *dom_yc, char *for_yc, double *time_fx,
                      long nb_time_fx, double *fx_fwd, double *fx_vol,
                      double max_time) {
  long i, j, k, k2;
  long index;
  double t, prev_t, dt, T;
  double logSpot;
  double *time2 = NULL;
  double *S = NULL;
  double volfx;
  double F, sigdom, sigfor, sigfx, x, bet;
  long nb_time;
  double last_mat;
  long last_index;

  Err err = NULL;

  logSpot = log(spot);

  /* discretise time */
  last_mat = time_fx[nb_time_fx - 1];
  nb_time = (long)(last_mat / max_time - 1.0E-08) + 1;
  if (nb_time < 2) {
    nb_time = 2;
  }

  time2 = calloc(nb_time, sizeof(double));
  if (!time2) {
    err = "Memory allocation failure (1) in fwd_approx";
    goto FREE_RETURN;
  }

  time2[0] = 0.0;

  for (i = 1; i < nb_time - 1; i++) {
    time2[i] = time2[i - 1] + max_time;
  }
  time2[i] = last_mat;

  num_f_concat_vector(&nb_time, &time2, nb_mat, maturity);
  num_f_concat_vector(&nb_time, &time2, nb_mat_fx, mat_fx);
  num_f_concat_vector(&nb_time, &time2, nb_time_fx, time_fx);
  num_f_sort_vector(nb_time, time2);
  num_f_unique_vector(&nb_time, time2);

  /* find the index of the last maturity */
  last_index = Get_Index(last_mat, time2, nb_time);

  S = (double *)calloc(last_index + 1, sizeof(double));

  if (!S) {
    err = "Memory allocation failure (2) in fwd_approx";
    goto FREE_RETURN;
  }

  /* diffuse */

  fx_fwd[0] = logSpot;
  S[0] = logSpot;

  for (i = 0; i < last_index; i++) {
    t = 0;
    T = time2[i + 1];
    /* initialisation at the fwd */
    F = logSpot + (swp_f_zr(today, today + T * 365.0, dom_yc) -
                   swp_f_zr(today, today + T * 365.0, for_yc)) *
                      T;

    k = 0;
    k2 = 0;
    sigdom = sig_dom[0];
    sigfor = sig_for[0];
    sigfx = sig_fx[0];
    bet = beta[0] - 1.0;

    for (j = 0; j < i + 1; j++) {
      prev_t = t;
      t = time2[j + 1];

      if (t > maturity[k]) {
        if (k < nb_mat - 1) {
          k += 1;
        }
        sigdom = sig_dom[k];
        sigfor = sig_for[k];
      }
      if (t > mat_fx[k2]) {
        if (k2 < nb_mat_fx - 1) {
          k2 += 1;
        }
        sigfx = sig_fx[k2];
        bet = beta[k2] - 1.0;
      }

      dt = t - prev_t;
      x = sigdom * (1.0 - exp(-lam_dom * (T - prev_t))) / lam_dom;

      F += x * dt *
           (corr_dom_fx * sigfx * exp(bet * S[j]) -
            corr_dom_for * sigfor * (1.0 - exp(-lam_for * (T - prev_t))) /
                lam_for +
            x);
    }
    S[i + 1] = F;
  }

  j = 0;
  /* fill the out structure */
  for (i = 0; i < nb_time_fx; i++) {
    t = time_fx[i];
    index = Get_Index(t, mat_fx, nb_mat_fx);
    volfx = sig_fx[index];
    bet = beta[index] - 1.0;
    fx_fwd[i] = S[Get_Index(t, time2, last_index + 1)];
    fx_vol[i] = exp(bet * fx_fwd[i]) * volfx;
  }

FREE_RETURN:

  if (time2) {
    free(time2);
  }

  if (S) {
    free(S);
  }

  return err;
}

/*	Calibration of a fx term structure to a set of fx options  */
Err Fx3DBetatsImpliedVol(long today, double opt_maturity, double start_date,
                         double end_date, double *maturity, long nbMat,
                         double *sig_curve_dom, double lda_dom,
                         double *sig_curve_for, double lda_for, double *mat_fx,
                         long nb_mat_fx, double *sig_curve_fx, double *beta,
                         double spot_fx, double corr_dom_for,
                         double corr_dom_fx, double corr_for_fx, char *dom_yc,
                         char *for_yc, double *fx_vol, double disc_dt,
                         double fx_dt)

{
  double *fx_vol_curve = NULL, *fx_fwd = NULL, *fx_time = NULL;

  long i;
  long nb_fx_time, nb_fx_time2;

  Err err = NULL;

  /* simple translation to make everything begin at 0 */

  nb_fx_time = (long)(end_date / fx_dt - 1.0E-08) + 1;
  if (nb_fx_time < 2) {
    nb_fx_time = 2;
  }

  fx_time = calloc(nb_fx_time, sizeof(double));

  if (!fx_time) {
    err = "Memory allocation error in Fx3DBetatsImpliedVol (1)";
    goto FREE_RETURN;
  }

  fx_time[0] = 0.0;

  /* fill the time where we want the volatility structure */
  for (i = 1; i < nb_fx_time - 1; i++) {
    fx_time[i] = fx_time[i - 1] + fx_dt;
  }

  fx_time[i] = end_date;

  num_f_concat_vector(&nb_fx_time, &fx_time, nbMat, maturity);
  num_f_concat_vector(&nb_fx_time, &fx_time, nb_mat_fx, mat_fx);
  num_f_sort_vector(nb_fx_time, fx_time);
  num_f_unique_vector(&nb_fx_time, fx_time);

  nb_fx_time2 = Get_Index(end_date, fx_time, nb_fx_time) + 1;

  fx_vol_curve = calloc(nb_fx_time2, sizeof(double));
  fx_fwd = calloc(nb_fx_time2, sizeof(double));

  if (!fx_vol || !fx_fwd) {
    err = "Memory allocation error in Fx3DBetatsImpliedVol (2)";
    goto FREE_RETURN;
  }

  err = Fxbeta_log_approx(today, maturity, nbMat, sig_curve_dom, sig_curve_for,
                          mat_fx, nb_mat_fx, sig_curve_fx, beta, lda_dom,
                          lda_for, corr_dom_for, corr_dom_fx, corr_for_fx,
                          spot_fx, dom_yc, for_yc, fx_time, nb_fx_time2, fx_fwd,
                          fx_vol_curve, disc_dt);
  if (err) {
    goto FREE_RETURN;
  }

  err = Fx3DtsImpliedVol(opt_maturity, start_date, end_date, maturity, nbMat,
                         sig_curve_dom, lda_dom, sig_curve_for, lda_for,
                         fx_time, fx_vol_curve, nb_fx_time2, corr_dom_for,
                         corr_dom_fx, corr_for_fx, fx_vol);

FREE_RETURN:

  if (fx_vol_curve)
    free(fx_vol_curve);
  if (fx_time)
    free(fx_time);
  if (fx_fwd)
    free(fx_fwd);

  return err;
}

/*	Implied vol direct from underlying */
Err Fx3DBetaImpliedVol(char *fx_underlying, double beta, double val_time,
                       double start_time, double end_time, double disc_dt,
                       double fx_dt, double *vol) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx,
      nb_merge_dates;
  long i;
  double *sigma_date_dom = NULL, *sigma_dom = NULL;
  double *tau_date_dom = NULL, *tau_dom = NULL;
  double *sigma_date_for = NULL, *sigma_for = NULL;
  double *tau_date_for = NULL, *tau_for = NULL;
  double *sigma_date_fx = NULL, *sigma_fx = NULL;
  double correl_dom_for, correl_dom_fx, correl_for_fx;
  double lda_dom, lda_for;
  double *sig_dom = NULL, *sig_for = NULL, *sig_fx = NULL, *merge_dates = NULL;

  SrtUndPtr fx_und, dom_und, for_und;
  char *domname, *forname;

  double *beta2 = NULL;

  long today;
  double spot_fx;
  char *dom_yc, *for_yc;

  Err err = NULL;

  err = Get_FX_StochRate_TermStructures(
      fx_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_dom_for, &correl_dom_fx, &correl_for_fx);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  /* now merge all the term structure */
  merge_dates = (double *)calloc(sigma_n_dom, sizeof(double));
  memcpy(merge_dates, sigma_date_dom, sigma_n_dom * sizeof(double));
  nb_merge_dates = sigma_n_dom;
  num_f_concat_vector(&nb_merge_dates, &merge_dates, sigma_n_for,
                      sigma_date_for);
  num_f_concat_vector(&nb_merge_dates, &merge_dates, sigma_n_fx, sigma_date_fx);
  num_f_sort_vector(nb_merge_dates, merge_dates);
  num_f_unique_vector(&nb_merge_dates, merge_dates);

  /*	Fill the new term structures */
  sig_dom = (double *)calloc(nb_merge_dates, sizeof(double));
  sig_for = (double *)calloc(nb_merge_dates, sizeof(double));
  sig_fx = (double *)calloc(nb_merge_dates, sizeof(double));
  beta2 = (double *)calloc(nb_merge_dates, sizeof(double));

  if (!sig_dom || !sig_for || !sig_fx || !beta2) {
    err = "Memory allocation error (1) in Fx3DBetaImpliedVol";
    goto FREE_RETURN;
  }

  for (i = nb_merge_dates - 1; i >= 0; i--) {
    sig_dom[i] =
        sigma_dom[Get_Index(merge_dates[i], sigma_date_dom, sigma_n_dom)];
    sig_for[i] =
        sigma_for[Get_Index(merge_dates[i], sigma_date_for, sigma_n_for)];
    sig_fx[i] = sigma_fx[Get_Index(merge_dates[i], sigma_date_fx, sigma_n_fx)];
    // beta2[i] = beta[Get_Index(merge_dates[i]        , sigma_date_fx        ,
    // sigma_n_fx)];
    beta2[i] = beta;
  }

  /* get today        , spot        , dom_yc        , for_yc */
  fx_und = lookup_und(fx_underlying);
  domname = get_domname_from_fxund(fx_und);
  forname = get_forname_from_fxund(fx_und);
  dom_und = lookup_und(domname);
  for_und = lookup_und(forname);

  err = get_underlying_discname(dom_und, &dom_yc);
  if (err) {
    goto FREE_RETURN;
  }
  err = get_underlying_discname(for_und, &for_yc);
  if (err) {
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(dom_und);

  err = Fx3DtsFwdFx(dom_und, for_und, fx_und, today, &spot_fx);

  err = Fx3DBetatsImpliedVol(
      today, val_time, start_time, end_time, merge_dates, nb_merge_dates,
      sig_dom, lda_dom, sig_for, lda_for, merge_dates, nb_merge_dates, sig_fx,
      beta2, spot_fx, correl_dom_for, correl_dom_fx, correl_for_fx, dom_yc,
      for_yc, vol, disc_dt, fx_dt);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sigma_date_fx) {
    free(sigma_date_fx);
  }

  if (sigma_fx) {
    free(sigma_fx);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (sig_fx) {
    free(sig_fx);
  }

  if (beta2) {
    free(beta2);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  return err;
}

Err Fx3DBetatsCalibration(long today, double *exercise_opt,
                          double *maturity_opt, double *vol_opt, long nbrOpt,
                          double *maturity, long nbrMat, double *sig_curve_dom,
                          double lda_dom, double *sig_curve_for, double lda_for,
                          double *beta, double spot_fx, double correl_dom_for,
                          double correl_dom_fx, double correl_for_fx,
                          char *dom_yc, char *for_yc, double **fx_vol_curve,
                          double disc_dt, double fx_dt, long nbIterMax)

{
  long i, nIter;
  double T, mat, Fwd;
  double bet;
  double *fx_vol_log = NULL;

  double vol_imp1, vol_imp2;
  double vol_tgt;
  double error;

  Err err = NULL;

  (*fx_vol_curve) = NULL;
  (*fx_vol_curve) = (double *)calloc(nbrOpt, sizeof(double));

  if (!(*fx_vol_curve)) {
    err = "Memory allocation error in Fx3DBetatsCalibration";
    goto FREE_RETURN;
  }

  /* First we calibrate a lognormal model */
  err = Fx3DtsCalibration(exercise_opt, maturity_opt, vol_opt, nbrOpt, maturity,
                          nbrMat, sig_curve_dom, lda_dom, sig_curve_for,
                          lda_for, correl_dom_for, correl_dom_fx, correl_for_fx,
                          &fx_vol_log);
  if (err) {
    goto FREE_RETURN;
  }

  for (i = 0; i < nbrOpt; i++) {
    /* Initialisation: SigBeta * F^beta = SigLn * F */
    mat = maturity_opt[i];
    T = exercise_opt[i];
    Fwd = spot_fx * swp_f_df(today, today + T * 365.0, for_yc) /
          swp_f_df(today, today + T * 365.0, dom_yc);
    bet = beta[i] - 1.0;
    vol_tgt = vol_opt[i];
    (*fx_vol_curve)[i] = fx_vol_log[i] * exp(-bet * log(Fwd));

    err = Fx3DBetatsImpliedVol(
        today, mat, 0.0, T, maturity, nbrMat, sig_curve_dom, lda_dom,
        sig_curve_for, lda_for, exercise_opt, nbrOpt, (*fx_vol_curve), beta,
        spot_fx, correl_dom_for, correl_dom_fx, correl_for_fx, dom_yc, for_yc,
        &vol_imp1, disc_dt, fx_dt);
    if (err) {
      goto FREE_RETURN;
    }

    error = fabs(vol_imp1 - vol_tgt);

    nIter = 0;

    while (nIter <= nbIterMax && error > MAX_ERROR) {
      (*fx_vol_curve)[i] += SHIFT_VOL;
      err = Fx3DBetatsImpliedVol(
          today, mat, 0.0, T, maturity, nbrMat, sig_curve_dom, lda_dom,
          sig_curve_for, lda_for, exercise_opt, nbrOpt, (*fx_vol_curve), beta,
          spot_fx, correl_dom_for, correl_dom_fx, correl_for_fx, dom_yc, for_yc,
          &vol_imp2, disc_dt, fx_dt);

      if (err) {
        goto FREE_RETURN;
      }

      (*fx_vol_curve)[i] -= SHIFT_VOL;
      (*fx_vol_curve)[i] +=
          (vol_tgt - vol_imp1) * SHIFT_VOL / (vol_imp2 - vol_imp1);

      if ((*fx_vol_curve)[i] < 0) {
        err = serror("Cannot calibrate option %d", i);
        goto FREE_RETURN;
      }

      err = Fx3DBetatsImpliedVol(
          today, mat, 0.0, T, maturity, nbrMat, sig_curve_dom, lda_dom,
          sig_curve_for, lda_for, exercise_opt, nbrOpt, (*fx_vol_curve), beta,
          spot_fx, correl_dom_for, correl_dom_fx, correl_for_fx, dom_yc, for_yc,
          &vol_imp1, disc_dt, fx_dt);

      if (err) {
        goto FREE_RETURN;
      }

      error = fabs(vol_imp1 - vol_tgt);

      nIter += 1;
    }
  }

FREE_RETURN:

  if (fx_vol_log)
    free(fx_vol_log);

  return err;
}

Err Fx3DBetaCalibration(char *dom_underlying, char *for_underlying,
                        double spot_fx, double beta, double correl_dom_for,
                        double correl_dom_fx, double correl_for_fx,
                        double *exercise_opt, double *maturity_opt,
                        double *vol_opt, long nbropt, double **fx_vol_curve,
                        double disc_dt, double fx_dt, long nbIterMax)

{
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL;

  double lda_dom, lda_for;

  SrtUndPtr dom_und, for_und;

  double *beta2 = NULL;

  long today, i;
  char *dom_yc, *for_yc;

  Err err = NULL;

  err =
      Get_LGM_TermStructure(dom_underlying, &sigma_date_dom, &sigma_dom,
                            &sigma_n_dom, &tau_date_dom, &tau_dom, &tau_n_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err =
      Get_LGM_TermStructure(for_underlying, &sigma_date_for, &sigma_for,
                            &sigma_n_for, &tau_date_for, &tau_for, &tau_n_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  dom_und = lookup_und(dom_underlying);
  for_und = lookup_und(for_underlying);

  err = get_underlying_discname(dom_und, &dom_yc);
  if (err) {
    goto FREE_RETURN;
  }
  err = get_underlying_discname(for_und, &for_yc);
  if (err) {
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(dom_und);

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  if (err) {
    goto FREE_RETURN;
  }

  beta2 = (double *)calloc(nbropt, sizeof(double));
  for (i = 0; i < nbropt; i++) {
    beta2[i] = beta;
  }

  err = Fx3DBetatsCalibration(
      today, exercise_opt, maturity_opt, vol_opt, nbropt, merge_dates,
      nb_merge_dates, sigma_dom, lda_dom, sigma_for, lda_for, beta2, spot_fx,
      correl_dom_for, correl_dom_fx, correl_for_fx, dom_yc, for_yc,
      fx_vol_curve, disc_dt, fx_dt, nbIterMax);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  if (beta2) {
    free(beta2);
  }

  return err;
}

Err FxCall_payoff_4_3dfxBeta_tree(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double spot_fx, void *dom_yc, double dom_lam, double dom_phi, void *for_yc,
    double for_lam, double for_phi,
    /* Nodes data */
    long n1, long n2, long n3,
    /* i: d1        , j: d2        , k: d3        , l = {0: xDom        , 1:
       xFor , 2: log (Fx/Fx0)} */
    double beta, double ****sv,
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
        Fx = exp((log((1.0 - beta) * sv[i][j][k][2] + 1.0)) / (1.0 - beta));

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

Err Fx3DBetatsTreeFxOptions(
    long today, long maturity_date, double *strikes, long nbrOpt,
    double *maturity, long nbrMat, double *sig_curve_dom, double dom_lam,
    double *sig_curve_for, double for_lam, double *maturity_fx, long nbrMat_fx,
    double *sig_curve_fx, double alpha, double beta, double spot_fx,
    double corr_dom_for, double corr_dom_fx, double corr_for_fx, char *dom_yc,
    char *for_yc, double *option_prices, long num_stp) {
  double maturity_opt;
  double *time = NULL, *date = NULL;
  long nstp = 1;
  long i, j;

  int *vol_change = NULL;
  double *dom_vol = NULL, *for_vol = NULL, *fx_vol = NULL,
         *corr_dom_for_ts = NULL, *corr_dom_fx_ts = NULL,
         *corr_for_fx_ts = NULL;

  double *dom_ifr = NULL, *dom_fwd = NULL, *dom_var = NULL, *for_ifr = NULL,
         *for_fwd = NULL, *for_var = NULL, *fx_fwd = NULL, *fx_var = NULL;

  double *sig_fx_approx = NULL, *fx_fwd_approx = NULL, *beta_tab = NULL;

  double **func_parm = NULL;
  int *has_evt = NULL;

  int *is_bar = NULL;
  double *bar_lvl = NULL;
  int *bar_cl = NULL;

  double beta2 = 1.0 - beta;

  Err err = NULL;

  smessage("Phase 1 -preprocessing        , time in sec: 0");

  maturity_opt = (maturity_date - today) / DAYS_IN_YEAR;

  time = (double *)calloc(nstp, sizeof(double));
  if (!time) {
    err = "Memory allocation error (1) in Fx3DBetatsTreeFxOptions";
    goto FREE_RETURN;
  }

  time[0] = maturity_opt;

  /*	Fill the time vector */
  err = fill_time_vector(&time, &nstp, 0, NULL, 0, NULL, num_stp);

  if (err) {
    goto FREE_RETURN;
  }

  date = (double *)calloc(nstp, sizeof(double));
  has_evt = (int *)calloc(nstp, sizeof(int));
  func_parm = dmatrix(0, nstp - 1, 0, nbrOpt - 1);

  if (!date || !func_parm || !has_evt) {
    err = "Memory allocation error (2) in Fx3DBetatsTreeFxOptions";
    goto FREE_RETURN;
  }

  for (j = 0; j < nbrOpt; j++) {
    func_parm[nstp - 1][j] = strikes[j];
  }

  for (i = 0; i < nstp; i++) {
    date[i] = today + DAYS_IN_YEAR * time[i];
    has_evt[i] = 0;

    if (i > 0 && date[i] - date[i - 1] >= 1) {
      date[i] = (long)(date[i] + 1.0e-08);
      time[i] = YEARS_IN_DAY * (date[i] - today);
    }
  }
  has_evt[nstp - 1] = 1;

  /*	Fill the model parameters as required by tree_main_3dfx */

  vol_change = (int *)calloc(nstp, sizeof(int));
  dom_vol = (double *)calloc(nstp, sizeof(double));
  for_vol = (double *)calloc(nstp, sizeof(double));
  fx_vol = (double *)calloc(nstp, sizeof(double));
  corr_dom_for_ts = (double *)calloc(nstp, sizeof(double));
  corr_dom_fx_ts = (double *)calloc(nstp, sizeof(double));
  corr_for_fx_ts = (double *)calloc(nstp, sizeof(double));

  beta_tab = (double *)calloc(nstp, sizeof(double));

  is_bar = (int *)calloc(nstp, sizeof(int));
  bar_lvl = (double *)calloc(nstp, sizeof(double));
  bar_cl = (int *)calloc(nstp, sizeof(int));

  if (!vol_change || !dom_vol || !for_vol || !fx_vol || !beta_tab || !is_bar ||
      !bar_lvl || !bar_cl) {
    err = "Memory allocation error (3) in Fx3DBetatsTreeFxOptions";
    goto FREE_RETURN;
  }

  vol_change[nstp - 1] = 1;

  dom_vol[nstp - 1] =
      sig_curve_dom[Get_Index(time[nstp - 1], maturity, nbrMat)];
  for_vol[nstp - 1] =
      sig_curve_for[Get_Index(time[nstp - 1], maturity, nbrMat)];
  fx_vol[nstp - 1] =
      sig_curve_fx[Get_Index(time[nstp - 1], maturity_fx, nbrMat_fx)] *
      exp(alpha * sqrt(time[nstp - 1]) - 0.5 * alpha * alpha * time[nstp - 1]);
  ;
  beta_tab[nstp - 1] = beta;

  corr_dom_for_ts[nstp - 1] = corr_dom_for;
  corr_dom_fx_ts[nstp - 1] = corr_dom_fx;
  corr_for_fx_ts[nstp - 1] = corr_for_fx;

  for (i = nstp - 2; i >= 0; i--) {
    dom_vol[i] = sig_curve_dom[Get_Index(time[i], maturity, nbrMat)];
    for_vol[i] = sig_curve_for[Get_Index(time[i], maturity, nbrMat)];
    fx_vol[i] = sig_curve_fx[Get_Index(time[i], maturity_fx, nbrMat_fx)] *
                exp(alpha * sqrt(time[i]) - 0.5 * alpha * alpha * time[i]);

    corr_dom_for_ts[i] = corr_dom_for;
    corr_dom_fx_ts[i] = corr_dom_fx;
    corr_for_fx_ts[i] = corr_for_fx;

    beta_tab[i] = beta;

    if (fabs(dom_vol[i] - dom_vol[i + 1]) + fabs(for_vol[i] - for_vol[i + 1]) +
            fabs(fx_vol[i] - fx_vol[i + 1]) >
        EPS) {
      vol_change[i] = 1;
    } else {
      vol_change[i] = 0;
    }
  }

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
      !corr_dom_for_ts || !corr_dom_fx_ts || !corr_for_fx_ts || !fx_fwd ||
      !fx_var || !sig_fx_approx || !fx_fwd_approx) {
    err = "Memory allocation error (3) in Fx3DBetatsTreeFxOptions";
    goto FREE_RETURN;
  }

  /* first get the coresponding lognormal volatilities */
  err = Fxbeta_log_approx(today, time, nstp, dom_vol, for_vol, time, nstp,
                          fx_vol, beta_tab, dom_lam, for_lam, corr_dom_for,
                          corr_dom_fx, corr_for_fx, spot_fx, dom_yc, for_yc,
                          time, nstp, fx_fwd_approx, sig_fx_approx, MAX_TIME);

  fill_fwd_var(nstp, time, date, dom_vol, for_vol, sig_fx_approx, dom_lam,
               for_lam, corr_dom_for, corr_dom_fx, corr_for_fx, dom_yc, for_yc,
               dom_ifr, dom_fwd, dom_var, for_ifr, for_fwd, for_var, fx_fwd,
               fx_var);

  for (i = nstp - 1; i >= 0; i--) {
    bar_lvl[i] = (exp(beta2 * (fx_fwd_approx[i])) - 1.0) / beta2;
  }

  err = treeBeta_main_3dfx(
      nstp, time, date, vol_change, dom_vol, for_vol, fx_vol, beta_tab, dom_ifr,
      dom_fwd, dom_var, for_ifr, for_fwd, for_var, fx_fwd, fx_var, func_parm,
      has_evt, bar_lvl, bar_cl, is_bar, dom_lam, for_lam, corr_dom_for_ts,
      corr_dom_fx_ts, corr_for_fx_ts, spot_fx, dom_yc, for_yc,
      FxCall_payoff_4_3dfxBeta_tree, nbrOpt, 1, option_prices);

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
  if (corr_dom_for_ts)
    free(corr_dom_for_ts);
  if (corr_dom_fx_ts)
    free(corr_dom_fx_ts);
  if (corr_for_fx_ts)
    free(corr_for_fx_ts);

  if (fx_fwd_approx)
    free(fx_fwd_approx);
  if (sig_fx_approx)
    free(sig_fx_approx);
  if (beta_tab)
    free(beta_tab);

  if (bar_lvl)
    free(bar_lvl);
  if (is_bar)
    free(is_bar);
  if (bar_cl)
    free(bar_cl);

  if (func_parm) {
    free_dmatrix(func_parm, 0, nstp - 1, 0, nbrOpt - 1);
  }

  return err;
}

Err Fx3DAlphaBetatsTreeFxOptions(
    long today, long maturity_date, double *strikes, long nbrOpt,
    double *maturity, long nbrMat, double *sig_curve_dom, double dom_lam,
    double *sig_curve_for, double for_lam, double *maturity_fx, long nbrMat_fx,
    double *sig_curve_fx, double alpha, double beta, double spot_fx,
    double corr_dom_for, double corr_dom_fx, double corr_for_fx, char *dom_yc,
    char *for_yc, double *option_prices, long num_stp) {
  char *cerr = NULL;
  double *pv1 = NULL, *pv2 = NULL;
  int i;

  pv1 = (double *)calloc(nbrOpt + 1, sizeof(double));
  if (!pv1) {
    cerr = "Memory allocation error in Fx3DAlphaBetatsTreeFxOptions";
    goto FREE_RETURN;
  }

  pv2 = (double *)calloc(nbrOpt + 1, sizeof(double));
  if (!pv2) {
    cerr = "Memory allocation error in Fx3DAlphaBetatsTreeFxOptions";
    goto FREE_RETURN;
  }

  if (fabs(alpha) <= 1.0e-03) {
    alpha = 0.0;
  }

  cerr = Fx3DBetatsTreeFxOptions(
      today, maturity_date, strikes, nbrOpt, maturity, nbrMat, sig_curve_dom,
      dom_lam, sig_curve_for, for_lam, maturity_fx, nbrMat_fx, sig_curve_fx,
      alpha, beta, spot_fx, corr_dom_for, corr_dom_fx, corr_for_fx, dom_yc,
      for_yc, pv1, num_stp);

  if (cerr) {
    goto FREE_RETURN;
  }

  if (fabs(alpha) > 1.0e-03) {
    cerr = Fx3DBetatsTreeFxOptions(
        today, maturity_date, strikes, nbrOpt, maturity, nbrMat, sig_curve_dom,
        dom_lam, sig_curve_for, for_lam, maturity_fx, nbrMat_fx, sig_curve_fx,
        -alpha, beta, spot_fx, corr_dom_for, corr_dom_fx, corr_for_fx, dom_yc,
        for_yc, pv2, num_stp);

    if (cerr) {
      goto FREE_RETURN;
    }
  } else {
    memcpy(pv2, pv1, nbrOpt * sizeof(double));
  }

  for (i = 0; i < nbrOpt; i++) {
    option_prices[i] = 0.5 * (pv1[i] + pv2[i]);
  }

FREE_RETURN:

  if (pv1) {
    free(pv1);
  }

  if (pv2) {
    free(pv2);
  }

  return cerr;
}

Err Fx3DBetatsCalibration2(long today, double *exercise_opt,
                           double *maturity_opt, double *vol_opt, long nbrOpt,
                           long nbrLong, double *maturity, long nbrMat,
                           double *sig_curve_dom, double lda_dom,
                           double *sig_curve_for, double lda_for, double beta,
                           double spot_fx, double corr_dom_for,
                           double corr_dom_fx, double corr_for_fx, char *dom_yc,
                           char *for_yc, double **fx_vol_curve, long nbSteps,
                           long nbNewton, double disc_dt, double fx_dt,
                           long nbIterMax)

{
  long i, j;
  long nbrShort;
  long maturity_date, exercise_date, spot_date;
  double Fwd;
  double strikes[1];
  double df;
  double price1[2], error;
  double shift, vol_imp;
  double *beta2 = NULL;
  double delta;

  Err err = NULL;

  if (fabs(beta - 1.0) < 1.0E-04) {
    return Fx3DtsCalibration(exercise_opt, maturity_opt, vol_opt, nbrOpt,
                             maturity, nbrMat, sig_curve_dom, lda_dom,
                             sig_curve_for, lda_for, corr_dom_for, corr_dom_fx,
                             corr_for_fx, fx_vol_curve);
  }

  nbrShort = nbrOpt - nbrLong;

  spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  /* calculation of the spot spot fx */
  spot_fx *=
      swp_f_df(today, spot_date, dom_yc) / swp_f_df(today, spot_date, for_yc);

  beta2 = (double *)calloc(nbrOpt, sizeof(double));
  if (!beta2) {
    goto FREE_RETURN;
  }

  for (i = 0; i < nbrOpt; i++) {
    beta2[i] = beta;
  }

  err = Fx3DBetatsCalibration(today, exercise_opt, maturity_opt, vol_opt,
                              nbrOpt, maturity, nbrMat, sig_curve_dom, lda_dom,
                              sig_curve_for, lda_for, beta2, spot_fx,
                              corr_dom_for, corr_dom_fx, corr_for_fx, dom_yc,
                              for_yc, fx_vol_curve, disc_dt, fx_dt, nbIterMax);
  if (err) {
    goto FREE_RETURN;
  }

  for (i = nbrShort; i < nbrOpt; i++) {
    smessage("Calibration of the option %d\n", i + 1);

    maturity_date = (long)(today + maturity_opt[i] * 365.0000000001);
    exercise_date = (long)(today + exercise_opt[i] * 365.0000000001);
    df = swp_f_df(today, exercise_date, dom_yc);
    Fwd = spot_fx * swp_f_df(today, exercise_date, for_yc) /
          swp_f_df(today, exercise_date, dom_yc);

    /* prices in the tree */
    strikes[0] = Fwd;

    error = 1.0E10;
    j = 0;

    while (j < nbNewton && (error > MAX_ERROR_VOL)) {
      err = Fx3DBetatsTreeFxOptions(
          today, exercise_date, strikes, 1, maturity, nbrMat, sig_curve_dom,
          lda_dom, sig_curve_for, lda_for, exercise_opt, nbrOpt, *fx_vol_curve,
          0.0, beta, spot_fx, corr_dom_for, corr_dom_fx, corr_for_fx, dom_yc,
          for_yc, &(price1[0]), nbSteps);
      if (err) {
        goto FREE_RETURN;
      }

      /* calculates the model implied vol */
      err = srt_f_optimpvol(price1[0], Fwd, Fwd, exercise_opt[i], df, 0, 0,
                            &vol_imp);

      if (err) {
        shift = 0.0;
        nbSteps = (long)(nbSteps * 1.2);
      } else {
        error = fabs(vol_imp - vol_opt[i]);

        if (error > MAX_ERROR_VOL) {
          if (i > 0) {
            delta = (*fx_vol_curve)[i] * (*fx_vol_curve)[i] *
                        exp(2 * (beta2[i] - 1) * log(Fwd)) -
                    (vol_imp * vol_imp - vol_opt[i] * vol_opt[i]) *
                        exercise_opt[i] /
                        (exercise_opt[i] - exercise_opt[i - 1]);
          } else {
            delta = (*fx_vol_curve)[i] * (*fx_vol_curve)[i] *
                        exp(2 * (beta2[i] - 1) * log(Fwd)) -
                    (vol_imp * vol_imp - vol_opt[i] * vol_opt[i]);
          }

          if (delta > 0) {
            shift = -(*fx_vol_curve)[i] +
                    sqrt(delta) * exp((1.0 - beta2[i]) * log(Fwd));
            /*
            for (k=i; k<nbrOpt; k++)
            {
                    (*fx_vol_curve)[k] += shift ;
            }
            */
            (*fx_vol_curve)[i] += shift;
          } else {
            shift = 0.0;
            nbSteps = (long)(nbSteps * 1.2);
          }
        }
      }
      j += 1;
    }
  }

FREE_RETURN:

  if (beta2)
    free(beta2);
  return err;
}

Err Fx3DBetaCalibration2(char *dom_underlying, char *for_underlying,
                         double spot_fx, double beta, double correl_dom_for,
                         double correl_dom_fx, double correl_for_fx,
                         double *exercise_opt, double *maturity_opt,
                         double *vol_opt, long nbropt, long nbLong,
                         double **fx_vol_curve, long nbSteps, long nbNewton,
                         double disc_dt, double fx_dt, long nbIterMax)

{
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL;

  double lda_dom, lda_for;

  SrtUndPtr dom_und, for_und;

  long today;
  char *dom_yc, *for_yc;

  Err err = NULL;

  err =
      Get_LGM_TermStructure(dom_underlying, &sigma_date_dom, &sigma_dom,
                            &sigma_n_dom, &tau_date_dom, &tau_dom, &tau_n_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err =
      Get_LGM_TermStructure(for_underlying, &sigma_date_for, &sigma_for,
                            &sigma_n_for, &tau_date_for, &tau_for, &tau_n_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  dom_und = lookup_und(dom_underlying);
  for_und = lookup_und(for_underlying);

  err = get_underlying_discname(dom_und, &dom_yc);
  if (err) {
    goto FREE_RETURN;
  }
  err = get_underlying_discname(for_und, &for_yc);
  if (err) {
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(dom_und);

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  if (err) {
    goto FREE_RETURN;
  }

  err = Fx3DBetatsCalibration2(
      today, exercise_opt, maturity_opt, vol_opt, nbropt, nbLong, merge_dates,
      nb_merge_dates, sigma_dom, lda_dom, sigma_for, lda_for, beta, spot_fx,
      correl_dom_for, correl_dom_fx, correl_for_fx, dom_yc, for_yc,
      fx_vol_curve, nbSteps, nbNewton, disc_dt, fx_dt, nbIterMax);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  return err;
}

Err Fx3DAlphaBetatsCalibration2(
    long today, double *exercise_opt, double *maturity_opt, double *vol_opt,
    long nbrOpt, long nbrLong, double *maturity, long nbrMat,
    double *sig_curve_dom, double lda_dom, double *sig_curve_for,
    double lda_for, double alpha, double beta, double spot_fx,
    double corr_dom_for, double corr_dom_fx, double corr_for_fx, char *dom_yc,
    char *for_yc, double **fx_vol_curve, long nbSteps, long nbNewton,
    double disc_dt, double fx_dt, long nbIterMax)

{
  long i, j;
  long nbrShort;
  long maturity_date, exercise_date, spot_date;
  double Fwd;
  double strikes[1];
  double df;
  double price1[2], error;
  double shift, vol_imp;
  double *beta2 = NULL;
  double delta;

  Err err = NULL;

  if (alpha < 1.0E-04) {
    return Fx3DBetatsCalibration2(
        today, exercise_opt, maturity_opt, vol_opt, nbrOpt, nbrLong, maturity,
        nbrMat, sig_curve_dom, lda_dom, sig_curve_for, lda_for, beta, spot_fx,
        corr_dom_for, corr_dom_fx, corr_for_fx, dom_yc, for_yc, fx_vol_curve,
        nbSteps, nbNewton, disc_dt, fx_dt, nbIterMax);
  }

  nbrShort = nbrOpt - nbrLong;

  spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  /* calculation of the spot spot fx */
  spot_fx *=
      swp_f_df(today, spot_date, dom_yc) / swp_f_df(today, spot_date, for_yc);

  beta2 = (double *)calloc(nbrOpt, sizeof(double));
  if (!beta2) {
    goto FREE_RETURN;
  }

  for (i = 0; i < nbrOpt; i++) {
    beta2[i] = beta;
  }

  err = Fx3DBetatsCalibration(today, exercise_opt, maturity_opt, vol_opt,
                              nbrOpt, maturity, nbrMat, sig_curve_dom, lda_dom,
                              sig_curve_for, lda_for, beta2, spot_fx,
                              corr_dom_for, corr_dom_fx, corr_for_fx, dom_yc,
                              for_yc, fx_vol_curve, disc_dt, fx_dt, nbIterMax);
  if (err) {
    goto FREE_RETURN;
  }

  for (i = nbrShort; i < nbrOpt; i++) {
    smessage("Calibration of the option %d\n", i + 1);

    maturity_date = (long)(today + maturity_opt[i] * 365.0000000001);
    exercise_date = (long)(today + exercise_opt[i] * 365.0000000001);
    df = swp_f_df(today, exercise_date, dom_yc);
    Fwd = spot_fx * swp_f_df(today, exercise_date, for_yc) /
          swp_f_df(today, exercise_date, dom_yc);

    /* prices in the tree */
    strikes[0] = Fwd;

    error = 1.0E10;
    j = 0;

    while (j < nbNewton && (error > MAX_ERROR_VOL)) {
      err = Fx3DAlphaBetatsTreeFxOptions(
          today, exercise_date, strikes, 1, maturity, nbrMat, sig_curve_dom,
          lda_dom, sig_curve_for, lda_for, exercise_opt, nbrOpt, *fx_vol_curve,
          alpha, beta, spot_fx, corr_dom_for, corr_dom_fx, corr_for_fx, dom_yc,
          for_yc, &(price1[0]), nbSteps);
      if (err) {
        goto FREE_RETURN;
      }

      /* calculates the model implied vol */
      err = srt_f_optimpvol(price1[0], Fwd, Fwd, exercise_opt[i], df, 0, 0,
                            &vol_imp);

      if (err) {
        shift = 0.0;
        nbSteps = (long)(nbSteps * 1.2);
      } else {
        error = fabs(vol_imp - vol_opt[i]);

        if (error > MAX_ERROR_VOL) {
          if (i > 0) {
            delta = (*fx_vol_curve)[i] * (*fx_vol_curve)[i] *
                        exp(2 * (beta2[i] - 1) * log(Fwd)) -
                    (vol_imp * vol_imp - vol_opt[i] * vol_opt[i]) *
                        exercise_opt[i] /
                        (exercise_opt[i] - exercise_opt[i - 1]);
          } else {
            delta = (*fx_vol_curve)[i] * (*fx_vol_curve)[i] *
                        exp(2 * (beta2[i] - 1) * log(Fwd)) -
                    (vol_imp * vol_imp - vol_opt[i] * vol_opt[i]);
          }

          if (delta > 0) {
            shift = -(*fx_vol_curve)[i] +
                    sqrt(delta) * exp((1.0 - beta2[i]) * log(Fwd));
            /*
            for (k=i; k<nbrOpt; k++)
            {
                    (*fx_vol_curve)[k] += shift ;
            }
            */
            (*fx_vol_curve)[i] += shift;
          } else {
            shift = 0.0;
            nbSteps = (long)(nbSteps * 1.2);
          }
        }
      }
      j += 1;
    }
  }

FREE_RETURN:

  if (beta2)
    free(beta2);
  return err;
}

Err Fx3DAlphaBetaCalibration2(char *dom_underlying, char *for_underlying,
                              double spot_fx, double alpha, double beta,
                              double correl_dom_for, double correl_dom_fx,
                              double correl_for_fx, double *exercise_opt,
                              double *maturity_opt, double *vol_opt,
                              long nbropt, long nbLong, double **fx_vol_curve,
                              long nbSteps, long nbNewton, double disc_dt,
                              double fx_dt, long nbIterMax)

{
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL;

  double lda_dom, lda_for;

  SrtUndPtr dom_und, for_und;

  long today;
  char *dom_yc, *for_yc;

  Err err = NULL;

  err =
      Get_LGM_TermStructure(dom_underlying, &sigma_date_dom, &sigma_dom,
                            &sigma_n_dom, &tau_date_dom, &tau_dom, &tau_n_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err =
      Get_LGM_TermStructure(for_underlying, &sigma_date_for, &sigma_for,
                            &sigma_n_for, &tau_date_for, &tau_for, &tau_n_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  dom_und = lookup_und(dom_underlying);
  for_und = lookup_und(for_underlying);

  err = get_underlying_discname(dom_und, &dom_yc);
  if (err) {
    goto FREE_RETURN;
  }
  err = get_underlying_discname(for_und, &for_yc);
  if (err) {
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(dom_und);

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  if (err) {
    goto FREE_RETURN;
  }

  err = Fx3DAlphaBetatsCalibration2(
      today, exercise_opt, maturity_opt, vol_opt, nbropt, nbLong, merge_dates,
      nb_merge_dates, sigma_dom, lda_dom, sigma_for, lda_for, alpha, beta,
      spot_fx, correl_dom_for, correl_dom_fx, correl_for_fx, dom_yc, for_yc,
      fx_vol_curve, nbSteps, nbNewton, disc_dt, fx_dt, nbIterMax);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  return err;
}

Err Fxbeta_log_approx_corr(long today, double *maturity, long nb_mat,
                           double *sig_dom, double *sig_for, double *mat_fx,
                           long nb_mat_fx, double *sig_fx, double *beta,
                           double lam_dom, double lam_for, double *mat_corr,
                           double *corr_dom_for_ts, double *corr_dom_fx_ts,
                           double *corr_for_fx_ts, long nb_mat_corr,
                           double spot, char *dom_yc, char *for_yc,
                           double *time_fx, long nb_time_fx, double *fx_fwd,
                           double *fx_vol, double max_time) {
  long i, j, k, k2, k3;
  long index;
  double t, prev_t, dt, T;
  double logSpot;
  double *time2 = NULL;
  double *S = NULL;
  double volfx;
  double F, sigdom, sigfor, sigfx, x, bet;
  double corr_dom_for, corr_dom_fx, corr_for_fx;
  long nb_time;
  double last_mat;
  long last_index;

  Err err = NULL;

  logSpot = log(spot);

  /* discretise time */
  last_mat = time_fx[nb_time_fx - 1];
  nb_time = (long)(last_mat / max_time - 1.0E-08) + 1;
  if (nb_time < 2) {
    nb_time = 2;
  }

  time2 = calloc(nb_time, sizeof(double));
  if (!time2) {
    err = "Memory allocation failure (1) in fwd_approx";
    goto FREE_RETURN;
  }

  time2[0] = 0.0;

  for (i = 1; i < nb_time - 1; i++) {
    time2[i] = time2[i - 1] + max_time;
  }
  time2[i] = last_mat;

  num_f_concat_vector(&nb_time, &time2, nb_mat, maturity);
  num_f_concat_vector(&nb_time, &time2, nb_mat_fx, mat_fx);
  num_f_concat_vector(&nb_time, &time2, nb_time_fx, time_fx);
  num_f_sort_vector(nb_time, time2);
  num_f_unique_vector(&nb_time, time2);

  /* find the index of the last maturity */
  last_index = Get_Index(last_mat, time2, nb_time);

  S = (double *)calloc(last_index + 1, sizeof(double));

  if (!S) {
    err = "Memory allocation failure (2) in fwd_approx";
    goto FREE_RETURN;
  }

  /* diffuse */

  fx_fwd[0] = logSpot;
  S[0] = logSpot;

  for (i = 0; i < last_index; i++) {
    t = 0;
    T = time2[i + 1];
    /* initialisation at the fwd */
    F = logSpot + (swp_f_zr(today, today + T * 365.0, dom_yc) -
                   swp_f_zr(today, today + T * 365.0, for_yc)) *
                      T;

    k = 0;
    k2 = 0;
    k3 = 0;

    sigdom = sig_dom[0];
    sigfor = sig_for[0];
    sigfx = sig_fx[0];
    corr_dom_for = corr_dom_for_ts[0];
    corr_dom_fx = corr_dom_fx_ts[0];
    corr_for_fx = corr_for_fx_ts[0];

    bet = beta[0] - 1.0;

    for (j = 0; j < i + 1; j++) {
      prev_t = t;
      t = time2[j + 1];

      if (t > maturity[k]) {
        if (k < nb_mat - 1) {
          k += 1;
          sigdom = sig_dom[k];
          sigfor = sig_for[k];
        }
      }
      if (t > mat_fx[k2]) {
        if (k2 < nb_mat_fx - 1) {
          k2 += 1;
          sigfx = sig_fx[k2];
          bet = beta[k2] - 1.0;
        }
      }
      if (t > mat_corr[k3]) {
        if (k3 < nb_mat_corr - 1) {
          k3 += 1;
          corr_dom_for = corr_dom_for_ts[k3];
          corr_dom_fx = corr_dom_fx_ts[k3];
          corr_for_fx = corr_for_fx_ts[k3];
        }
      }

      dt = t - prev_t;
      x = sigdom * (1.0 - exp(-lam_dom * (T - prev_t))) / lam_dom;

      F += x * dt *
           (corr_dom_fx * sigfx * exp(bet * S[j]) -
            corr_dom_for * sigfor * (1.0 - exp(-lam_for * (T - prev_t))) /
                lam_for +
            x);
    }
    S[i + 1] = F;
  }

  j = 0;
  /* fill the out structure */
  for (i = 0; i < nb_time_fx; i++) {
    t = time_fx[i];
    index = Get_Index(t, mat_fx, nb_mat_fx);
    volfx = sig_fx[index];
    bet = beta[index] - 1.0;
    fx_fwd[i] = S[Get_Index(t, time2, last_index + 1)];
    fx_vol[i] = exp(bet * fx_fwd[i]) * volfx;
  }

FREE_RETURN:

  if (time2) {
    free(time2);
  }

  if (S) {
    free(S);
  }

  return err;
}

/*	Calibration of a fx term structure to a set of fx options Corr TS
 * available */
Err Fx3DBetatsImpliedVol_corr(
    long today, double opt_maturity, double start_date, double end_date,
    double *maturity, long nbMat, double *sig_curve_dom, double lda_dom,
    double *sig_curve_for, double lda_for, double *mat_fx, long nb_mat_fx,
    double *sig_curve_fx, double *beta, double spot_fx, double *mat_corr,
    double *corr_dom_for_ts, double *corr_dom_fx_ts, double *corr_for_fx_ts,
    long nb_mat_corr, char *dom_yc, char *for_yc, double *fx_vol,
    double disc_dt, double fx_dt)

{
  double *fx_vol_curve = NULL, *fx_fwd = NULL, *fx_time = NULL;

  long i;
  long nb_fx_time, nb_fx_time2;

  Err err = NULL;

  /* simple translation to make everything begin at 0 */

  nb_fx_time = (long)(end_date / fx_dt - 1.0E-08) + 1;
  if (nb_fx_time < 2) {
    nb_fx_time = 2;
  }

  fx_time = calloc(nb_fx_time, sizeof(double));

  if (!fx_time) {
    err = "Memory allocation error in Fx3DBetatsImpliedVol (1)";
    goto FREE_RETURN;
  }

  fx_time[0] = 0.0;

  /* fill the time where we want the volatility structure */
  for (i = 1; i < nb_fx_time - 1; i++) {
    fx_time[i] = fx_time[i - 1] + fx_dt;
  }

  fx_time[i] = end_date;

  num_f_concat_vector(&nb_fx_time, &fx_time, nbMat, maturity);
  num_f_concat_vector(&nb_fx_time, &fx_time, nb_mat_fx, mat_fx);
  num_f_sort_vector(nb_fx_time, fx_time);
  num_f_unique_vector(&nb_fx_time, fx_time);

  nb_fx_time2 = Get_Index(end_date, fx_time, nb_fx_time) + 1;

  fx_vol_curve = calloc(nb_fx_time2, sizeof(double));
  fx_fwd = calloc(nb_fx_time2, sizeof(double));

  if (!fx_vol || !fx_fwd) {
    err = "Memory allocation error in Fx3DBetatsImpliedVol (2)";
    goto FREE_RETURN;
  }

  err = Fxbeta_log_approx_corr(
      today, maturity, nbMat, sig_curve_dom, sig_curve_for, mat_fx, nb_mat_fx,
      sig_curve_fx, beta, lda_dom, lda_for, mat_corr, corr_dom_for_ts,
      corr_dom_fx_ts, corr_for_fx_ts, nb_mat_corr, spot_fx, dom_yc, for_yc,
      fx_time, nb_fx_time2, fx_fwd, fx_vol_curve, disc_dt);
  if (err) {
    goto FREE_RETURN;
  }

  err = Fx3DtsImpliedVol_corr(opt_maturity, start_date, end_date, maturity,
                              nbMat, sig_curve_dom, lda_dom, sig_curve_for,
                              lda_for, fx_time, fx_vol_curve, nb_fx_time2,
                              mat_corr, corr_dom_for_ts, corr_dom_fx_ts,
                              corr_for_fx_ts, nb_mat_corr, fx_vol);

FREE_RETURN:

  if (fx_vol_curve)
    free(fx_vol_curve);
  if (fx_time)
    free(fx_time);
  if (fx_fwd)
    free(fx_fwd);

  return err;
}

/*	Implied vol direct from underlying */
Err Fx3DBetaImpliedVol_corr(char *fx_underlying, double beta, double val_time,
                            double start_time, double end_time, double disc_dt,
                            double fx_dt, double *vol) {
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for, sigma_n_fx,
      nb_merge_dates, nb_corr;
  long i;
  double *sigma_date_dom = NULL, *sigma_dom = NULL;
  double *tau_date_dom = NULL, *tau_dom = NULL;
  double *sigma_date_for = NULL, *sigma_for = NULL;
  double *tau_date_for = NULL, *tau_for = NULL;
  double *sigma_date_fx = NULL, *sigma_fx = NULL;
  double *correl_date = NULL, *correl_dom_for = NULL, *correl_dom_fx = NULL,
         *correl_for_fx = NULL, *correl_dom_for_ts = NULL,
         *correl_dom_fx_ts = NULL, *correl_for_fx_ts = NULL;
  double lda_dom, lda_for;
  double *sig_dom = NULL, *sig_for = NULL, *sig_fx = NULL, *merge_dates = NULL;

  SrtUndPtr fx_und, dom_und, for_und;
  char *domname, *forname;

  double *beta2 = NULL;

  long today;
  double spot_fx;
  char *dom_yc, *for_yc;

  Err err = NULL;

  err = Get_FX_StochRate_TermStructures_corr(
      fx_underlying, &sigma_date_dom, &sigma_dom, &sigma_n_dom, &tau_date_dom,
      &tau_dom, &tau_n_dom, &sigma_date_for, &sigma_for, &sigma_n_for,
      &tau_date_for, &tau_for, &tau_n_for, &sigma_date_fx, &sigma_fx,
      &sigma_n_fx, &correl_date, &correl_dom_for, &correl_dom_fx,
      &correl_for_fx, &nb_corr);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = merge_rates_fx_corr_ts(
      sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for, sigma_for,
      sigma_n_for, sigma_date_fx, sigma_fx, sigma_n_fx, correl_date,
      correl_dom_for, correl_dom_fx, correl_for_fx, nb_corr, &merge_dates,
      &sig_dom, &sig_for, &sig_fx, &correl_dom_for_ts, &correl_dom_fx_ts,
      &correl_for_fx_ts, &nb_merge_dates);

  if (err) {
    goto FREE_RETURN;
  }

  beta2 = (double *)calloc(nb_merge_dates, sizeof(double));

  if (!beta2) {
    err = "Memory allocation error (1) in Fx3DBetaImpliedVol";
    goto FREE_RETURN;
  }

  for (i = nb_merge_dates - 1; i >= 0; i--) {
    beta2[i] = beta;
  }

  /* get today        , spot        , dom_yc        , for_yc */
  fx_und = lookup_und(fx_underlying);
  domname = get_domname_from_fxund(fx_und);
  forname = get_forname_from_fxund(fx_und);
  dom_und = lookup_und(domname);
  for_und = lookup_und(forname);

  err = get_underlying_discname(dom_und, &dom_yc);
  if (err) {
    goto FREE_RETURN;
  }
  err = get_underlying_discname(for_und, &for_yc);
  if (err) {
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(dom_und);

  err = Fx3DtsFwdFx(dom_und, for_und, fx_und, today, &spot_fx);

  err = Fx3DBetatsImpliedVol_corr(
      today, val_time, start_time, end_time, merge_dates, nb_merge_dates,
      sig_dom, lda_dom, sig_for, lda_for, merge_dates, nb_merge_dates, sig_fx,
      beta2, spot_fx, merge_dates, correl_dom_for, correl_dom_fx, correl_for_fx,
      nb_merge_dates, dom_yc, for_yc, vol, disc_dt, fx_dt);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sigma_date_fx) {
    free(sigma_date_fx);
  }

  if (sigma_fx) {
    free(sigma_fx);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (sig_fx) {
    free(sig_fx);
  }

  if (beta2) {
    free(beta2);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  if (correl_dom_for) {
    free(correl_dom_for);
  }

  if (correl_dom_fx) {
    free(correl_dom_fx);
  }

  if (correl_for_fx) {
    free(correl_for_fx);
  }

  if (correl_date) {
    free(correl_date);
  }

  if (correl_dom_for_ts) {
    free(correl_dom_for_ts);
  }

  if (correl_dom_fx_ts) {
    free(correl_dom_fx_ts);
  }

  if (correl_for_fx_ts) {
    free(correl_for_fx_ts);
  }

  return err;
}

Err Fx3DBetatsCalibration_corr(
    long today, double *exercise_opt, double *maturity_opt, double *vol_opt,
    long nbrOpt, double *maturity, long nbrMat, double *sig_curve_dom,
    double lda_dom, double *sig_curve_for, double lda_for, double *beta,
    double spot_fx, double *mat_corr, double *corr_dom_for_ts,
    double *corr_dom_fx_ts, double *corr_for_fx_ts, long nb_mat_corr,
    char *dom_yc, char *for_yc, double **fx_vol_curve, double disc_dt,
    double fx_dt, long nbIterMax)

{
  long i, nIter;
  double T, mat, Fwd;
  double bet;
  double *fx_vol_log = NULL;

  double vol_imp1, vol_imp2;
  double vol_tgt;
  double error;

  Err err = NULL;

  (*fx_vol_curve) = NULL;
  (*fx_vol_curve) = (double *)calloc(nbrOpt, sizeof(double));

  if (!(*fx_vol_curve)) {
    err = "Memory allocation error in Fx3DBetatsCalibration";
    goto FREE_RETURN;
  }

  /* First we calibrate a lognormal model */
  err = Fx3DtsCalibration_corr(
      exercise_opt, maturity_opt, vol_opt, nbrOpt, maturity, nbrMat,
      sig_curve_dom, lda_dom, sig_curve_for, lda_for, mat_corr, corr_dom_for_ts,
      corr_dom_fx_ts, corr_for_fx_ts, nb_mat_corr, &fx_vol_log);
  if (err) {
    goto FREE_RETURN;
  }

  for (i = 0; i < nbrOpt; i++) {
    /* Initialisation: SigBeta * F^beta = SigLn * F */
    mat = maturity_opt[i];
    T = exercise_opt[i];
    Fwd = spot_fx * swp_f_df(today, today + T * 365.0, for_yc) /
          swp_f_df(today, today + T * 365.0, dom_yc);
    bet = beta[i] - 1.0;
    vol_tgt = vol_opt[i];
    (*fx_vol_curve)[i] = fx_vol_log[i] * exp(-bet * log(Fwd));

    err = Fx3DBetatsImpliedVol_corr(
        today, mat, 0.0, T, maturity, nbrMat, sig_curve_dom, lda_dom,
        sig_curve_for, lda_for, exercise_opt, nbrOpt, (*fx_vol_curve), beta,
        spot_fx, mat_corr, corr_dom_for_ts, corr_dom_fx_ts, corr_for_fx_ts,
        nb_mat_corr, dom_yc, for_yc, &vol_imp1, disc_dt, fx_dt);
    if (err) {
      goto FREE_RETURN;
    }

    error = fabs(vol_imp1 - vol_tgt);

    nIter = 0;

    while (nIter <= nbIterMax && error > MAX_ERROR) {
      (*fx_vol_curve)[i] += SHIFT_VOL;
      err = Fx3DBetatsImpliedVol_corr(
          today, mat, 0.0, T, maturity, nbrMat, sig_curve_dom, lda_dom,
          sig_curve_for, lda_for, exercise_opt, nbrOpt, (*fx_vol_curve), beta,
          spot_fx, mat_corr, corr_dom_for_ts, corr_dom_fx_ts, corr_for_fx_ts,
          nb_mat_corr, dom_yc, for_yc, &vol_imp2, disc_dt, fx_dt);

      if (err) {
        goto FREE_RETURN;
      }

      (*fx_vol_curve)[i] -= SHIFT_VOL;
      (*fx_vol_curve)[i] +=
          (vol_tgt - vol_imp1) * SHIFT_VOL / (vol_imp2 - vol_imp1);

      if ((*fx_vol_curve)[i] < 0) {
        err = serror("Cannot calibrate option %d", i);
        goto FREE_RETURN;
      }

      err = Fx3DBetatsImpliedVol_corr(
          today, mat, 0.0, T, maturity, nbrMat, sig_curve_dom, lda_dom,
          sig_curve_for, lda_for, exercise_opt, nbrOpt, (*fx_vol_curve), beta,
          spot_fx, mat_corr, corr_dom_for_ts, corr_dom_fx_ts, corr_for_fx_ts,
          nb_mat_corr, dom_yc, for_yc, &vol_imp1, disc_dt, fx_dt);

      if (err) {
        goto FREE_RETURN;
      }

      error = fabs(vol_imp1 - vol_tgt);

      nIter += 1;
    }
  }

FREE_RETURN:

  if (fx_vol_log)
    free(fx_vol_log);

  return err;
}

Err Fx3DBetaCalibration_corr(char *dom_underlying, char *for_underlying,
                             double spot_fx, double beta, double *mat_corr,
                             double *corr_dom_for_ts, double *corr_dom_fx_ts,
                             double *corr_for_fx_ts, long nb_mat_corr,
                             double *exercise_opt, double *maturity_opt,
                             double *vol_opt, long nbropt,
                             double **fx_vol_curve, double disc_dt,
                             double fx_dt, long nbIterMax)

{
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL;

  double lda_dom, lda_for;

  SrtUndPtr dom_und, for_und;

  double *beta2 = NULL;

  long today, i;
  char *dom_yc, *for_yc;

  Err err = NULL;

  err =
      Get_LGM_TermStructure(dom_underlying, &sigma_date_dom, &sigma_dom,
                            &sigma_n_dom, &tau_date_dom, &tau_dom, &tau_n_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err =
      Get_LGM_TermStructure(for_underlying, &sigma_date_for, &sigma_for,
                            &sigma_n_for, &tau_date_for, &tau_for, &tau_n_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  dom_und = lookup_und(dom_underlying);
  for_und = lookup_und(for_underlying);

  err = get_underlying_discname(dom_und, &dom_yc);
  if (err) {
    goto FREE_RETURN;
  }
  err = get_underlying_discname(for_und, &for_yc);
  if (err) {
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(dom_und);

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  if (err) {
    goto FREE_RETURN;
  }

  beta2 = (double *)calloc(nbropt, sizeof(double));
  for (i = 0; i < nbropt; i++) {
    beta2[i] = beta;
  }

  err = Fx3DBetatsCalibration_corr(
      today, exercise_opt, maturity_opt, vol_opt, nbropt, merge_dates,
      nb_merge_dates, sigma_dom, lda_dom, sigma_for, lda_for, beta2, spot_fx,
      mat_corr, corr_dom_for_ts, corr_dom_fx_ts, corr_for_fx_ts, nb_mat_corr,
      dom_yc, for_yc, fx_vol_curve, disc_dt, fx_dt, nbIterMax);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  if (beta2) {
    free(beta2);
  }

  return err;
}

Err Fx3DBetatsCalibration2_corr(
    long today, double *exercise_opt, double *maturity_opt, double *vol_opt,
    long nbrOpt, long nbrLong, double *maturity, long nbrMat,
    double *sig_curve_dom, double lda_dom, double *sig_curve_for,
    double lda_for, double beta, double spot_fx, double *corr_mat,
    double *corr_dom_for_ts, double *corr_dom_fx_ts, double *corr_for_fx_ts,
    long nb_corr_mat, char *dom_yc, char *for_yc, double **fx_vol_curve,
    long nbSteps, long nbNewton, double disc_dt, double fx_dt, long nbIterMax) {
  long i, j;
  long nbrShort;
  long maturity_date, exercise_date, spot_date;
  double Fwd;
  double strikes[1];
  double df;
  double price1[2], error;
  double shift, vol_imp;
  double *beta2 = NULL;
  double delta;

  Err err = NULL;

  if (fabs(beta - 1.0) < 1.0E-04) {
    return Fx3DtsCalibration_corr(exercise_opt, maturity_opt, vol_opt, nbrOpt,
                                  maturity, nbrMat, sig_curve_dom, lda_dom,
                                  sig_curve_for, lda_for, corr_mat,
                                  corr_dom_for_ts, corr_dom_fx_ts,
                                  corr_for_fx_ts, nb_corr_mat, fx_vol_curve);
  }

  nbrShort = nbrOpt - nbrLong;

  spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  /* calculation of the spot spot fx */
  spot_fx *=
      swp_f_df(today, spot_date, dom_yc) / swp_f_df(today, spot_date, for_yc);

  beta2 = (double *)calloc(nbrOpt, sizeof(double));
  if (!beta2) {
    goto FREE_RETURN;
  }

  for (i = 0; i < nbrOpt; i++) {
    beta2[i] = beta;
  }

  err = Fx3DBetatsCalibration_corr(
      today, exercise_opt, maturity_opt, vol_opt, nbrOpt, maturity, nbrMat,
      sig_curve_dom, lda_dom, sig_curve_for, lda_for, beta2, spot_fx, corr_mat,
      corr_dom_for_ts, corr_dom_fx_ts, corr_for_fx_ts, nb_corr_mat, dom_yc,
      for_yc, fx_vol_curve, disc_dt, fx_dt, nbIterMax);
  if (err) {
    goto FREE_RETURN;
  }

  for (i = nbrShort; i < nbrOpt; i++) {
    smessage("Calibration of the option %d\n", i + 1);

    maturity_date = (long)(today + maturity_opt[i] * 365.0000000001);
    exercise_date = (long)(today + exercise_opt[i] * 365.0000000001);
    df = swp_f_df(today, exercise_date, dom_yc);
    Fwd = spot_fx * swp_f_df(today, exercise_date, for_yc) /
          swp_f_df(today, exercise_date, dom_yc);

    /* prices in the tree */
    strikes[0] = Fwd;

    error = 1.0E10;
    j = 0;

    while (j < nbNewton && (error > MAX_ERROR_VOL)) {
      err = Fx3DBetatsTreeFxOptions_corr(
          today, exercise_date, strikes, 1, maturity, nbrMat, sig_curve_dom,
          lda_dom, sig_curve_for, lda_for, exercise_opt, nbrOpt, *fx_vol_curve,
          0.0, beta, spot_fx, corr_mat, corr_dom_for_ts, corr_dom_fx_ts,
          corr_for_fx_ts, nb_corr_mat, dom_yc, for_yc, &(price1[0]), nbSteps);
      if (err) {
        goto FREE_RETURN;
      }

      /* calculates the model implied vol */
      err = srt_f_optimpvol(price1[0], Fwd, Fwd, exercise_opt[i], df, 0, 0,
                            &vol_imp);

      if (err) {
        shift = 0.0;
        nbSteps = (long)(nbSteps * 1.2);
      } else {
        error = fabs(vol_imp - vol_opt[i]);

        if (error > MAX_ERROR_VOL) {
          if (i > 0) {
            delta = (*fx_vol_curve)[i] * (*fx_vol_curve)[i] *
                        exp(2 * (beta2[i] - 1) * log(Fwd)) -
                    (vol_imp * vol_imp - vol_opt[i] * vol_opt[i]) *
                        exercise_opt[i] /
                        (exercise_opt[i] - exercise_opt[i - 1]);
          } else {
            delta = (*fx_vol_curve)[i] * (*fx_vol_curve)[i] *
                        exp(2 * (beta2[i] - 1) * log(Fwd)) -
                    (vol_imp * vol_imp - vol_opt[i] * vol_opt[i]);
          }

          if (delta > 0) {
            shift = -(*fx_vol_curve)[i] +
                    sqrt(delta) * exp((1.0 - beta2[i]) * log(Fwd));
            /*
            for (k=i; k<nbrOpt; k++)
            {
                    (*fx_vol_curve)[k] += shift ;
            }
            */
            (*fx_vol_curve)[i] += shift;
          } else {
            shift = 0.0;
            nbSteps = (long)(nbSteps * 1.2);
          }
        }
      }
      j += 1;
    }
  }

FREE_RETURN:

  if (beta2)
    free(beta2);
  return err;
}

Err Fx3DBetatsTreeFxOptions_corr(
    long today, long maturity_date, double *strikes, long nbrOpt,
    double *maturity, long nbrMat, double *sig_curve_dom, double dom_lam,
    double *sig_curve_for, double for_lam, double *maturity_fx, long nbrMat_fx,
    double *sig_curve_fx, double alpha, double beta, double spot_fx,
    double *corr_mat, double *corr_dom_for_ts, double *corr_dom_fx_ts,
    double *corr_for_fx_ts, long nb_corr_mat, char *dom_yc, char *for_yc,
    double *option_prices, long num_stp) {
  double maturity_opt;
  double *time = NULL, *date = NULL;
  long nstp = 1;
  long i, j;

  int *vol_change = NULL;
  double *dom_vol = NULL, *for_vol = NULL, *fx_vol = NULL;

  double *dom_ifr = NULL, *dom_fwd = NULL, *dom_var = NULL, *for_ifr = NULL,
         *for_fwd = NULL, *for_var = NULL, *fx_fwd = NULL, *fx_var = NULL,
         *corr_dom_for = NULL, *corr_dom_fx = NULL, *corr_for_fx = NULL;

  double *sig_fx_approx = NULL, *fx_fwd_approx = NULL, *beta_tab = NULL;

  double **func_parm = NULL;
  int *has_evt = NULL;

  int *is_bar = NULL;
  double *bar_lvl = NULL;
  int *bar_cl = NULL;

  double beta2 = 1.0 - beta;

  Err err = NULL;

  smessage("Phase 1 -preprocessing        , time in sec: 0");

  maturity_opt = (maturity_date - today) / DAYS_IN_YEAR;

  time = (double *)calloc(nstp, sizeof(double));
  if (!time) {
    err = "Memory allocation error (1) in Fx3DBetatsTreeFxOptions";
    goto FREE_RETURN;
  }

  time[0] = maturity_opt;

  /*	Fill the time vector */
  err = fill_time_vector(&time, &nstp, 0, NULL, 0, NULL, num_stp);

  if (err) {
    goto FREE_RETURN;
  }

  date = (double *)calloc(nstp, sizeof(double));
  has_evt = (int *)calloc(nstp, sizeof(int));
  func_parm = dmatrix(0, nstp - 1, 0, nbrOpt - 1);

  if (!date || !func_parm || !has_evt) {
    err = "Memory allocation error (2) in Fx3DBetatsTreeFxOptions";
    goto FREE_RETURN;
  }

  for (j = 0; j < nbrOpt; j++) {
    func_parm[nstp - 1][j] = strikes[j];
  }

  for (i = 0; i < nstp; i++) {
    date[i] = today + DAYS_IN_YEAR * time[i];
    has_evt[i] = 0;

    if (i > 0 && date[i] - date[i - 1] >= 1) {
      date[i] = (long)(date[i] + 1.0e-08);
      time[i] = YEARS_IN_DAY * (date[i] - today);
    }
  }
  has_evt[nstp - 1] = 1;

  /*	Fill the model parameters as required by tree_main_3dfx */

  vol_change = (int *)calloc(nstp, sizeof(int));
  dom_vol = (double *)calloc(nstp, sizeof(double));
  for_vol = (double *)calloc(nstp, sizeof(double));
  fx_vol = (double *)calloc(nstp, sizeof(double));
  corr_dom_for = (double *)calloc(nstp, sizeof(double));
  corr_dom_fx = (double *)calloc(nstp, sizeof(double));
  corr_for_fx = (double *)calloc(nstp, sizeof(double));

  beta_tab = (double *)calloc(nstp, sizeof(double));

  is_bar = (int *)calloc(nstp, sizeof(int));
  bar_lvl = (double *)calloc(nstp, sizeof(double));
  bar_cl = (int *)calloc(nstp, sizeof(int));

  if (!vol_change || !dom_vol || !for_vol || !fx_vol || !beta_tab || !is_bar ||
      !bar_lvl || !bar_cl || !corr_dom_for || !corr_dom_fx || !corr_for_fx) {
    err = "Memory allocation error (3) in Fx3DBetatsTreeFxOptions_corr";
    goto FREE_RETURN;
  }

  vol_change[nstp - 1] = 1;

  dom_vol[nstp - 1] =
      sig_curve_dom[Get_Index(time[nstp - 1], maturity, nbrMat)];
  for_vol[nstp - 1] =
      sig_curve_for[Get_Index(time[nstp - 1], maturity, nbrMat)];
  fx_vol[nstp - 1] =
      sig_curve_fx[Get_Index(time[nstp - 1], maturity_fx, nbrMat_fx)] *
      exp(alpha * sqrt(time[nstp - 1]) - 0.5 * alpha * alpha * time[nstp - 1]);
  ;
  corr_dom_for[nstp - 1] =
      corr_dom_for_ts[Get_Index(time[nstp - 1], corr_mat, nb_corr_mat)];
  corr_dom_fx[nstp - 1] =
      corr_dom_fx_ts[Get_Index(time[nstp - 1], corr_mat, nb_corr_mat)];
  corr_for_fx[nstp - 1] =
      corr_for_fx_ts[Get_Index(time[nstp - 1], corr_mat, nb_corr_mat)];

  beta_tab[nstp - 1] = beta;

  for (i = nstp - 2; i >= 0; i--) {
    dom_vol[i] = sig_curve_dom[Get_Index(time[i], maturity, nbrMat)];
    for_vol[i] = sig_curve_for[Get_Index(time[i], maturity, nbrMat)];
    fx_vol[i] = sig_curve_fx[Get_Index(time[i], maturity_fx, nbrMat_fx)] *
                exp(alpha * sqrt(time[i]) - 0.5 * alpha * alpha * time[i]);

    corr_dom_for[i] =
        corr_dom_for_ts[Get_Index(time[i], corr_mat, nb_corr_mat)];
    corr_dom_fx[i] = corr_dom_fx_ts[Get_Index(time[i], corr_mat, nb_corr_mat)];
    corr_for_fx[i] = corr_for_fx_ts[Get_Index(time[i], corr_mat, nb_corr_mat)];

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
    err = "Memory allocation error (3) in Fx3DBetatsTreeFxOptions";
    goto FREE_RETURN;
  }

  /* first get the coresponding lognormal volatilities */
  err = Fxbeta_log_approx_corr(today, time, nstp, dom_vol, for_vol, time, nstp,
                               fx_vol, beta_tab, dom_lam, for_lam, corr_mat,
                               corr_dom_for_ts, corr_dom_fx_ts, corr_for_fx_ts,
                               nb_corr_mat, spot_fx, dom_yc, for_yc, time, nstp,
                               fx_fwd_approx, sig_fx_approx, MAX_TIME);

  fill_fwd_var_corr(nstp, time, date, dom_vol, for_vol, sig_fx_approx, dom_lam,
                    for_lam, corr_dom_for, corr_dom_fx, corr_for_fx, dom_yc,
                    for_yc, dom_ifr, dom_fwd, dom_var, for_ifr, for_fwd,
                    for_var, fx_fwd, fx_var);

  for (i = nstp - 1; i >= 0; i--) {
    bar_lvl[i] = (exp(beta2 * (fx_fwd_approx[i])) - 1.0) / beta2;
  }

  err = treeBeta_main_3dfx(
      nstp, time, date, vol_change, dom_vol, for_vol, fx_vol, beta_tab, dom_ifr,
      dom_fwd, dom_var, for_ifr, for_fwd, for_var, fx_fwd, fx_var, func_parm,
      has_evt, bar_lvl, bar_cl, is_bar, dom_lam, for_lam, corr_dom_for,
      corr_dom_fx, corr_for_fx, spot_fx, dom_yc, for_yc,
      FxCall_payoff_4_3dfxBeta_tree, nbrOpt, 1, option_prices);

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

  if (bar_lvl)
    free(bar_lvl);
  if (is_bar)
    free(is_bar);
  if (bar_cl)
    free(bar_cl);

  if (func_parm) {
    free_dmatrix(func_parm, 0, nstp - 1, 0, nbrOpt - 1);
  }

  return err;
}

Err Fx3DBetaCalibration2_corr(char *dom_underlying, char *for_underlying,
                              double spot_fx, double beta, double *correl_mat,
                              double *correl_dom_for, double *correl_dom_fx,
                              double *correl_for_fx, long nb_correl,
                              double *exercise_opt, double *maturity_opt,
                              double *vol_opt, long nbropt, long nbLong,
                              double **fx_vol_curve, long nbSteps,
                              long nbNewton, double disc_dt, double fx_dt,
                              long nbIterMax)

{
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL;

  double lda_dom, lda_for;

  SrtUndPtr dom_und, for_und;

  long today;
  char *dom_yc, *for_yc;

  Err err = NULL;

  err =
      Get_LGM_TermStructure(dom_underlying, &sigma_date_dom, &sigma_dom,
                            &sigma_n_dom, &tau_date_dom, &tau_dom, &tau_n_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err =
      Get_LGM_TermStructure(for_underlying, &sigma_date_for, &sigma_for,
                            &sigma_n_for, &tau_date_for, &tau_for, &tau_n_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  dom_und = lookup_und(dom_underlying);
  for_und = lookup_und(for_underlying);

  err = get_underlying_discname(dom_und, &dom_yc);
  if (err) {
    goto FREE_RETURN;
  }
  err = get_underlying_discname(for_und, &for_yc);
  if (err) {
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(dom_und);

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  if (err) {
    goto FREE_RETURN;
  }

  err = Fx3DBetatsCalibration2_corr(
      today, exercise_opt, maturity_opt, vol_opt, nbropt, nbLong, merge_dates,
      nb_merge_dates, sig_dom, lda_dom, sig_for, lda_for, beta, spot_fx,
      correl_mat, correl_dom_for, correl_dom_fx, correl_for_fx, nb_correl,
      dom_yc, for_yc, fx_vol_curve, nbSteps, nbNewton, disc_dt, fx_dt,
      nbIterMax);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  return err;
}

Err Fx3DAlphaBetatsTreeFxOptions_corr(
    long today, long maturity_date, double *strikes, long nbrOpt,
    double *maturity, long nbrMat, double *sig_curve_dom, double dom_lam,
    double *sig_curve_for, double for_lam, double *maturity_fx, long nbrMat_fx,
    double *sig_curve_fx, double alpha, double beta, double spot_fx,
    double *corr_mat, double *corr_dom_for, double *corr_dom_fx,
    double *corr_for_fx, long nb_corr, char *dom_yc, char *for_yc,
    double *option_prices, long num_stp) {
  char *cerr = NULL;
  double *pv1 = NULL, *pv2 = NULL;
  int i;

  pv1 = (double *)calloc(nbrOpt + 1, sizeof(double));
  if (!pv1) {
    cerr = "Memory allocation error in Fx3DAlphaBetatsTreeFxOptions";
    goto FREE_RETURN;
  }

  pv2 = (double *)calloc(nbrOpt + 1, sizeof(double));
  if (!pv2) {
    cerr = "Memory allocation error in Fx3DAlphaBetatsTreeFxOptions";
    goto FREE_RETURN;
  }

  if (fabs(alpha) <= 1.0e-03) {
    alpha = 0.0;
  }

  cerr = Fx3DBetatsTreeFxOptions_corr(
      today, maturity_date, strikes, nbrOpt, maturity, nbrMat, sig_curve_dom,
      dom_lam, sig_curve_for, for_lam, maturity_fx, nbrMat_fx, sig_curve_fx,
      alpha, beta, spot_fx, corr_mat, corr_dom_for, corr_dom_fx, corr_for_fx,
      nb_corr, dom_yc, for_yc, pv1, num_stp);

  if (cerr) {
    goto FREE_RETURN;
  }

  if (fabs(alpha) > 1.0e-03) {
    cerr = Fx3DBetatsTreeFxOptions_corr(
        today, maturity_date, strikes, nbrOpt, maturity, nbrMat, sig_curve_dom,
        dom_lam, sig_curve_for, for_lam, maturity_fx, nbrMat_fx, sig_curve_fx,
        -alpha, beta, spot_fx, corr_mat, corr_dom_for, corr_dom_fx, corr_for_fx,
        nb_corr, dom_yc, for_yc, pv2, num_stp);

    if (cerr) {
      goto FREE_RETURN;
    }
  } else {
    memcpy(pv2, pv1, nbrOpt * sizeof(double));
  }

  for (i = 0; i < nbrOpt; i++) {
    option_prices[i] = 0.5 * (pv1[i] + pv2[i]);
  }

FREE_RETURN:

  if (pv1) {
    free(pv1);
  }

  if (pv2) {
    free(pv2);
  }

  return cerr;
}

Err Fx3DAlphaBetatsCalibration2_corr(
    long today, double *exercise_opt, double *maturity_opt, double *vol_opt,
    long nbrOpt, long nbrLong, double *maturity, long nbrMat,
    double *sig_curve_dom, double lda_dom, double *sig_curve_for,
    double lda_for, double alpha, double beta, double spot_fx, double *corr_mat,
    double *corr_dom_for, double *corr_dom_fx, double *corr_for_fx,
    long nb_corr, char *dom_yc, char *for_yc, double **fx_vol_curve,
    long nbSteps, long nbNewton, double disc_dt, double fx_dt, long nbIterMax)

{
  long i, j;
  long nbrShort;
  long maturity_date, exercise_date, spot_date;
  double Fwd;
  double strikes[1];
  double df;
  double price1[2], error;
  double shift, vol_imp;
  double *beta2 = NULL;
  double delta;

  Err err = NULL;

  if (alpha < 1.0E-04) {
    return Fx3DBetatsCalibration2_corr(
        today, exercise_opt, maturity_opt, vol_opt, nbrOpt, nbrLong, maturity,
        nbrMat, sig_curve_dom, lda_dom, sig_curve_for, lda_for, beta, spot_fx,
        corr_mat, corr_dom_for, corr_dom_fx, corr_for_fx, nb_corr, dom_yc,
        for_yc, fx_vol_curve, nbSteps, nbNewton, disc_dt, fx_dt, nbIterMax);
  }

  nbrShort = nbrOpt - nbrLong;

  spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  /* calculation of the spot spot fx */
  spot_fx *=
      swp_f_df(today, spot_date, dom_yc) / swp_f_df(today, spot_date, for_yc);

  beta2 = (double *)calloc(nbrOpt, sizeof(double));
  if (!beta2) {
    goto FREE_RETURN;
  }

  for (i = 0; i < nbrOpt; i++) {
    beta2[i] = beta;
  }

  err = Fx3DBetatsCalibration_corr(
      today, exercise_opt, maturity_opt, vol_opt, nbrOpt, maturity, nbrMat,
      sig_curve_dom, lda_dom, sig_curve_for, lda_for, beta2, spot_fx, corr_mat,
      corr_dom_for, corr_dom_fx, corr_for_fx, nb_corr, dom_yc, for_yc,
      fx_vol_curve, disc_dt, fx_dt, nbIterMax);
  if (err) {
    goto FREE_RETURN;
  }

  for (i = nbrShort; i < nbrOpt; i++) {
    smessage("Calibration of the option %d\n", i + 1);

    maturity_date = (long)(today + maturity_opt[i] * 365.0000000001);
    exercise_date = (long)(today + exercise_opt[i] * 365.0000000001);
    df = swp_f_df(today, exercise_date, dom_yc);
    Fwd = spot_fx * swp_f_df(today, exercise_date, for_yc) /
          swp_f_df(today, exercise_date, dom_yc);

    /* prices in the tree */
    strikes[0] = Fwd;

    error = 1.0E10;
    j = 0;

    while (j < nbNewton && (error > MAX_ERROR_VOL)) {
      err = Fx3DAlphaBetatsTreeFxOptions_corr(
          today, exercise_date, strikes, 1, maturity, nbrMat, sig_curve_dom,
          lda_dom, sig_curve_for, lda_for, exercise_opt, nbrOpt, *fx_vol_curve,
          alpha, beta, spot_fx, corr_mat, corr_dom_for, corr_dom_fx,
          corr_for_fx, nb_corr, dom_yc, for_yc, &(price1[0]), nbSteps);
      if (err) {
        goto FREE_RETURN;
      }

      /* calculates the model implied vol */
      err = srt_f_optimpvol(price1[0], Fwd, Fwd, exercise_opt[i], df, 0, 0,
                            &vol_imp);

      if (err) {
        shift = 0.0;
        nbSteps = (long)(nbSteps * 1.2);
      } else {
        error = fabs(vol_imp - vol_opt[i]);

        if (error > MAX_ERROR_VOL) {
          if (i > 0) {
            delta = (*fx_vol_curve)[i] * (*fx_vol_curve)[i] *
                        exp(2 * (beta2[i] - 1) * log(Fwd)) -
                    (vol_imp * vol_imp - vol_opt[i] * vol_opt[i]) *
                        exercise_opt[i] /
                        (exercise_opt[i] - exercise_opt[i - 1]);
          } else {
            delta = (*fx_vol_curve)[i] * (*fx_vol_curve)[i] *
                        exp(2 * (beta2[i] - 1) * log(Fwd)) -
                    (vol_imp * vol_imp - vol_opt[i] * vol_opt[i]);
          }

          if (delta > 0) {
            shift = -(*fx_vol_curve)[i] +
                    sqrt(delta) * exp((1.0 - beta2[i]) * log(Fwd));
            /*
            for (k=i; k<nbrOpt; k++)
            {
                    (*fx_vol_curve)[k] += shift ;
            }
            */
            (*fx_vol_curve)[i] += shift;
          } else {
            shift = 0.0;
            nbSteps = (long)(nbSteps * 1.2);
          }
        }
      }
      j += 1;
    }
  }

FREE_RETURN:

  if (beta2)
    free(beta2);
  return err;
}

Err Fx3DAlphaBetaCalibration2_corr(
    char *dom_underlying, char *for_underlying, double spot_fx, double alpha,
    double beta, double *corr_mat, double *correl_dom_for,
    double *correl_dom_fx, double *correl_for_fx, long nb_corr,
    double *exercise_opt, double *maturity_opt, double *vol_opt, long nbropt,
    long nbLong, double **fx_vol_curve, long nbSteps, long nbNewton,
    double disc_dt, double fx_dt, long nbIterMax)

{
  long sigma_n_dom, tau_n_dom, sigma_n_for, tau_n_for;
  long nb_merge_dates;
  double *sigma_date_dom = NULL, *sigma_dom = NULL, *tau_date_dom = NULL,
         *tau_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
         *tau_date_for = NULL, *tau_for = NULL, *merge_dates = NULL,
         *sig_dom = NULL, *sig_for = NULL;

  double lda_dom, lda_for;

  SrtUndPtr dom_und, for_und;

  long today;
  char *dom_yc, *for_yc;

  Err err = NULL;

  err =
      Get_LGM_TermStructure(dom_underlying, &sigma_date_dom, &sigma_dom,
                            &sigma_n_dom, &tau_date_dom, &tau_dom, &tau_n_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err =
      Get_LGM_TermStructure(for_underlying, &sigma_date_for, &sigma_for,
                            &sigma_n_for, &tau_date_for, &tau_for, &tau_n_for);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_dom, tau_n_dom, &lda_dom);
  if (err) {
    goto FREE_RETURN;
  }

  err = get_unique_lambda(tau_for, tau_n_for, &lda_for);
  if (err) {
    goto FREE_RETURN;
  }

  dom_und = lookup_und(dom_underlying);
  for_und = lookup_und(for_underlying);

  err = get_underlying_discname(dom_und, &dom_yc);
  if (err) {
    goto FREE_RETURN;
  }
  err = get_underlying_discname(for_und, &for_yc);
  if (err) {
    goto FREE_RETURN;
  }

  today = get_today_from_underlying(dom_und);

  err = merge_rates_ts(sigma_date_dom, sigma_dom, sigma_n_dom, sigma_date_for,
                       sigma_for, sigma_n_for, &merge_dates, &sig_dom, &sig_for,
                       &nb_merge_dates);

  if (err) {
    goto FREE_RETURN;
  }

  err = Fx3DAlphaBetatsCalibration2_corr(
      today, exercise_opt, maturity_opt, vol_opt, nbropt, nbLong, merge_dates,
      nb_merge_dates, sigma_dom, lda_dom, sigma_for, lda_for, alpha, beta,
      spot_fx, corr_mat, correl_dom_for, correl_dom_fx, correl_for_fx, nb_corr,
      dom_yc, for_yc, fx_vol_curve, nbSteps, nbNewton, disc_dt, fx_dt,
      nbIterMax);

FREE_RETURN:

  if (sigma_date_dom) {
    free(sigma_date_dom);
  }

  if (sigma_dom) {
    free(sigma_dom);
  }

  if (tau_date_dom) {
    free(tau_date_dom);
  }

  if (tau_dom) {
    free(tau_dom);
  }

  if (sigma_date_for) {
    free(sigma_date_for);
  }

  if (sigma_for) {
    free(sigma_for);
  }

  if (tau_date_for) {
    free(tau_date_for);
  }

  if (tau_for) {
    free(tau_for);
  }

  if (sig_dom) {
    free(sig_dom);
  }

  if (sig_for) {
    free(sig_for);
  }

  if (merge_dates) {
    free(merge_dates);
  }

  return err;
}