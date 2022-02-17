/* ------------------------------------------------------------------------
   FILENAME:  	srt_f_ts_irm_fct.cxx

   PURPOSE:     Provide a few utilities when dealing with a TermStruct:
                                Provides the interpolation functions for the TS
                                          components
   ------------------------------------------------------------------------ */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_ts.h"
#include "srt_h_ts_irm.h"
#include "srtaccess.h"
#include <opfnctns.H>

#define ZERO_SHIFT 1.0e-6
#define PROP_SHIFT 1.0e-3

Err srt_f_get_vasicek_mean_sr(double time, TermStruct *ts,
                              double *vasicek_mean_sr) {
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;
  Err err = NULL;

  ls = ts->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls == NULL)
    ls = ts->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != ts->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  (*vasicek_mean_sr) =
      tsval->vasicek_mean_sr +
      (tsval->mean_rev_level / tsval->F) * (exp(time / tsval->tau) - 1);

  return err;
}

Err srt_f_get_vasicek_mean_int_sr(double time, TermStruct *ts,
                                  double *vasicek_mean_int_sr) {
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;
  Err err = NULL;
  double exp_term;

  ls = ts->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls == NULL)
    ls = ts->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != ts->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  exp_term = (1 - exp(-time / tsval->tau)) * tsval->tau;

  (*vasicek_mean_int_sr) =
      tsval->vasicek_mean_int_sr +
      (tsval->vasicek_mean_sr * tsval->F - tsval->mean_rev_level) * exp_term +
      tsval->mean_rev_level * time;

  return err;
}
Err srt_f_get_vasicek_var_int_sr(double time, TermStruct *ts,
                                 double *vasicek_var_int_sr) {
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;
  double dom_psi;

  dom_psi = Psi_func(time, ts);

  ls = ts->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = ts->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != ts->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  (*vasicek_var_int_sr) = 0.0;

  (*vasicek_var_int_sr) =
      (tsval->L +
       tsval->F * tsval->G * (1 - exp(-time / tsval->tau)) * tsval->tau +
       0.5 * pow(tsval->sig * tsval->tau, 2) * (exp(-time / tsval->tau) - 2) /
           tsval->F);

  (*vasicek_var_int_sr) *= dom_psi;

  (*vasicek_var_int_sr) -=
      (tsval->O +
       tsval->F * tsval->G * tsval->Psi * (1 - exp(-time / tsval->tau)) *
           tsval->tau +
       tsval->F * tsval->F * tsval->G * tsval->tau *
           ((1 - exp(-time / tsval->tau)) * tsval->tau -
            0.5 * (1 - exp(-2 * time / tsval->tau)) * tsval->tau) -
       0.5 * (tsval->Psi / tsval->F) * tsval->tau * tsval->sig * tsval->sig *
           (1 - exp(-time / tsval->tau)) * tsval->tau -
       0.5 * pow(tsval->sig * tsval->tau, 2) * (1 - exp(-time / tsval->tau)) *
           tsval->tau +
       0.5 * pow(tsval->sig * tsval->tau, 2) *
           (0.5 * (1 - exp(-2 * time / tsval->tau)) * tsval->tau - time));

  (*vasicek_var_int_sr) +=
      0.5 * pow(tsval->sig * tsval->tau, 2) * tsval->Psi / tsval->F;

  (*vasicek_var_int_sr) *= 2.0;

  return NULL;
}

Err srt_f_get_vasicek_init_cond(double time, TermStruct *ts,
                                double *vasicek_init_cond) {

  SrtLst *ls;

  ls = ts->head;

  if (ls != NULL)
    (*vasicek_init_cond) =
        ((IrmTermStructVal *)ls->element->val.pval)->vasicek_init_cond;

  else
    return serror("the vasicek term struct. is not initialised ");

  return NULL;
}

static Err srt_f_get_lgm_mean_int_sr(double time, TermStruct *ts,
                                     double *mean_int_nom_sr) {
  (*mean_int_nom_sr) = 0.0;
  (*mean_int_nom_sr) =
      -Psi_func(time, ts) * L_func(time, ts) + O_func(time, ts);

  return NULL;
}

static Err srt_f_get_lgm_var_int_sr(double time, TermStruct *ts,
                                    double *var_int_nom_sr) {
  (*var_int_nom_sr) =
      2 * (Psi_func(time, ts) * L_func(time, ts) - O_func(time, ts));

  return NULL;
}

Err srt_f_vasicek_discount_factor(Date discount_mat_date, SrtUndPtr und,
                                  double *discount_factor) {
  TermStruct *ts;
  Date date_today;
  double discount_factor_mat;
  double vasicek_init_cond, var_int_sr, vasicek_mean_int_sr;
  double Psi_at_T;
  Err err = NULL;

  /* get today from the underlying */
  date_today = get_today_from_underlying(und);

  /* get the discount factor maturity */
  discount_factor_mat = (discount_mat_date - date_today) * YEARS_IN_DAY;

  /* get the underlying term structure */
  err = get_underlying_ts(und, &ts);
  if (err)
    return err;

  err = srt_f_get_vasicek_var_int_sr(discount_factor_mat, ts, &var_int_sr);
  if (err)
    return err;

  /* first compute the variance of the integral of the short rate */
  Psi_at_T = Psi_func(discount_factor_mat, ts);

  /* compute the mean of the integral of the short rate */

  err =
      srt_f_get_vasicek_init_cond(discount_factor_mat, ts, &vasicek_init_cond);
  if (err)
    return err;

  err = srt_f_get_vasicek_mean_int_sr(discount_factor_mat, ts,
                                      &vasicek_mean_int_sr);
  if (err)
    return err;

  vasicek_mean_int_sr += vasicek_init_cond * Psi_at_T;

  (*discount_factor) = exp(-vasicek_mean_int_sr + 0.5 * var_int_sr);

  return err;
}

Err srt_f_vasicek_fra(Date start_date, Date end_date, SrtUndPtr und,
                      double *fra) {
  Err err = NULL;
  double df_end, df_start;
  double fra_cash;

  err = srt_f_vasicek_discount_factor(start_date, und, &df_start);
  if (err)
    return err;

  err = srt_f_vasicek_discount_factor(end_date, und, &df_end);
  if (err)
    return err;

  fra_cash = df_start / df_end - 1.0;
  fra_cash /= coverage(DTOL(start_date), DTOL(end_date), BASIS_ACT_360);

  (*fra) = fra_cash;

  return err;
}

Err srt_f_vasicek_swap(Date start_date, Date theo_end_date, char *swap_freq_str,
                       char *swap_basis_str, SrtUndPtr und, double *swap) {
  Err err = NULL;
  SwapDP swp_dp;
  SrtCompounding swap_freq_code;
  BasisCode swap_basis_code;
  Date date_today, *pay_dates, *start_dates, *end_dates;
  long num_pay_dates, num_dates, i;
  double *coverages, swp_level, df_start, df_end, df;

  /* Get the compounding and basis codes */
  err = interp_compounding(swap_freq_str, &swap_freq_code);
  if (err)
    return err;

  err = interp_basis(swap_basis_str, &swap_basis_code);
  if (err)
    return err;

  /* Set the SwapDP */
  err = swp_f_setSwapDP(start_date, theo_end_date, swap_freq_code,
                        swap_basis_code, &swp_dp);
  if (err)
    return err;

  /* Get today from underlying */
  date_today = get_today_from_underlying(und);

  /* Make the fixed leg and the coverages */
  err = swp_f_make_FixedLegDatesAndCoverages(
      &swp_dp, date_today, &pay_dates, &num_pay_dates, &start_dates, &end_dates,
      &coverages, &num_dates);
  if (err)
    return err;

  /* Get the dicount factor at start and end date */
  err = srt_f_vasicek_discount_factor(start_date, und, &df_start);
  if (err)
    return err;

  err = srt_f_vasicek_discount_factor(pay_dates[num_dates], und, &df_end);
  if (err)
    return err;

  swp_level = 0;
  for (i = 1; i <= num_dates; i++) {
    err = srt_f_vasicek_discount_factor(pay_dates[i], und, &df);

    if (err)
      return err;

    swp_level += coverages[i - 1] * df;
  }

  (*swap) = (df_start - df_end) / swp_level;

  return err;
}

Err srt_f_vasicek_fut(double start_date, double pay_date,
                      char *vasicek_und_name, double *fut_pr) {

  Err err = NULL;
  SrtUndPtr und;
  TermStruct *und_ts;
  double date_today, start_time, pay_time;
  double G_start, H_start, Psi_start_pay;
  double df_start, df_pay;
  double adjustement;

  und = lookup_und(vasicek_und_name);
  if (!und)
    return serror("can not find underlying");

  err = get_underlying_ts(und, &und_ts);
  if (err)
    return err;

  date_today = get_today_from_underlying(und);

  start_time = (start_date - date_today) * YEARS_IN_DAY;
  pay_time = (pay_date - date_today) * YEARS_IN_DAY;

  G_H_func(start_time, und_ts, &G_start, &H_start);

  Psi_start_pay = Psi_func(pay_time, und_ts) - Psi_func(start_time, und_ts);

  adjustement = 0.0;
  adjustement += Psi_start_pay * Psi_start_pay * G_start +
                 Psi_start_pay * L_func(start_time, und_ts);

  err = srt_f_vasicek_discount_factor((long)start_date, und, &df_start);
  if (err)
    return err;

  err = srt_f_vasicek_discount_factor((long)pay_date, und, &df_pay);
  if (err)
    return err;

  (*fut_pr) = ((df_start / df_pay) * exp(adjustement) - 1) /
              coverage((Date)start_date, (Date)pay_date, BASIS_ACT_360);

  return err;
}

Err srt_f_vasicek_market_libor_expectation(double start_date, double pay_date,
                                           double risk_premuim,
                                           char *vasicek_und_name,
                                           double *market_libor_expectation) {
  Err err = NULL;
  SrtUndPtr und;
  TermStruct *und_ts;
  double date_today, start_time, pay_time;
  double G_start, H_start, Psi_start_pay;
  double df_start, df_pay;
  double adjustement;

  und = lookup_und(vasicek_und_name);
  if (!und)
    return serror("can not find underlying");

  err = get_underlying_ts(und, &und_ts);
  if (err)
    return err;

  date_today = get_today_from_underlying(und);

  start_time = (start_date - date_today) * YEARS_IN_DAY;
  pay_time = (pay_date - date_today) * YEARS_IN_DAY;

  G_H_func(start_time, und_ts, &G_start, &H_start);

  Psi_start_pay = Psi_func(pay_time, und_ts) - Psi_func(start_time, und_ts);

  adjustement = 0.0;
  adjustement += Psi_start_pay * Psi_start_pay * G_start +
                 Psi_start_pay * L_func(start_time, und_ts);

  adjustement += risk_premuim * Psi_start_pay * J_func(start_time, und_ts);

  err = srt_f_vasicek_discount_factor((long)start_date, und, &df_start);
  if (err)
    return err;

  err = srt_f_vasicek_discount_factor((long)pay_date, und, &df_pay);
  if (err)
    return err;

  (*market_libor_expectation) =
      ((df_start / df_pay) * exp(adjustement) - 1) /
      coverage((Date)start_date, (Date)pay_date, BASIS_ACT_360);

  return err;
}

Err srt_f_vasicek_cont_fwd_zr(double fixing_date, double pay_date,
                              SrtUndPtr und, double *cont_fwd_zr) {
  Err err = NULL;
  double df_fixing, df_pay;

  err = srt_f_vasicek_discount_factor((long)fixing_date, und, &df_fixing);
  if (err)
    return err;

  err = srt_f_vasicek_discount_factor((long)pay_date, und, &df_pay);
  if (err)
    return err;

  (*cont_fwd_zr) =
      -log(df_pay / df_fixing) / ((pay_date - fixing_date) * YEARS_IN_DAY);

  return NULL;
}

Err srt_f_vasicek_dirty_pr(SrtUndPtr und, long num_dates, double settlt_date,
                           double *pay_dates, double *cpns,
                           SRT_Boolean discount_fwd_dirty_pr,
                           double *dirty_price) {
  Err err = NULL;
  long i;
  double df;

  if (num_dates == 0)
    return serror("bond does not bear coupon");

  (*dirty_price) = 0;
  for (i = 0; i < num_dates; i++) {
    err = srt_f_vasicek_discount_factor((long)pay_dates[i], und, &df);
    if (err)
      return err;

    (*dirty_price) += cpns[i] * df;
  }

  if (discount_fwd_dirty_pr) {
    err = srt_f_vasicek_discount_factor((long)settlt_date, und, &df);
    if (err)
      return err;

    (*dirty_price) /= df;
  }

  return NULL;
}

Err srt_f_vasicek_overnight_interest_rate_swap(
    char *und_name, long start_date, long end_date,

    double *overnight_interest_rate_swap) {
  Err err = NULL;
  Date date_today;
  double start_t, pay_T;
  double Psi_t, vasicek_mean_sr_t, vasicek_mean_int_sr_t,
      vasicek_mean_int_sr_pay_T;
  double G_t, H_t, G_pay_T, H_pay_T;
  double start_t_cond_var, start_t_cond_mean, fwd_measure_sr_mean,
      vasicek_init_cond, sr_coeff;
  double sr_variance;
  TermStruct *ts;
  SrtUndPtr und;

  und = lookup_und(und_name);
  if (!und)
    return serror("can not find the underlying");

  /* get the underlying term struct */
  err = get_underlying_ts(und, &ts);
  if (err)
    return err;

  /* get the date today from the underlying */
  date_today = get_today_from_underlying(und);

  start_t = (start_date - date_today) * YEARS_IN_DAY;
  pay_T = (end_date - date_today) * YEARS_IN_DAY;

  err = srt_f_get_vasicek_mean_int_sr(start_t, ts, &vasicek_mean_int_sr_t);
  if (err)
    return err;

  err = srt_f_get_vasicek_mean_int_sr(pay_T, ts, &vasicek_mean_int_sr_pay_T);
  if (err)
    return err;

  err = srt_f_get_vasicek_mean_sr(start_t, ts, &vasicek_mean_sr_t);
  if (err)
    return err;

  Psi_t = Psi_func(start_t, ts);
  G_H_func(start_t, ts, &G_t, &H_t);

  /* start_t_cond_mean is the start_date cond mean of the integral of the sr
   * (forward measure) */

  start_t_cond_mean = 0.0;

  start_t_cond_mean = (vasicek_mean_int_sr_pay_T - vasicek_mean_int_sr_t);
  start_t_cond_mean -= (Psi_func(pay_T, ts) - Psi_t) * vasicek_mean_sr_t;

  /* extra term du to the forward measure expectation */
  start_t_cond_mean -=
      Psi_func(pay_T, ts) * (I_func(pay_T, ts) - I_func(start_t, ts));
  start_t_cond_mean +=
      Psi_func(pay_T, ts) * G_t * (Psi_func(pay_T, ts) - Psi_func(start_t, ts));
  start_t_cond_mean -=
      S_func(start_t, ts) * (Psi_func(pay_T, ts) - Psi_func(start_t, ts));
  start_t_cond_mean += (T_func(pay_T, ts) - T_func(start_t, ts));

  /* start_t_cond_var is the start_date cond variance of the integral of the sr
   */
  start_t_cond_var = 0.0;
  G_H_func(pay_T, ts, &G_pay_T, &H_pay_T);

  start_t_cond_var =
      Psi_func(pay_T, ts) * Psi_func(pay_T, ts) * (G_pay_T - G_t);

  start_t_cond_var -=
      2 * Psi_func(pay_T, ts) * (S_func(pay_T, ts) - S_func(start_t, ts));
  start_t_cond_var += (G_pay_T * Psi_func(pay_T, ts) * Psi_func(pay_T, ts) -
                       2 * O_func(pay_T, ts));
  start_t_cond_var -= (G_t * Psi_t * Psi_t - 2 * O_func(start_t, ts));

  /* compute the mean of the sr (fwd measure) */

  err = srt_f_get_vasicek_init_cond(start_t, ts, &vasicek_init_cond);
  if (err)
    return err;

  fwd_measure_sr_mean = 0.0;

  fwd_measure_sr_mean += vasicek_init_cond * F_func(start_t, ts);
  fwd_measure_sr_mean += F_func(start_t, ts) * vasicek_mean_sr_t;
  fwd_measure_sr_mean -=
      F_func(start_t, ts) * Psi_func(start_t, ts) * (G_pay_T - G_t);
  fwd_measure_sr_mean +=
      F_func(start_t, ts) * (S_func(pay_T, ts) - S_func(start_t, ts));

  /* compute the short rate coeff. */

  sr_coeff = 0.0;
  sr_coeff =
      (Psi_func(pay_T, ts) - Psi_func(start_t, ts)) / F_func(start_t, ts);

  sr_variance = 0.0;
  sr_variance = F_func(start_t, ts) * F_func(start_t, ts) * G_t;

  (*overnight_interest_rate_swap) =
      exp(start_t_cond_mean + 0.5 * start_t_cond_var);

  (*overnight_interest_rate_swap) *= exp(
      sr_coeff * fwd_measure_sr_mean + 0.5 * sr_coeff * sr_coeff * sr_variance);

  (*overnight_interest_rate_swap) -= 1;

  (*overnight_interest_rate_swap) *=
      1.0 / coverage(start_date, end_date, BASIS_ACT_360);

  return NULL;
}

static Err one_bond_calib_func(
    char *model, SrtMdlType mdl_type, SrtMdlDim mdl_dim, char *vasicek_und_name,
    char *yc_name, double **sig_data, long num_sig, double **tau_data,
    long num_tau, double **vasicek_mean_rev_level_data,
    double vasicek_init_cond,

    double settlt_date, long num_bond, long *num_cpns, double **pay_dates,
    double **cpns,

    long bond_index,

    double new_mean_rev_level, SRT_Boolean discount_fwd_dirty_pr,

    SRT_Boolean price_floor, double *floor_gearing, double *pre_factor_fwd,
    double *floor_strike, double *floor_implied_vol,

    SrtUndPtr *new_und, double *dirty_price) {
  Err err = NULL;
  SrtUndPtr und;
  Date date_today;
  char *yield_curve_name;
  double real_df_fixing, real_df_spot;
  Date spot_date;
  double forward;
  double floor_price;

  und = lookup_und(vasicek_und_name);
  if (!und)
    return serror("can not find underlying");

  date_today = get_today_from_underlying(und);

  yield_curve_name = (String)malloc(strlen(yc_name) + 1);
  strcpy(yield_curve_name, yc_name);

  err = srt_f_destroy_und(vasicek_und_name);
  if (err)
    return err;

  vasicek_mean_rev_level_data[1][bond_index] = new_mean_rev_level;

  err =
      SrtInitIRUnd(vasicek_und_name, yield_curve_name, model, num_sig, 2,
                   sig_data, num_tau, 2, tau_data, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   vasicek_init_cond, num_bond, 2, vasicek_mean_rev_level_data);

  if (err)
    return err;

  (*new_und) = lookup_und(vasicek_und_name);
  if (!(*new_und))
    return serror("can not find underlying");

  err = srt_f_vasicek_dirty_pr((*new_und), num_cpns[bond_index], settlt_date,
                               pay_dates[bond_index], cpns[bond_index],
                               discount_fwd_dirty_pr, dirty_price);
  if (err)
    return err;

  if (price_floor == SRT_YES) {
    err = srt_f_vasicek_discount_factor(
        (long)(pay_dates[bond_index][num_cpns[bond_index] - 1]), (*new_und),
        &real_df_fixing);

    if (err)
      return err;

    spot_date = get_spotdate_from_underlying((*new_und));

    err = srt_f_vasicek_discount_factor(spot_date, (*new_und), &real_df_spot);

    if (err)
      return err;

    forward = pre_factor_fwd[bond_index] * (real_df_fixing / real_df_spot);
    floor_price = srt_f_optblksch(
        forward, floor_strike[bond_index], floor_implied_vol[bond_index],
        (pay_dates[bond_index][num_cpns[bond_index] - 1] - date_today) *
            YEARS_IN_DAY,
        1.0, SRT_PUT, PREMIUM);

    floor_price *= floor_gearing[bond_index];

    (*dirty_price) += floor_price;
  }

  srt_free(yield_curve_name);

  return NULL;
}

Err srt_f_vasicek_calibrate(double settlt_date, double vasicek_init_cond,
                            long num_bonds, long *num_cpns, double **pay_dates,
                            double **cpns,

                            double *mkt_dirty_prices,
                            SRT_Boolean discount_fwd_dirty_pr,

                            SRT_Boolean price_floor, double *floor_gearing,
                            double *pre_factor_fwd, double *floor_strike,
                            double *floor_implied_vol,

                            char **vasicek_und_name) {
  Err err = NULL;
  SrtUndPtr und, new_und;
  SrtMdlType mdl_type;
  SrtMdlDim mdl_dim;
  Date date_today;
  TermStruct *ts;
  SrtLst *ls;
  IrmTermStructVal *tsval;
  long num_sig, num_tau;
  double **sig_data, **tau_data, **vasicek_mean_rev_level_data;
  char *yc_name;

  long k, i, sig_index, tau_index;
  double nstop, a[3], b[3];

  /* get the und from the underlying name */
  und = lookup_und((*vasicek_und_name));
  if (!und)
    return serror("fatal: can not get underlying");

  /* get today from the underlying */
  date_today = get_today_from_underlying(und);

  yc_name = get_discname_from_underlying(und);

  /* get the underlying term structure */
  err = get_underlying_ts(und, &ts);
  if (err)
    return err;

  ls = ts->head;
  num_sig = num_tau = 0;
  while (ls != NULL) {
    tsval = (IrmTermStructVal *)ls->element->val.pval;
    switch (tsval->val_origin) {
    case SIGMA_DATE:
      num_sig += 1;
      break;
    case TAU_DATE:
      num_tau += 1;
      break;
    case BOTH_DATE:
      num_sig += 1;
      num_tau += 1;
      break;
    }

    ls = ls->next;
  }

  sig_data = dmatrix(0, 1, 0, num_sig - 1);
  tau_data = dmatrix(0, 1, 0, num_tau - 1);
  vasicek_mean_rev_level_data = dmatrix(0, 1, 0, num_bonds - 1);

  ls = ts->head;
  sig_index = tau_index = 0;
  while (ls != NULL) {
    tsval = (IrmTermStructVal *)ls->element->val.pval;
    switch (tsval->val_origin) {
    case SIGMA_DATE:
      sig_index += 1;
      sig_data[0][sig_index - 1] = tsval->date;
      sig_data[1][sig_index - 1] = find_sig(tsval->time, ts);
      break;
    case TAU_DATE:
      tau_index += 1;
      tau_data[0][tau_index - 1] = tsval->date;
      tau_data[1][tau_index - 1] = find_tau(tsval->time, ts);
      break;
    case BOTH_DATE:
      tau_index += 1;
      sig_index += 1;
      sig_data[0][sig_index - 1] = tau_data[0][tau_index - 1] = tsval->date;

      sig_data[1][sig_index - 1] = find_sig(tsval->time, ts);
      tau_data[1][tau_index - 1] = find_tau(tsval->time, ts);

      break;
    }

    ls = ls->next;
  }

  for (i = 0; i < num_bonds; i++) {
    vasicek_mean_rev_level_data[0][i] = pay_dates[i][num_cpns[i] - 1];
    vasicek_mean_rev_level_data[1][i] = vasicek_init_cond;
  }

  if (err)
    return err;

  err = get_underlying_mdltype(und, &mdl_type);
  if (err)
    return err;

  err = get_underlying_mdldim(und, &mdl_dim);
  if (err)
    return err;

  for (i = 0; i < num_bonds; i++) {
    nstop = 0.0;

    a[0] = vasicek_init_cond;
    if (i > 0)
      yc_name = get_discname_from_underlying(new_und);

    err = one_bond_calib_func("VASICEK", mdl_type, mdl_dim, (*vasicek_und_name),
                              yc_name, sig_data, num_sig, tau_data, num_tau,
                              vasicek_mean_rev_level_data, vasicek_init_cond,
                              settlt_date, num_bonds, num_cpns, pay_dates, cpns,
                              i, a[0], discount_fwd_dirty_pr,

                              price_floor, floor_gearing, pre_factor_fwd,
                              floor_strike, floor_implied_vol,

                              &new_und, &b[0]);

    if (err)
      return err;

    a[1] = vasicek_init_cond + 1.0 / 100;
    yc_name = get_discname_from_underlying(new_und);

    err = one_bond_calib_func("VASICEK", mdl_type, mdl_dim, (*vasicek_und_name),
                              yc_name, sig_data, sig_index, tau_data, tau_index,
                              vasicek_mean_rev_level_data, vasicek_init_cond,
                              settlt_date, num_bonds, num_cpns, pay_dates, cpns,
                              i, a[1], discount_fwd_dirty_pr,

                              price_floor, floor_gearing, pre_factor_fwd,
                              floor_strike, floor_implied_vol,

                              &new_und, &b[1]);
    if (err)
      return err;

    a[2] = a[1] + 1.0 / 100;
    yc_name = get_discname_from_underlying(new_und);

    err = one_bond_calib_func("VASICEK", mdl_type, mdl_dim, (*vasicek_und_name),
                              yc_name, sig_data, sig_index, tau_data, tau_index,
                              vasicek_mean_rev_level_data, vasicek_init_cond,
                              settlt_date, num_bonds, num_cpns, pay_dates, cpns,
                              i, a[2], discount_fwd_dirty_pr,

                              price_floor, floor_gearing, pre_factor_fwd,
                              floor_strike, floor_implied_vol,

                              &new_und, &b[2]);
    if (err)
      return err;

    nstop = 0.0;
    k = 0;
    while (nstop < 1 && k < MAX_ITER) {
      newton(mkt_dirty_prices[i], 2.0, a, b, &nstop);

      yc_name = get_discname_from_underlying(new_und);
      err = one_bond_calib_func(
          "VASICEK", mdl_type, mdl_dim, (*vasicek_und_name), yc_name, sig_data,
          sig_index, tau_data, tau_index, vasicek_mean_rev_level_data,
          vasicek_init_cond, settlt_date, num_bonds, num_cpns, pay_dates, cpns,
          i, a[2], discount_fwd_dirty_pr,

          price_floor, floor_gearing, pre_factor_fwd, floor_strike,
          floor_implied_vol,

          &new_und, &b[2]);

      if (err)
        return err;
      k++;
    }
  }

  if (sig_data)
    free_dmatrix(sig_data, 0, 1, 0, num_sig - 1);
  sig_data = NULL;
  if (tau_data)
    free_dmatrix(tau_data, 0, 1, 0, num_tau - 1);
  tau_data = NULL;
  if (vasicek_mean_rev_level_data)
    free_dmatrix(vasicek_mean_rev_level_data, 0, 1, 0, num_bonds - 1);
  vasicek_mean_rev_level_data = NULL;

  return NULL;
}

Err srt_f_vasicek_calibrated_mean_rev_lvl(
    double settlt_date, double vasicek_init_cond, long num_bonds,
    long *num_cpns, double **pay_dates, double **cpns,

    double *bond_dirty_prices, SRT_Boolean discount_fwd_dirty_pr,

    SRT_Boolean price_floor, double *floor_gearing, double *pre_factor_fwd,
    double *floor_strike, double *floor_implied_vol,

    char *vasicek_und_name, double ***calibrated_mean_rev_lvl) {
  Err err = NULL;
  SrtUndPtr und;
  TermStruct *ts;
  long i;
  double date_today;

  err = srt_f_vasicek_calibrate(settlt_date, vasicek_init_cond, num_bonds,
                                num_cpns, pay_dates, cpns, bond_dirty_prices,
                                discount_fwd_dirty_pr,

                                price_floor, floor_gearing, pre_factor_fwd,
                                floor_strike, floor_implied_vol,

                                &vasicek_und_name);
  if (err)
    return err;

  und = lookup_und(vasicek_und_name);
  if (!und)
    return serror("can not find underlying");

  date_today = get_today_from_underlying(und);

  err = get_underlying_ts(und, &ts);
  if (err)
    return err;

  for (i = 0; i < num_bonds; i++) {
    (*calibrated_mean_rev_lvl)[i][0] = pay_dates[i][num_cpns[i] - 1];
    (*calibrated_mean_rev_lvl)[i][1] = find_mean_rev_level(
        (pay_dates[i][num_cpns[i] - 1] - date_today) * YEARS_IN_DAY, ts);
  }

  return NULL;
}

static struct calibration_struct {

  Date today;
  Date spot_date;

  int num_tau;
  int tau_cols;
  double **tau_data;

  int num_sig;
  int sig_cols;
  double **sig_data;

  double vasicek_init_cond;

  int num_levels;
  int levels_cols;
  double **mean_rev_levels_data;

  SrtMdlType mdl_type;
  SrtMdlDim mdl_dim;

  char *und_name;
};

static struct calibration_struct calibration_parms;

static Err init_static_params_struct(Date *overnight_swap_end_dates,
                                     long num_overnight_swap_rates,
                                     double *market_overnight_swap_rates,
                                     double *overnight_swap_rates_weights,
                                     char *und_name) {
  Err err = NULL;

  Date today, spot_date;
  SrtUndPtr und;
  TermStruct *ts;
  SrtLst *ls;
  IrmTermStructVal *tsval;
  long num_sig, num_tau;
  double **sig_data, **tau_data, **mean_rev_levels_data;
  double vasicek_init_cond;
  long sig_index, tau_index;
  long i;
  SrtMdlType mdl_type;
  SrtMdlDim mdl_dim;

  /* get the sort underlying pointer */
  und = lookup_und(und_name);
  if (!und)
    return serror("can not find underlying");

  /* get the date today from the sort undelying pointer */
  today = get_today_from_underlying(und);
  spot_date = get_spotdate_from_underlying(und);

  calibration_parms.today = today;
  calibration_parms.spot_date = spot_date;

  mean_rev_levels_data = dmatrix(0, 1, 0, num_overnight_swap_rates - 1);

  /* get the underlying term structure */
  err = get_underlying_ts(und, &ts);
  if (err)
    return err;

  ls = ts->head;
  num_sig = num_tau = 0;
  while (ls != NULL) {
    tsval = (IrmTermStructVal *)ls->element->val.pval;
    switch (tsval->val_origin) {
    case SIGMA_DATE:
      num_sig += 1;
      break;
    case TAU_DATE:
      num_tau += 1;
      break;
    case BOTH_DATE:
      num_sig += 1;
      num_tau += 1;
      break;
    }

    ls = ls->next;
  }

  calibration_parms.num_tau = num_tau;
  calibration_parms.num_sig = num_sig;
  calibration_parms.levels_cols = calibration_parms.sig_cols =
      calibration_parms.tau_cols = 2;

  sig_data = dmatrix(0, 1, 0, num_sig - 1);
  tau_data = dmatrix(0, 1, 0, num_tau - 1);

  ls = ts->head;
  sig_index = tau_index = 0;
  while (ls != NULL) {
    tsval = (IrmTermStructVal *)ls->element->val.pval;
    switch (tsval->val_origin) {
    case SIGMA_DATE:
      sig_index += 1;
      sig_data[0][sig_index - 1] = tsval->date;
      sig_data[1][sig_index - 1] = find_sig(tsval->time, ts);
      break;
    case TAU_DATE:
      tau_index += 1;
      tau_data[0][tau_index - 1] = tsval->date;
      tau_data[1][tau_index - 1] = find_tau(tsval->time, ts);
      break;
    case BOTH_DATE:
      tau_index += 1;
      sig_index += 1;
      sig_data[0][sig_index - 1] = tau_data[0][tau_index - 1] = tsval->date;

      sig_data[1][sig_index - 1] = find_sig(tsval->time, ts);
      tau_data[1][tau_index - 1] = find_tau(tsval->time, ts);

      break;
    }

    ls = ls->next;
  }

  calibration_parms.sig_data = sig_data;
  calibration_parms.tau_data = tau_data;

  err = srt_f_get_vasicek_init_cond(0.0, ts, &vasicek_init_cond);
  if (err)
    return err;

  calibration_parms.vasicek_init_cond = vasicek_init_cond;

  for (i = 0; i < num_overnight_swap_rates; i++) {
    mean_rev_levels_data[0][i] = overnight_swap_end_dates[i];
    mean_rev_levels_data[1][i] = vasicek_init_cond;
  }

  calibration_parms.mean_rev_levels_data = mean_rev_levels_data;

  err = get_underlying_mdltype(und, &mdl_type);
  if (err)
    return err;

  calibration_parms.mdl_type = mdl_type;

  err = get_underlying_mdldim(und, &mdl_dim);
  if (err)
    return err;

  calibration_parms.mdl_dim = mdl_dim;

  return NULL;
}

static Err
delete_static_params_struct(struct calibration_struct *calibration_parms) {
  int num_levels, num_sig, num_tau;
  struct calibration_struct tmp_calibration_parms;

  tmp_calibration_parms = (*calibration_parms);
  num_levels = tmp_calibration_parms.num_levels;
  num_sig = tmp_calibration_parms.num_sig;
  num_tau = tmp_calibration_parms.num_tau;

  if (tmp_calibration_parms.mean_rev_levels_data) {
    free_dmatrix(tmp_calibration_parms.mean_rev_levels_data, 0, 1, 0,
                 num_levels - 1);
    tmp_calibration_parms.mean_rev_levels_data = NULL;
  }

  if (tmp_calibration_parms.sig_data) {
    free_dmatrix(tmp_calibration_parms.sig_data, 0, 1, 0, num_sig - 1);
    tmp_calibration_parms.sig_data = NULL;
  }

  if (tmp_calibration_parms.tau_data) {
    free_dmatrix(tmp_calibration_parms.tau_data, 0, 1, 0, num_tau - 1);
    tmp_calibration_parms.tau_data = NULL;
  }

  return NULL;
}

static Err from_optparam_to_calibration_params(
    double *optim_params, long num_parms,
    struct calibration_struct *calibration_parms) {
  long i;

  for (i = 1; i <= num_parms; i++) {
    if (i == 1) {
      (*calibration_parms).vasicek_init_cond = optim_params[i];

    } else {
      (*calibration_parms).mean_rev_levels_data[1][i] = optim_params[i];
    }
  }

  return NULL;
}

static Err
levenberg_price_funcs(double data_index,
                      double optim_params[], /* from [1] to [num_parms] */
                      double *value, int num_params) {
  Err err;
  TermStruct *ts;
  SrtUndPtr und;

  err = from_optparam_to_calibration_params(optim_params, num_params,
                                            &calibration_parms);
  if (err)
    return err;

  /* initialise the term structure */

  err = srt_f_init_IRM_TermStruct(
      calibration_parms.today, calibration_parms.sig_data,
      calibration_parms.sig_cols, calibration_parms.num_sig,
      calibration_parms.tau_data, calibration_parms.tau_cols,
      calibration_parms.num_tau,

      calibration_parms.mdl_type, calibration_parms.mdl_dim,

      0.0, /*beta*/

      0.0, /*alpha*/
      0.0, /*gamma*/
      0.0, /*rho*/

      0.0, /*vovol*/

      0.0, /*etabeta*/

      calibration_parms.vasicek_init_cond, calibration_parms.num_levels,
      calibration_parms.levels_cols, calibration_parms.mean_rev_levels_data,

      &ts);

  if (err)
    return err;

  /* attaches the term structure to the underlying */

  und = lookup_und(calibration_parms.und_name);
  if (und == NULL)
    return serror("can not find underlying %s", calibration_parms.und_name);
  set_irund_ts(und, ts);

  err = srt_f_vasicek_overnight_interest_rate_swap(
      calibration_parms.und_name, calibration_parms.spot_date,
      (long)(calibration_parms.mean_rev_levels_data[0][(long)(data_index - 1)]),
      value);
  if (err)
    return err;

  err = free_underlying_ts(und);
  if (err)
    return err;

  return NULL;
}

static Err levenberg_calib_funcs(double instr_index, double optim_parms[],
                                 double *value, double deriv[],
                                 int num_params) {
  long i;
  double shift;
  Err err;

  /* computes the overnight swa rate */
  err = levenberg_price_funcs(instr_index, optim_parms, value, num_params);
  if (err)
    return err;
  /* computes the derivatives of its value with respect to each parameter */
  for (i = 1; i <= num_params; i++) {
    if (optim_parms[i] == 0.0)
      shift = ZERO_SHIFT;
    else
      shift = PROP_SHIFT * optim_parms[i];

    optim_parms[i] += shift;

    err = levenberg_price_funcs(instr_index, optim_parms, &(deriv[i]),
                                num_params);
    if (err)
      return err;

    deriv[i] -= *value;
    deriv[i] /= shift;

    /* resets parameters */

    optim_parms[i] -= shift;
  }

  return NULL;
}

Err srt_f_vasicek_overnights_calibrate(long num_overnight_swap_rates,
                                       Date *overnights_swap_end_dates,
                                       double *market_overnight_swap_rates,
                                       double *overnight_swap_rates_weights,
                                       char **und_name) {
  Err err = NULL;
  long i;
  double *data, *params;
  int niter = 10;
  double chisq;

  err = init_static_params_struct(
      overnights_swap_end_dates, num_overnight_swap_rates,
      market_overnight_swap_rates, overnight_swap_rates_weights, *und_name);
  if (err)
    return err;

  data = dvector(1, num_overnight_swap_rates);

  for (i = 1; i <= num_overnight_swap_rates; i++)
    data[i] = (double)(i);

  params = dvector(1, (num_overnight_swap_rates + 1));

  for (i = 2; i <= (num_overnight_swap_rates + 1); i++)
    params[i] = calibration_parms.vasicek_init_cond;

  err = levenberg_marquardt(
      data, market_overnight_swap_rates, overnight_swap_rates_weights,
      num_overnight_swap_rates, params, (num_overnight_swap_rates + 1), niter,
      levenberg_calib_funcs, &chisq);
  if (err)
    return err;

  /* get the calibrated underlying name */

  err = delete_static_params_struct(&calibration_parms);
  if (err)
    return err;

  if (data)
    free_dvector(data, 1, num_overnight_swap_rates);
  data = NULL;
  if (params)
    free_dvector(params, 1, num_overnight_swap_rates + 1);
  params = NULL;

  return NULL;
}

Err srt_f_vasicek_money_market_account_vol(Date date, char *vasicek_und_name,
                                           double *money_market_account_vol) {
  Err err = NULL;
  double date_today;
  SrtUndPtr und;
  TermStruct *ts;
  double var_int_sr, date_exp;

  und = lookup_und(vasicek_und_name);
  if (!und)
    return serror("fatal: can not get the underlying in "
                  "srt_f_vasicek_money_market_account_vol");

  date_today = get_today_from_underlying(und);

  err = get_underlying_ts(und, &ts);
  if (err)
    return err;

  date_exp = (date - date_today) * YEARS_IN_DAY;

  err = srt_f_get_vasicek_var_int_sr(date_exp, ts, &var_int_sr);
  if (err)
    return err;

  var_int_sr /= date_exp;

  (*money_market_account_vol) = sqrt(var_int_sr);

  return NULL;
}

Err srt_f_inflation_fwd(long start_date, long end_date,
                        char *price_index_und_name, double *inflation_fwd) {
  Err err = NULL;
  SrtUndPtr price_index_und, real_rate_und, nom_rate_und;
  char *real_rate_und_name, *nom_rate_und_name, *nom_rate_yc;
  TermStruct *nom_rate_ts, *real_rate_ts, *price_index_ts;
  double date_today, start_date_time, end_date_time;
  double df_nom_end_date, df_nom_start_date, df_real_end_date,
      df_real_start_date;
  double adjustment, G_real_start, H_real_start, Psi_real_start_end,
      Psi_nom_start_end;

  /* get the price index sort pointer */
  price_index_und = lookup_und(price_index_und_name);
  if (!price_index_und) {
    return serror("fatal: can find price index und.");
  }

  /* get the real rate underlying name & srt pointer */
  real_rate_und_name = get_forname_from_fxund(price_index_und);
  real_rate_und = lookup_und(real_rate_und_name);
  if (!real_rate_und) {
    return serror("fatal: can find real rate und.");
  }

  /* get the nominal rate underlying */
  nom_rate_und_name = get_domname_from_fxund(price_index_und);
  nom_rate_und = lookup_und(nom_rate_und_name);
  if (!nom_rate_und) {
    return serror("can not get the nom. rate underlying");
  }
  /* get the nominal rate yield curve */
  err = get_underlying_discname(nom_rate_und, &nom_rate_yc);
  if (err) {
    return err;
  }

  date_today = get_today_from_underlying(price_index_und);

  df_nom_end_date = swp_f_df(date_today, end_date, nom_rate_yc);

  df_nom_start_date = swp_f_df(date_today, start_date, nom_rate_yc);

  err =
      srt_f_vasicek_discount_factor(end_date, real_rate_und, &df_real_end_date);
  if (err)
    return err;

  err = srt_f_vasicek_discount_factor(start_date, real_rate_und,
                                      &df_real_start_date);
  if (err)
    return err;

  (*inflation_fwd) = 0.0;

  /* get the underlyings term structure */
  err = get_underlying_ts(real_rate_und, &real_rate_ts);
  if (err)
    return err;

  err = get_underlying_ts(nom_rate_und, &nom_rate_ts);
  if (err)
    return err;

  err = get_underlying_ts(price_index_und, &price_index_ts);
  if (err)
    return err;

  /* compute the different adjustment */
  start_date_time = (start_date - date_today) * YEARS_IN_DAY;
  end_date_time = (end_date - date_today) * YEARS_IN_DAY;

  adjustment = 0.0;

  G_H_func(start_date_time, real_rate_ts, &G_real_start, &H_real_start);

  Psi_real_start_end = Psi_func(end_date_time, real_rate_ts) -
                       Psi_func(start_date_time, real_rate_ts);
  Psi_nom_start_end = Psi_func(end_date_time, nom_rate_ts) -
                      Psi_func(start_date_time, nom_rate_ts);

  adjustment += Psi_real_start_end *
                (Psi_func(start_date_time, real_rate_ts) * G_real_start -
                 L_func(start_date_time, real_rate_ts) -
                 Psi_func(end_date_time, real_rate_ts) * G_real_start);
  adjustment += Psi_real_start_end * Psi_real_start_end * G_real_start;
  adjustment += Psi_real_start_end * M_fx_func(start_date_time, price_index_ts);

  adjustment +=
      Psi_real_start_end * (Psi_func(end_date_time, nom_rate_ts) *
                                O_fd_func(start_date_time, price_index_ts) -
                            R_fd_func(start_date_time, price_index_ts));
  adjustment -= Psi_real_start_end * Psi_nom_start_end *
                O_fd_func(start_date_time, price_index_ts);

  /*compute the fwd inflation */

  (*inflation_fwd) = (df_nom_start_date / df_nom_end_date) /
                     (df_real_start_date / df_real_end_date) * exp(adjustment);

  (*inflation_fwd) -= 1.0;

  return err;
}

/* ------------------------------------------------------------------------- */
double **find_M_eta_beta(double time, TermStruct *l) {
  SrtLst *ls;
  double **M;

  ls = l->head;

  /*
          M matrix is identical at all nodes of ts so just point
          to the first node
  */

  if (ls != NULL) {
    M = ((IrmTermStructVal *)ls->element->val.pval)->M_beta_eta;
  }
  return (M);
}
/* ------------------------------------------------------------------------- */

double Lambda_func(double starttime, double endtime, TermStruct *l) {

  return (Psi_func(endtime, l) - Psi_func(starttime, l)) / F_func(starttime, l);
}

/* ------------------------------------------------------------------------- */
/* Zeta corresponds to the variance of the state variable in the EtaBeta
   Exactly        , as the volatility of the state variable depends on the
   mean-reversion        , Zeta is the integral of sig^2 * exp(+2lambda*s) ds
   i.e. it is the G function
 */
double Zeta_func(double time, TermStruct *l) {
  double zeta;
  double rubbish;

  G_H_func(time, l, &zeta, &rubbish);

  return zeta;
  /*SrtLst          * ls;
  IrmTermStructVal   * tsval        , *tsval_p;
  double          zeta;

          ls = l->head;

          while ( (ls!=NULL) &&
          ( ((IrmTermStructVal*)ls->element->val.pval)->time < time) ) ls =
  ls->next; if ( ls == NULL ) ls = l->tail;

          tsval = (IrmTermStructVal*)ls->element->val.pval;

          if (ls != l->head)
          {
                  tsval_p = (IrmTermStructVal*)ls->previous->element->val.pval;
                  time -= tsval_p->time;
          }

          zeta = tsval->Zeta +  tsval->sig * tsval->sig * time;
          return ( zeta);
  */
}

/* ------------------------------------------------------------------------- */

double M_eta_beta_func(double s, double theta, double dOmega, double **M) {

  double M_stheta;
  double *s_arr, *theta_arr;
  int i;

  if (M == NULL)
    return 0.;

  M_stheta = 0.;

#if 1 /* use numeric M's for the moment */

  /* For 0        ,.5        ,1 used closed-form soln's */
  if (dOmega == 0.) {
    M_stheta = .5 * s * theta * theta;
    return (M_stheta);
  }
  if (dOmega == .5)
    if (theta == 0) {
      M_stheta = 0.;
      return (M_stheta);
    } else {
      M_stheta = (.5 * s * theta * theta) / (1. + .5 * s * theta);
      return (M_stheta);
    }
  if (dOmega == 1.) {
    M_stheta = .5 * theta * theta * (exp(s) - 1.) *
               (1. - .5 * theta / 3. * ((exp(s) - 1.) * (exp(s) + 2.)) +
                .25 * theta * theta / 12. * (exp(s) - 1.) * (exp(s) - 1.) *
                    (exp(3. * s) + 3. * exp(2. * s) + 6. * exp(s) + 6.));
    return (M_stheta);
  }

#endif

  /*     Interpolate M      */

  /* allocate memory Probably should think of a better way
     allocating and deallocating on each call will slow things down
     Most likely would like to carry along s        ,theta arrays which are
     precomputed        , stored on the termstruct ts        , and accessed with
     a func
  */
  s_arr = dvector(0, SMAX - 1);
  theta_arr = dvector(0, THETAMAX - 1);

  /*
     find nodal values for s        , theta. Currently hardwired so 0<s<1
     delta_s = .01 and 0<theta<10        , delta_theta = .1

  */
  i = 0;
  while (i < SMAX) {
    s_arr[i] = (double)i / 100.;
    i++;
  }

  i = 0;
  while (i < THETAMAX) {
    theta_arr[i] = (double)i / 10.;
    ;
    i++;
  }

  /*	M_stheta =  lq_inter_2d (s        , theta        , s_arr        ,
   * theta_arr        , M        , SMAX , THETAMAX);
   */
  /* deallocate memory */
  free_dvector(s_arr, 0, SMAX);
  free_dvector(theta_arr, 0, THETAMAX);

  return (M_stheta);
}

/* ------------------------------------------------------------------------- */

double find_beta(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->beta;
  } else {
    return (((IrmTermStructVal *)l->tail->element->val.pval)->beta);
  }
}

/* ------------------------------------------------------------------------- */

double find_eta(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->eta;
  } else {
    return (((IrmTermStructVal *)l->tail->element->val.pval)->eta);
  }
}

/* ------------------------------------------------------------------------- */

double find_omega(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->omega;
  } else {
    return (((IrmTermStructVal *)l->tail->element->val.pval)->omega);
  }
}

/* ------------------------------------------------------------------------- */
double find_tau(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->tau;
  } else {
    return (((IrmTermStructVal *)l->tail->element->val.pval)->tau);
  }
}

/* ------------------------------------------------------------------------- */

double find_sig(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->sig;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->sig;
  }
}

/* ----------------------------------------------------------------------------
 */

double find_vovol(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->vovol;
  } else {
    return (((IrmTermStructVal *)l->tail->element->val.pval)->vovol);
  }
}

/* ----------------------------------------------------------------------------
 */

double find_rho(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->rho;
  } else {
    return (((IrmTermStructVal *)l->tail->element->val.pval)->rho);
  }
}

double find_meanvol(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->meanvol;
  } else {
    return (((IrmTermStructVal *)l->tail->element->val.pval)->meanvol);
  }
}

double find_vasicek_init_cond(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->vasicek_init_cond;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->vasicek_init_cond;
  }
}

double find_mean_rev_level(double time, TermStruct *l) {
  SrtLst *ls;
  IrmTermStructVal *tsval;

  ls = l->head;
  tsval = (IrmTermStructVal *)ls->element->val.pval;

  while ((ls != NULL) && tsval->time < time) {
    ls = ls->next;
    tsval = (IrmTermStructVal *)ls->element->val.pval;
  }
  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->mean_rev_level;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->mean_rev_level;
  }
}

/* --------------- INTEPOLATION OF TERM STRUCTURE FUNCTIONS
 * -------------------*/
/* ------------------------------------------------------------------------- */

/* F_func def: exp(-int(s=0        ,s=t) l*s) */
double F_func(double time, TermStruct *l) {
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;

  /*	if (l == NULL) ENTER_DEBUG;  [MN] Not needed */

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = l->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  return (tsval->F * exp(-time / tsval->tau));
}

/* Psi_func def: int (s = 0        , s = t) [F(s) ds] */
double Psi_func(double time, TermStruct *l) {
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;

  /*	if (l == NULL) ENTER_DEBUG;  [MN] Not needed */

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = l->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }
  return (tsval->Psi + tsval->F * tsval->tau * (1 - exp((-time) / tsval->tau)));
}

/* ------------------------------------------------------------------------- */

double dln_F_func_dt(double time, TermStruct *l) {
  double val;

  /*	if (l==NULL) ENTER_DEBUG; [MN] Not needed */

  val = find_tau(time, l);
  return (-1.0 / val);
}

/* ------------------------------------------------------------------------- */

/* Descript: interpolates on (sig^2)*t        , computing the variance from the
                required dates up to the relevant IrmTermStructVal        , and
   then computing the variance from the first IrmTermStructVal to the second one
*/
double find_sig2_interp(double time1, double time2, TermStruct *l) {
  SrtLst *ls1, *ls2;
  IrmTermStructVal *tsval1, *tsval2;
  double var, var2, var1, sig1, sig2, prev_time;
  int ind; /* to index the cells of the 2linked ts */

  /*	if (l == NULL) ENTER_DEBUG;  [MN] Not needed */

  ls1 = l->head;
  ind = 0;
  var1 = 0.0;
  var2 = 0.0;

  /* If only 1 element in the ts return sig[0] */
  if (l->head == l->tail) {
    sig1 = ((IrmTermStructVal *)ls1->element->val.pval)->sig;
    return (sig1 * sig1);
  }
  /* Gets the closest (relevant) IrmTermStructVal that has a date bigger than t1
   */
  while ((ls1 != NULL) &&
         (((IrmTermStructVal *)ls1->element->val.pval)->time < time1)) {
    ls1 = ls1->next;
  }
  if (ls1 == NULL)
    ls1 = l->tail;

  ls2 = ls1;
  /* Gets the closest (relevant) IrmTermStructVal that has a date bigger than t2
   */
  while ((ls2 != NULL) &&
         (((IrmTermStructVal *)ls2->element->val.pval)->time < time2)) {
    /* ind represents the number of IrmTermStructVal that separates time1 and
     * 2*/
    ind++;
    ls2 = ls2->next;
  }
  /* If ls2 is NULL        , we are one IrmTermStructVal too far */
  if (ls2 == NULL) {
    ls2 = l->tail;
    ind--;
  }

  /* if time1 and time lead to the same time        , we return it */
  if (ind == 0) {
    sig1 = ((IrmTermStructVal *)ls1->element->val.pval)->sig;
    return (sig1 * sig1);
  } else {
    /* Gets the variance from time2 to the closest IrmTermStructVal */
    tsval2 = (IrmTermStructVal *)ls2->element->val.pval;
    sig2 = tsval2->sig;
    var2 = sig2 * sig2 * (tsval2->time - time2);
    var1 = 0.0;

    /* Gets the variance from time1 to the closest IrmTermStructVal */
    tsval1 = (IrmTermStructVal *)ls1->element->val.pval;
    sig1 = tsval1->sig;
    var1 = sig1 * sig1 * (tsval1->time - time1);

    prev_time = tsval1->time;
    /* Loops on all the IrmTermStructVal between 1 and 2 (thanks to ind) */
    while (ind > 0) {
      ls1 = ls1->next;
      tsval1 = (IrmTermStructVal *)ls1->element->val.pval;
      sig1 = tsval1->sig;
      var1 += sig1 * sig1 * (tsval1->time - prev_time);
      prev_time = tsval1->time;
      ind--;
    }
    /* Retuns the result: variance 1 minus variance 2 */
    var = var1 - var2;
    return var / (time2 - time1);
  }
}

/* G_func def: int(s = 0        , s = t) [s(s)*s(s)/(F(s)*F(s)] ds */
/* H_func def: int(s = 0        , s = t) [-s_ifr(s        ,t)*s_bond(s ,t)] ds =
 * F(t)*I(t)
 */

void G_H_func(double t, TermStruct *l, double *val_G, double *val_H) {
  SrtLst *ls;
  double temp, temp2, temp3;
  IrmTermStructVal *tsval, *tsval_p; /* _p for previous */

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < t))
    ls = ls->next;
  if (ls == NULL)
    ls = l->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) /*Not first element in the ts */
  {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval; /* [i-1] */
    t -= tsval_p->time;
  }
  temp = exp(-1 / tsval->tau * t);
  if (fabs(1 / tsval->tau) > EPS) {
    temp2 = (1 - temp) * tsval->tau;
    temp3 = (1.0 / (temp * temp) - 1.0) * tsval->tau;
  } else {
    temp2 = t;
    temp3 = 2.0 * temp2;
  }

  *val_G = tsval->G +
           0.5 * temp3 * (tsval->sig * tsval->sig) / (tsval->F * tsval->F);

  *val_H = temp * tsval->H + temp2 * temp * tsval->F * tsval->F * tsval->G +
           0.5 * tsval->sig * tsval->sig * temp2 * temp2;
}

/* J_func def: int(s = 0        , s = t) [s(s)/F(s)] ds */
double J_func(double t, TermStruct *l) {
  SrtLst *ls;
  double temp, temp2, val_J;
  IrmTermStructVal *tsval, *tsval_p; /* _p for previous */

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < t))
    ls = ls->next;
  if (ls == NULL)
    ls = l->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) /*Not first element in the ts */
  {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval; /* [i-1] */
    t -= tsval_p->time;
  }
  temp = exp(1 / tsval->tau * t);
  if (fabs(1 / tsval->tau) > EPS)
    temp2 = (temp - 1) * tsval->tau;
  else
    temp2 = t;

  val_J = tsval->J + temp2 * tsval->sig / tsval->F;

  return (val_J);
}

/* ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   Computes the forward LGM cumulative volatility between two dates t1 and t2 ,
   defined by:
                CumVol(t        , T1        , T2) = [Psi(T2)-Psi(T1)]^2 * G(t)
   ------------------------------------------------------------------------- */

double srt_f_lgm_cum_vol(TermStruct *ts, double fixing_time, double start_time,
                         double end_time, double *answer) {
  double G_t, Psi_t1, Psi_t2, tmp;
  double cv;

  Psi_t1 = Psi_func(start_time, ts);
  Psi_t2 = Psi_func(end_time, ts);
  G_H_func(fixing_time, ts, &G_t, &tmp);
  cv = (Psi_t2 - Psi_t1);
  cv *= cv;
  cv *= G_t;

  if (answer)
    *answer = cv;

  return (cv);
}

/* I_func def: int(s = 0        , s = t) [G(s)*F(s)] ds */
double I_func(double time, TermStruct *l) {
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;
  double result;
  double temp, temp1;
  result = 0.0;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = l->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  temp = (1 - exp(-time / tsval->tau)) * tsval->tau;
  temp1 = exp(time / tsval->tau) + exp(-time / tsval->tau) - 2;

  result = tsval->I + tsval->G * tsval->F * temp +
           (tsval->sig * tsval->tau) * (tsval->sig * tsval->tau) * temp1 /
               (2 * tsval->F);

  return result;
}

/* K_func def: int(s = 0        , s = t) [F(s)*J(s)] ds */
double K_func(double time, TermStruct *l)

{
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;
  double result;
  double temp;
  result = 0.0;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = l->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  temp = (1 - exp(-time / tsval->tau)) * tsval->tau;
  result = tsval->K + tsval->F * tsval->J * temp +
           tsval->sig * tsval->tau * (time - temp);

  return result;
}

/* L_func def: int(s = 0        , s = t) [phi(s)/F(s)] ds */
double L_func(double time, TermStruct *l)

{
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;
  double result;
  result = 0.0;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = l->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  result = tsval->L +
           tsval->F * tsval->G * (1 - exp(-time / tsval->tau)) * tsval->tau +
           0.5 * pow(tsval->sig * tsval->tau, 2) *
               (exp(-time / tsval->tau) + exp(time / tsval->tau) - 2) /
               tsval->F;

  return result;
}

/* O_func def: int(s = 0        , s = t) [phi(s)*psi(s)/F(s)] ds */
double O_func(double time, TermStruct *l)

{
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;
  double result;
  result = 0.0;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = l->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  result = tsval->O +
           tsval->F * tsval->G * tsval->Psi * (1 - exp(-time / tsval->tau)) *
               tsval->tau +
           tsval->F * tsval->F * tsval->G * tsval->tau *
               ((1 - exp(-time / tsval->tau)) * tsval->tau -
                0.5 * (1 - exp(-2 * time / tsval->tau)) * tsval->tau) +
           0.5 * (tsval->Psi / tsval->F) * tsval->tau * tsval->sig *
               tsval->sig *
               ((exp(time / tsval->tau) - 1) * tsval->tau -
                (1 - exp(-time / tsval->tau)) * tsval->tau) +
           0.5 * pow(tsval->sig * tsval->tau, 2) *
               ((exp(time / tsval->tau) - 1) * tsval->tau -
                (1 - exp(-time / tsval->tau)) * tsval->tau) +
           0.5 * pow(tsval->sig * tsval->tau, 2) *
               (0.5 * (1 - exp(-2 * time / tsval->tau)) * tsval->tau - time);

  return result;
}

double Q_func(double time, TermStruct *l)

{
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;
  double result;
  double temp, temp1;

  result = 0.0;
  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = l->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  temp = (1.0 - exp(-time / tsval->tau)) * tsval->tau;
  temp1 = 0.5 * (1.0 - exp(-2 * time / tsval->tau)) * tsval->tau;

  result = tsval->Q + tsval->H * temp +
           tsval->F * tsval->F * tsval->G * tsval->tau * (temp - temp1) +
           0.5 * tsval->sig * tsval->sig * tsval->tau * tsval->tau *
               (time - 2 * temp + temp1);

  return result;
}

/* S_func def: int(s = 0        , s = t) [s(s)*s(s)*Psi(s)/(F(s)*F(s)] */
double S_func(double time, TermStruct *l)

{
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;
  double temp, temp1, result;

  result = 0.0;
  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = l->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  temp = (exp(time / tsval->tau) - 1) * tsval->tau;
  temp1 = 0.5 * (exp(2 * time / tsval->tau) - 1) * tsval->tau;

  result =
      tsval->S +
      (tsval->sig * tsval->sig / (tsval->F * tsval->F)) * tsval->Psi * temp1 +
      tsval->sig * tsval->sig * tsval->tau * (temp1 - temp) / (tsval->F);

  return result;
}

/* T_func def: int(s = 0        , s = t) [F(s)*S(s)ds] */
double T_func(double time, TermStruct *l) {
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;
  double temp, temp1, temp2, result;

  result = 0.0;
  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = l->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != l->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  temp = (1 - exp(-time / tsval->tau)) * tsval->tau;
  temp1 = (exp(time / tsval->tau) - 1) * tsval->tau;
  temp2 = 0.5 * (exp(2 * time / tsval->tau) - 1) * tsval->tau;

  result =
      tsval->T + tsval->F * tsval->S * temp +
      0.5 * tsval->sig * tsval->sig * tsval->Psi / tsval->F * tsval->tau *
          (temp1 - temp) +
      0.5 * tsval->sig * tsval->sig * tsval->tau * tsval->tau * (temp1 - temp) -
      tsval->sig * tsval->sig * tsval->tau * tsval->tau * (time - temp);

  return result;
}

/* -------------------------------------------------------------------------------
 */

/* For Quanto Adjustment in LGM Jumping : computes the M function */

Err srt_f_extend_lgm_jumping_ts_for_quanto(TermStruct *lgm_ts,
                                           TermStruct *fx_ts,
                                           char *szLgmUndName,
                                           char *szFxUndName) {
  SrtListAtom *lgm_lc, *lgm_lp, *fx_lc, *corr_lc;
  IrmTermStructVal *lgm_tsval, *lgm_tsval_p;
  FxTermStructVal *fx_tsval;
  SrtCorrLstVal *corrval;
  SrtCorrLstPtr corrlist;
  double dCorrelation;
  long lTicker;
  double dFxVolatility;
  double time;
  double temp;
  Err err;

  /* ------------------------------ Fx Dates ------------------------------ */

  /* Starts at the top of the FX Term Struct */
  fx_lc = fx_ts->head;

  /* Starts at the top of the LGM Term Struct */
  lgm_lc = lgm_ts->head;

  /* Loops on the Fx dates to insert them into the IRM Term Struct */
  while (fx_lc != NULL) {
    /* Gets the Term Struct Val attached */
    fx_tsval = (FxTermStructVal *)fx_lc->element->val.pval;

    /* Moves on to the IRM date just after (or on) this date */
    while ((lgm_lc != NULL) &&
           (((IrmTermStructVal *)lgm_lc->element->val.pval)->time <
            fx_tsval->time))
      lgm_lc = lgm_lc->next;
    if (lgm_lc == NULL)
      lgm_lc = lgm_ts->tail;
    lgm_tsval = ((IrmTermStructVal *)lgm_lc->element->val.pval);

    /* If not on the same date        , we add the date onto the LGM ts */
    if (lgm_tsval->time != fx_tsval->time) {
      /* Creation of a new IrmTermStructVal (simple element of the Term Struct)
       */
      lgm_tsval = (IrmTermStructVal *)srt_calloc(1, sizeof(IrmTermStructVal));

      /* Transfer of the information for the LGM term struct */
      lgm_tsval->time = fx_tsval->time;
      lgm_tsval->date = fx_tsval->date;
      lgm_tsval->val_origin = BOTH_DATE;
      lgm_tsval->sig = find_sig(fx_tsval->time, lgm_ts);
      lgm_tsval->F = find_F(fx_tsval->time, lgm_ts);
      lgm_tsval->tau = find_tau(fx_tsval->time, lgm_ts);
      lgm_tsval->beta = 0.0;
      lgm_tsval->vovol = 0.0;
      lgm_tsval->rho = 0.0;

      /* Insert the TermStructAtom in the TS (==linked ts sorted by key = date)
       */
      srt_f_lstins(lgm_ts, "OneFacTsAtom", lgm_tsval->date,
                   OBJ_PTR_IRM_TermStruct, (void *)lgm_tsval,
                   &srt_f_irmtsvalfree, &lTicker);

    } /* END if lgm time != fx time */

    /* Moves on to the next FX Term Struct date */
    fx_lc = fx_lc->next;

  } /* END of loop on FX Term Structure Dates */

  /* ------------------------------ Correlation Dates
   * ------------------------------ */

  /* Gets the Global Correlation List */
  corrlist = srt_f_GetTheCorrelationList();
  if (!corrlist->head->element)
    return serror("Correlation list improperly initialised...");

  /* Starts at the top of the Correlation Term Struct */
  corr_lc = corrlist->head;

  /* Starts at the top of the LGM Term Struct */
  lgm_lc = lgm_ts->head;

  /* Loops on the Correlation dates to insert them into the IRM Term Struct */
  while (corr_lc != NULL) {
    /* Gets the Term Struct Val attached */
    corrval = (SrtCorrLstVal *)corr_lc->element->val.pval;

    /* Moves on to the IRM date just after (or on) this date */
    while (
        (lgm_lc != NULL) &&
        (((IrmTermStructVal *)lgm_lc->element->val.pval)->time < corrval->time))
      lgm_lc = lgm_lc->next;
    if (lgm_lc == NULL)
      lgm_lc = lgm_ts->tail;
    lgm_tsval = ((IrmTermStructVal *)lgm_lc->element->val.pval);

    /* If not on the same date        , we add the date onto the LGM ts */
    if (lgm_tsval->time != corrval->time) {
      /* Creation of a new IrmTermStructVal (simple element of the Term Struct)
       */
      lgm_tsval = (IrmTermStructVal *)srt_calloc(1, sizeof(IrmTermStructVal));

      /* Transfer of the information for the LGM term struct */
      lgm_tsval->time = corrval->time;
      lgm_tsval->date = corrval->date;
      lgm_tsval->val_origin = BOTH_DATE;
      lgm_tsval->F = find_F(fx_tsval->time, lgm_ts);
      lgm_tsval->sig = find_sig(corrval->time, lgm_ts);
      lgm_tsval->tau = find_tau(corrval->time, lgm_ts);
      lgm_tsval->beta = 0.0;
      lgm_tsval->vovol = 0.0;
      lgm_tsval->rho = 0.0;

      /* Insert the TermStructAtom in the TS (==linked ts sorted by key = date)
       */
      srt_f_lstins(lgm_ts, "OneFacTsAtom", lgm_tsval->date,
                   OBJ_PTR_IRM_TermStruct, (void *)lgm_tsval,
                   &srt_f_irmtsvalfree, &lTicker);

    } /* END if lgm time != fx time */

    /* Moves on to the next Correlation Term Struct date */
    corr_lc = corr_lc->next;

  } /* END of loop on Correlation Term Structure Dates */

  /* --------------------- Stores Rho and M ----------------------------- */

  /* Starts at the top of the LGM Term Struct */
  lgm_lc = lgm_ts->head;

  /* Gets the Term Struct Val attached */
  lgm_tsval = ((IrmTermStructVal *)lgm_lc->element->val.pval);

  /* Gets the Correlation between the LGM und and the FX one */
  err = srt_f_get_corr_from_CorrList(corrlist, szLgmUndName, szFxUndName,
                                     lgm_tsval->time, &dCorrelation);
  if (err)
    return err;

  /* Gets the FX volatility on that date */
  dFxVolatility = find_fx_sig(lgm_tsval->time, fx_ts);

  /* Stores the correlation of the FX times its vol */
  lgm_tsval->dFxVolTimesCor = dCorrelation * dFxVolatility;

  /* The Value of the M function on the first point: 0.0 */
  lgm_tsval->M = 0.0;

  /* Now that all FX and Correlation dates have been inserted        , compute
   * the M function on the LGM Ts Dates */
  while (lgm_lc->next != NULL) {
    if (lgm_lc->previous != NULL) {
      lgm_lp = lgm_lc->previous;
      lgm_tsval_p = ((IrmTermStructVal *)lgm_lp->element->val.pval);
      time = lgm_tsval->time - lgm_tsval_p->time;
    } else
      time = lgm_tsval->time;

    lgm_tsval_p = lgm_tsval;
    lgm_lc = lgm_lc->next;
    lgm_tsval = ((IrmTermStructVal *)lgm_lc->element->val.pval);

    /* Gets the Correlation between the LGM und and the FX from the previous
     * date to the current one*/
    err = srt_f_get_corr_from_CorrList(corrlist, szLgmUndName, szFxUndName,
                                       lgm_tsval->time, &dCorrelation);
    if (err)
      return err;

    /* Gets the FX volatility from the previous date to the current one */
    dFxVolatility = find_fx_sig(lgm_tsval->time, fx_ts);

    /* Computes and stores the correlation */
    lgm_tsval->dFxVolTimesCor = dCorrelation * dFxVolatility;

    /* Computes and stores the M function (quanto adjustment term) (on the end
     * of the period) */
    if (fabs(1 / lgm_tsval_p->tau) > EPS)
      temp = lgm_tsval_p->tau * (exp(time / lgm_tsval_p->tau) - 1.0);
    else
      temp = time;

    lgm_tsval->M = lgm_tsval_p->M + lgm_tsval_p->dFxVolTimesCor *
                                        (lgm_tsval_p->sig / lgm_tsval_p->F) *
                                        temp;
  }

  /* Returns a success message */
  return NULL;

} /* END srt_f_extend_lgm_jumping_ts_for_quanto(...) */

/* ---------------------------------------------------------------------------------
 */

/* FOr Quanto Adjustments in LGM Jumping */

double M_func(double time, TermStruct *lgm_ts) {
  SrtLst *ls;
  IrmTermStructVal *tsval, *tsval_p;
  double temp;

  ls = lgm_ts->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;
  if (ls == NULL)
    ls = lgm_ts->tail;

  tsval = (IrmTermStructVal *)ls->element->val.pval;

  if (ls != lgm_ts->head) {
    tsval_p = (IrmTermStructVal *)ls->previous->element->val.pval;
    time -= tsval_p->time;
  }

  if (fabs(1 / tsval->tau) > EPS)
    temp = (exp(time / tsval->tau) - 1) * tsval->tau;
  else
    temp = time;

  return tsval->M + tsval->dFxVolTimesCor * (tsval->sig / tsval->F) * temp;
}

/* ---------------------------------------------------------------------------------
 */

double find_F(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->F;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->F;
  }
}

double find_G(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->G;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->G;
  }
}

double find_H(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->H;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->H;
  }
}

double find_Psi(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->Psi;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->Psi;
  }
}

double find_I(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->I;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->I;
  }
}

double find_J(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->J;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->J;
  }
}

double find_K(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->K;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->K;
  }
}

double find_L(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->L;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->L;
  }
}

double find_O(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->O;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->O;
  }
}

double find_Q(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->Q;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->Q;
  }
}

double find_Phi(double time, TermStruct *l) {
  SrtLst *ls;

  ls = l->head;

  while ((ls != NULL) &&
         (((IrmTermStructVal *)ls->element->val.pval)->time < time))
    ls = ls->next;

  if (ls != NULL) {
    return ((IrmTermStructVal *)ls->element->val.pval)->Phi;
  } else {
    return ((IrmTermStructVal *)l->tail->element->val.pval)->Phi;
  }
}

/* ------------------------------------------------------------------------- */
/* FOR THE NEW LGM MODEL WE NEED TO KNOW
/* the functions :
/*                 H(t) and G(t)  beware we have use the former label but there
/*                 are not the same. We lineary interpolate those function
between two dates
/*                 for t > tn   H(t) = Hn-1 + (Hn - Hn-1) (t - tn-1) / (tn -
tn-1)
/*                 for t < t1   H(t) = H1 + (H2 - H1) (t - t1) / (t2 - t1)
/*                 same thing for G(t)
/* ------------------------------------------------------------------------- */

static double get_struct(SrtLst *ls, LabelStruct ind) {
  if (ind == 0)
    return (((IrmTermStructVal *)ls->element->val.pval)->G);
  else if (ind == 1)
    return (((IrmTermStructVal *)ls->element->val.pval)->H);
  else
    return -1;
}

double find_struct_interp(double time, LabelStruct label, TermStruct *l) {
  SrtLst *ls1;
  double time1, time2;
  double G1, G2;

  /*	if (l == NULL) ENTER_DEBUG;  [MN] Not needed */

  ls1 = l->head;

  /* If only 1 element in the ts return G1 * t / t1 */
  if (l->head == l->tail) {
    G1 = get_struct(ls1, label);
    time1 = ((IrmTermStructVal *)ls1->element->val.pval)->time;
    return (G1 * time / time1);
  }
  /* Gets the closest (relevant) IrmTermStructVal that has a date bigger than t1
   */
  while ((ls1 != NULL) &&
         (((IrmTermStructVal *)ls1->element->val.pval)->time < time)) {
    ls1 = ls1->next;
  }

  if (ls1 == NULL) {
    ls1 = l->tail;
  }

  G2 = get_struct(ls1, label);
  time2 = ((IrmTermStructVal *)ls1->element->val.pval)->time;

  if (ls1->previous == NULL) {
    return G2 * time / time2;
  } else {
    ls1 = ls1->previous;
  }
  G1 = get_struct(ls1, label);
  time1 = ((IrmTermStructVal *)ls1->element->val.pval)->time;

  return G1 + (G2 - G1) * (time - time1) / (time2 - time1);
}

double find_struct_interp_der(double time, LabelStruct label, TermStruct *l) {
  double result;
  double shift = YEARS_IN_DAY;

  result = find_struct_interp(time + shift, label, l);
  if (time - shift > 0) {
    result -= find_struct_interp(time - shift, label, l);
    result /= 2 * shift;
  } else {
    result -= find_struct_interp(time, label, l);
    result /= shift;
  }

  return result;
}

#undef MAX_ITER
#undef ZERO_SHIFT
#undef PROP_SHIFT
