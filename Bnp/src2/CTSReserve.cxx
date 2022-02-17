
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"

#include "LGMSVClosedForm.h"
#include "LGMSVClosedFormApprox.h"
#include "LGMSVPDE.h"
#include "LGMSVUtil.h"

#include "AmortMidatCalib.h"
#include "CTSProdStruct.h"
#include "CTSReserve.h"

#define MAX_FIRST_VOL 0.10
#define MAX_FIRST_MAT 0.40

#define MAX_ONETIME 200

#define LGMSV_MINVOL 0.0001
#define LGMSV_MAXVOL 0.1
#define MAX_EURO 200
#define BUMP_VOL 0.0005

#define OPT_EPS 1.0e-15
#define VOL_MULT 5.0
#define VOL_MULT_ITER 5

/*	 Static variables for the Levenberg Marquard */
static CTS_UND stat_und;
static CTS stat_market_cts, stat_model_cts;

static int stat_for_fund, stat_calc_fwd_iv, stat_adj_fee, stat_nb_euro;

static CTS_ADI_ARG stat_adi_arg;

static LGMSV_NUMERPARAMS stat_num_params;

static double stat_price[MAX_EURO], *stat_target,
    stat_sensi[MAX_EURO][MAX_EURO];

static double stat_call, stat_numer_tstar, stat_euro[MAX_EURO];

Err static_cts_price_european(double index, double vol[], double *price,
                              double *gradient, int nvol) {
  static int i, j, ind;
  static double sum2;
  Err err = NULL;

  ind = (int)(index);

  if (ind == 0) {
    /* We calculate everything */

    /* Initialise the modell */

    for (i = 0; i < nvol; i++) {
      stat_und->model.dSigma[i] = vol[i + 1];

      if (stat_und->model.dSigma[i] < LGMSV_MINVOL) {
        stat_und->model.dSigma[i] = LGMSV_MINVOL;
      }
      if (stat_und->model.dSigma[i] > LGMSV_MAXVOL) {
        stat_und->model.dSigma[i] = LGMSV_MAXVOL;
      }
    }

    /* Calculate Sensi */
    for (j = 0; j < nvol; j++) {
      stat_und->model.dSigma[j] += BUMP_VOL;

      Convert_Tstar_model(&(stat_und->model), stat_numer_tstar);

      /* prices */
      err = cts_calc_model_fwd_iv(stat_und, stat_model_cts, stat_calc_fwd_iv,
                                  stat_adj_fee, stat_num_params);

      if (err)
        goto FREE_RETURN;

      err = cts_adjust_model_fwd_iv(stat_und, stat_market_cts, stat_model_cts,
                                    stat_for_fund, stat_calc_fwd_iv,
                                    stat_adj_fee, 0, 0, NULL);

      if (err)
        goto FREE_RETURN;

      err = cts_launch_algo(stat_model_cts, stat_und, stat_adi_arg, &stat_call,
                            &(stat_euro[0]), NULL);

      if (err)
        goto FREE_RETURN;

      for (i = 0; i < stat_nb_euro; i++) {
        stat_sensi[i][j] = stat_euro[i];
      }

      Convert_Tstar_model(&(stat_und->model), stat_und->model.dInitTStar);

      stat_und->model.dSigma[j] -= BUMP_VOL;
    }

    err = cts_calc_model_fwd_iv(stat_und, stat_model_cts, stat_calc_fwd_iv,
                                stat_adj_fee, stat_num_params);

    if (err)
      goto FREE_RETURN;

    err = cts_adjust_model_fwd_iv(stat_und, stat_market_cts, stat_model_cts,
                                  stat_for_fund, stat_calc_fwd_iv, stat_adj_fee,
                                  0, 0, NULL);

    if (err)
      goto FREE_RETURN;

    Convert_Tstar_model(&(stat_und->model), stat_numer_tstar);

    err = cts_launch_algo(stat_model_cts, stat_und, stat_adi_arg, &stat_call,
                          &(stat_euro[0]), NULL);

    Convert_Tstar_model(&(stat_und->model), stat_und->model.dInitTStar);

    if (err)
      goto FREE_RETURN;

    /* Calculate Prices */

    sum2 = 0.0;

    for (i = 0; i < stat_nb_euro; i++) {
      for (j = 0; j < nvol; j++) {
        stat_sensi[i][j] = (stat_sensi[i][j] - stat_euro[i]) / BUMP_VOL;
      }

      stat_price[i] = stat_euro[i] - stat_target[i];

      sum2 += stat_price[i] * stat_price[i];
    }

    sum2 = sqrt(sum2);
  }

  /* Return the asked result */

  (*price) = stat_price[ind];

  for (j = 0; j < nvol; j++) {
    gradient[j + 1] = stat_sensi[ind][j];
  }

FREE_RETURN:

  return err;
}

Err cts_calibrate_european_and_price(CTS_UND und, CTS market_cts, CTS model_cts,
                                     int for_fund, int calc_fwd_iv, int adj_fee,

                                     CTS_ADI_ARG adi_arg,
                                     LGMSV_NUMERPARAMS numer_params,
                                     double numer_tstar,

                                     int nb_european, double *european_values,

                                     int nb_iter,

                                     double *call, double *fitting_error) {
  Err err = NULL;

  double data[MAX_CPN], weight[MAX_CPN], target[MAX_CPN], vol[MAX_CPN];

  int i;

  for (i = 0; i < nb_european; i++) {
    data[i + 1] = i;
    weight[i + 1] = 1.0 / european_values[i];
    target[i + 1] = 0.0;
  }

  for (i = 0; i < und->model.iNbPWTime; i++) {
    vol[i + 1] = und->model.dSigma[i];
  }

  stat_und = und;
  stat_market_cts = market_cts;
  stat_model_cts = model_cts;
  stat_for_fund = for_fund;
  stat_calc_fwd_iv = calc_fwd_iv;
  stat_adj_fee = adj_fee;
  stat_adi_arg = adi_arg;
  stat_num_params = numer_params;
  stat_numer_tstar = numer_tstar;

  stat_nb_euro = nb_european;
  stat_target = european_values;

  err = levenberg_marquardt(data, target, weight, nb_european, vol,
                            und->model.iNbPWTime, nb_iter,
                            static_cts_price_european, fitting_error);

  *call = stat_call;

  return err;
}

static double optval(double fwd, double strike, double vol, double mat,
                     double disc_fact, SrtCallPutType call_put,
                     SrtGreekType greek, SrtDiffusionType log_or_norm) {
  double premium;

  if (log_or_norm == SRT_LOGNORMAL) {
    premium =
        srt_f_optblksch(fwd, strike, vol, mat, disc_fact, call_put, greek);
  } else {
    premium =
        srt_f_optblknrm(fwd, strike, vol, mat, disc_fact, call_put, greek);
  }

  return premium;
}

Err cts_optimpvol(double premium, double fwd_price, double strike, double mat,
                  double disc, SrtCallPutType call_put,
                  SrtDiffusionType lognormal_normal, double *implied_vol) {
  int i;
  double vol_up, vol_down, vol_middle, vol_temp, vol_acc, vol_shift;
  double vol_diffold, vol_diff;
  double store_fwd_price;
  double deriv;
  double prem_new, prem_up, prem_down;
  double intrinsic;

  /* value by default */
  *implied_vol = 0.0;

  /* Compute option intrinsic value */
  intrinsic = optval(fwd_price, strike, NULL_VOL, mat, disc, call_put, PREMIUM,
                     lognormal_normal);

  /* Checks target premium is above intrinsic value */
  if (intrinsic > premium + OPT_EPS) {
    return serror("Intrinsic higher than option premium");
  }

  /* Renormalise the option price by the spot value and factorise df */
  if (fwd_price != 0.00) {
    store_fwd_price = fabs(fwd_price);
    strike /= store_fwd_price;
    premium /= store_fwd_price;
    fwd_price = fwd_price / store_fwd_price;

    premium /= disc;
    disc = 1.0;
  }

  /* Sets initial guess for vol */
  vol_up = 20;
  vol_down = 0.000000001;
  vol_shift = 0.0000000001;
  vol_acc = 0.000001;

  prem_up = optval(fwd_price, strike, vol_up, mat, disc, call_put, PREMIUM,
                   lognormal_normal);

  prem_down = optval(fwd_price, strike, vol_down, mat, disc, call_put, PREMIUM,
                     lognormal_normal);

  if ((fabs(prem_up - premium) < DBL_EPSILON) &&
      (fabs(prem_down - premium) < DBL_EPSILON)) {
    return serror("Too many solution in range");
  }

  if (fabs(prem_down - premium) < DBL_EPSILON) {
    *implied_vol = vol_down;
    if (lognormal_normal == SRT_NORMAL)
      *implied_vol *= store_fwd_price;
    return NULL;
  }

  if (fabs(prem_up - premium) < DBL_EPSILON) {
    *implied_vol = vol_up;
    if (lognormal_normal == SRT_NORMAL)
      *implied_vol *= store_fwd_price;
    return NULL;
  }

  /* same sign for the interval */
  if ((prem_down - premium) * (prem_up - premium) > 0.0) {
    /* modification for very low forward in JPY */
    if (prem_up < premium) {
      /* one more try with vol_up * 5.0 */
      i = 0;
      while (i < VOL_MULT_ITER && prem_up < premium) {
        vol_up *= VOL_MULT;
        prem_up = optval(fwd_price, strike, vol_up, mat, disc, call_put,
                         PREMIUM, lognormal_normal);
        i++;
      }

      if (prem_up < premium) {
        return serror("No solution in range");
      }
    } else {
      return serror("No solution in range");
    }
  }

  /* orientation of the search */
  if (prem_up < 0) {
    vol_middle = vol_up;
    vol_up = vol_down;
    vol_down = vol_middle;
  }

  vol_middle = 0.5 * (vol_up + vol_down);
  vol_diffold = fabs(vol_up - vol_down);
  vol_diff = vol_diffold;

  prem_new = optval(fwd_price, strike, vol_middle, mat, disc, call_put, PREMIUM,
                    lognormal_normal);

  deriv = (optval(fwd_price, strike, vol_middle + vol_shift, mat, disc,
                  call_put, PREMIUM, lognormal_normal) -
           prem_new) /
          vol_shift;

  for (i = 0; i <= MAX_ITER; i++) {
    if ((((vol_middle - vol_up) * deriv - (prem_new - premium)) *
             ((vol_middle - vol_down) * deriv - (prem_new - premium)) >=
         0.0) ||
        (fabs(2.0 * (prem_new - premium)) > fabs(vol_diffold * deriv))) {
      /* bissection if Newton is out of range        , or not decreasing fast
       * enough
       */
      vol_diffold = vol_diff;
      vol_diff = 0.5 * (vol_up - vol_down);
      vol_middle = vol_down + vol_diff;
      if (vol_down == vol_middle) /* The change is negligible */
      {
        *implied_vol = vol_middle;
        if (lognormal_normal == SRT_NORMAL)
          *implied_vol *= store_fwd_price;
        return NULL;
      }
    } else {
      /* the change in newton is acceptable        , take it */
      vol_diffold = vol_diff;
      vol_diff = (prem_new - premium) / deriv;
      vol_temp = vol_middle;
      vol_middle -= vol_diff;
      if (vol_temp == vol_middle) {
        *implied_vol = vol_middle;
        if (lognormal_normal == SRT_NORMAL)
          *implied_vol *= store_fwd_price;
        return NULL;
      }
    }

    if (fabs(vol_diff) < vol_acc) {
      *implied_vol = vol_middle;
      if (lognormal_normal == SRT_NORMAL)
        *implied_vol *= store_fwd_price;
      return NULL;
    }

    prem_new = optval(fwd_price, strike, vol_middle, mat, disc, call_put,
                      PREMIUM, lognormal_normal);

    deriv = (optval(fwd_price, strike, vol_middle + vol_shift, mat, disc,
                    call_put, PREMIUM, lognormal_normal) -
             prem_new) /
            vol_shift;

    /* maintain the bracket on the root */
    if ((prem_new - premium) < 0.0) {
      vol_down = vol_middle;
    } else {
      vol_up = vol_middle;
    }
  }

  *implied_vol = 0.0;

  return NULL;

} /* srt_f_optimpvol() */

Err cts_get_one_time_infos(CTS_MKT mkt, CTS_ADI_ARG adi_arg_for_midat, CTS cts,
                           double vol_shift, double *onetime_call,
                           int *one_time_calc, double *european_prices,
                           double *european_prices_bumped,
                           double *european_fees, int *nb_one_time_calc,
                           double *most_expensive_pv, CTS_IV fwd_iv_info) {
  Err err = NULL;
  int i, j, k, last_k;
  int index_call, index_onetime;
  long today;
  CTS_PAY_ARG cts_prm;
  SrtCallPutType onetime_type;
  double onetime_fwd, onetime_strike, onetime_level, onetime_notional,
      onetime_vol, onetime_pv, onetime_mat;

  double normal_vol;

  double temp;

  today = mkt->today;

  /* Memory allocation */
  fwd_iv_info->has_one_time = 1;

  fwd_iv_info->one_time_call =
      (double *)calloc(fwd_iv_info->ncall, sizeof(double));
  fwd_iv_info->one_time_coupon =
      (double *)calloc(fwd_iv_info->ncall, sizeof(double));
  fwd_iv_info->one_time_log_vol =
      (double *)calloc(fwd_iv_info->ncall, sizeof(double));
  fwd_iv_info->one_time_norm_vol =
      (double *)calloc(fwd_iv_info->ncall, sizeof(double));
  fwd_iv_info->one_time_strike =
      (double *)calloc(fwd_iv_info->ncall, sizeof(double));
  fwd_iv_info->one_time_notional =
      (double *)calloc(fwd_iv_info->ncall, sizeof(double));

  if (!fwd_iv_info->one_time_call || !fwd_iv_info->one_time_coupon ||
      !fwd_iv_info->one_time_log_vol || !fwd_iv_info->one_time_norm_vol ||
      !fwd_iv_info->one_time_strike || !fwd_iv_info->one_time_notional) {
    err = "Memory allocation faillue in cts_get_one_time_infos";
    goto FREE_RETURN;
  }

  if (cts->call[0].pay_rec == 1) {
    onetime_type = SRT_PUT;
  } else {
    onetime_type = SRT_CALL;
  }

  for (i = 0; i < cts->exo_leg->num_cpn; i++) {
    one_time_calc[i] = 0;
  }

  index_onetime = 0;
  index_call = 0;
  *nb_one_time_calc = 0;

  *most_expensive_pv = 0.0;

  for (i = 0; i < adi_arg_for_midat->nstp; i++) {
    if (adi_arg_for_midat->is_event[i]) {
      cts_prm = (CTS_PAY_ARG)(adi_arg_for_midat->void_prm[i]);

      if (cts_prm->is_call) {
        fwd_iv_info->one_time_coupon[index_call] =
            cts->exo_leg->cpn[cts_prm->call->exo_idx].cxxpn_coupon;
        fwd_iv_info->one_time_notional[index_call] =
            cts->fund_leg->cpn[cts_prm->call->fund_idx].not ;

        if (cts_prm->eval_const_call.eval_one_time) {
          one_time_calc[cts_prm->call->exo_idx] = 1;
          european_prices[cts_prm->call->exo_idx] = onetime_call[index_onetime];
          european_fees[cts_prm->call->exo_idx] = cts->call[index_call].fee;

          fwd_iv_info->one_time_call[index_call] = onetime_call[index_onetime];

          if (onetime_call[index_onetime] > *most_expensive_pv) {
            *most_expensive_pv = onetime_call[index_onetime];
          }

          onetime_level = 0.0;
          onetime_strike = 0.0;
          onetime_notional = 0.0;
          onetime_fwd = 0.0;

          k = cts_prm->call->fund_idx;

          for (j = cts_prm->call->exo_idx; j < cts->exo_leg->num_cpn; j++) {
            temp =
                cts->exo_leg->cpn[j].cxxpn_cvg *
                swp_f_df(today, cts->exo_leg->cpn[j].cxxpn_pay_date, mkt->yc);

            last_k = k;

            while (k < cts->fund_leg->num_cpn - 1 &&
                   cts->fund_leg->cpn[k].end_date <
                       cts->exo_leg->cpn[j].cxxpn_end_date - 10) {
              k++;
            }

            onetime_level += temp * cts->fund_leg->cpn[k].not ;
            onetime_strike += temp * cts->exo_leg->cpn[j].cxxpn_coupon *
                              cts->exo_leg->cpn[j].cxxpn_not;
            onetime_fwd +=
                (swp_f_df(today, cts->fund_leg->cpn[last_k].start_date,
                          mkt->yc) -
                 swp_f_df(today, cts->fund_leg->cpn[k].pay_date, mkt->yc)) *
                cts->fund_leg->cpn[k].not ;

            k++;
          }

          onetime_strike /= onetime_level;
          onetime_notional = 1.0;
          onetime_pv = onetime_call[index_onetime] / onetime_notional;
          onetime_fwd /= onetime_level;
          onetime_mat = cts->call[index_call].ex_time;

          if (onetime_fwd < 0.0 && onetime_strike < 0.0 &&
              onetime_level < 0.0) {
            onetime_fwd *= -1.0;
            onetime_strike *= -1.0;
            onetime_level *= -1.0;
          }

          fwd_iv_info->one_time_strike[index_call] = onetime_strike;

          if (onetime_type == SRT_CALL) {
            /* we need to transform into put */
            onetime_pv =
                onetime_pv + (onetime_strike - onetime_fwd) * onetime_level;
            onetime_call[index_onetime] = onetime_pv * onetime_notional;
            european_prices[cts_prm->call->exo_idx] =
                onetime_call[index_onetime];
          }

          err = cts_optimpvol(onetime_pv, onetime_fwd, onetime_strike,
                              onetime_mat, onetime_level, SRT_PUT, SRT_NORMAL,
                              &onetime_vol);

          if (err) {
            onetime_vol = 0.0;
            /*
            one_time_calc[cts_prm->call->exo_idx] = 0;
            */
            err = NULL;
          }

          fwd_iv_info->one_time_norm_vol[index_call] = onetime_vol;

          if (onetime_fwd > 0.0 && onetime_strike > 0.0 && onetime_vol > 0.0) {
            err = cts_optimpvol(onetime_pv, onetime_fwd, onetime_strike,
                                onetime_mat, onetime_level, SRT_PUT,
                                SRT_LOGNORMAL, &onetime_vol);

            if (err && fwd_iv_info->one_time_norm_vol[index_call] > 1.0E-10) {
              onetime_vol = 0.0;
              normal_vol = fwd_iv_info->one_time_norm_vol[index_call] /
                           (onetime_fwd + onetime_strike) * 2.0;
              normal_vol = (normal_vol + vol_shift) *
                           (onetime_fwd + onetime_strike) / 2.0;
              /*
              one_time_calc[cts_prm->call->exo_idx] = 1;
              */

              european_prices_bumped[cts_prm->call->exo_idx] = srt_f_optblknrm(
                  onetime_fwd, onetime_strike, normal_vol, onetime_mat,
                  onetime_level, SRT_PUT, SRT_PREMIUM);
              err = NULL;
            } else {
              european_prices_bumped[cts_prm->call->exo_idx] = srt_f_optblksch(
                  onetime_fwd, onetime_strike, onetime_vol + vol_shift,
                  onetime_mat, onetime_level, SRT_PUT, SRT_PREMIUM);
            }
          } else {
            if (fwd_iv_info->one_time_norm_vol[index_call] > 1.0E-10) {
              onetime_vol = 0.0;
              normal_vol = fwd_iv_info->one_time_norm_vol[index_call] /
                           (onetime_fwd + onetime_strike) * 2.0;
              normal_vol = (normal_vol + vol_shift) *
                           (onetime_fwd + onetime_strike) / 2.0;
              /*
              one_time_calc[cts_prm->call->exo_idx] = 1;
              */

              european_prices_bumped[cts_prm->call->exo_idx] = srt_f_optblknrm(
                  onetime_fwd, onetime_strike, normal_vol, onetime_mat,
                  onetime_level, SRT_PUT, SRT_PREMIUM);
              err = NULL;
            } else {
              onetime_vol = 0.0;
              european_prices_bumped[cts_prm->call->exo_idx] = onetime_pv;
              /*
              one_time_calc[cts_prm->call->exo_idx] = 0;
              */
            }
          }

          *nb_one_time_calc += one_time_calc[cts_prm->call->exo_idx];

          fwd_iv_info->one_time_log_vol[index_call] = onetime_vol;

          index_onetime++;
        } else {
          fwd_iv_info->one_time_call[index_call] = 0.0;
          fwd_iv_info->one_time_strike[index_call] = 0.0;
          fwd_iv_info->one_time_norm_vol[index_call] = 0.0;
          fwd_iv_info->one_time_log_vol[index_call] = 0.0;
        }

        index_call++;
      }
    }
  }

FREE_RETURN:

  return err;
}

Err cts_convert_to_equi_midat(CTS_MKT mkt, CTS mkt_cts, CTS mdl_cts) {
  Err err = NULL;
  int i, j, last_index_funding, index_funding;
  double fund_pv, fund_cash_pv, coupon;
  double float_notional, float_ratio;
  long today;
  char *yc;
  CTS_CPN cpn;
  CTS_FUND_CPN fund_cpn;
  long fix;

  today = mkt->today;
  yc = mkt->yc;

  index_funding = mdl_cts->call[0].fund_idx;

  for (i = mdl_cts->call[0].exo_idx; i < mdl_cts->exo_leg->num_cpn; i++) {
    cpn = &(mdl_cts->exo_leg->cpn[i]);

    fund_pv = 0.0;

    last_index_funding = index_funding;

    while (index_funding < mkt_cts->fund_leg->num_cpn &&
           mkt_cts->fund_leg->cpn[index_funding].pay_date <
               mkt_cts->exo_leg->cpn[i].cxxpn_pay_date + 5) {
      fund_pv += mkt_cts->fund_leg->cpn[index_funding].mkt_val;

      /* Change the funding */

      fund_cpn = &(mdl_cts->fund_leg->cpn[index_funding]);

      index_funding++;
    }

    if (index_funding == mkt_cts->fund_leg->num_cpn) {
      fund_pv += swp_f_df(today, mkt_cts->exo_leg->cpn[i].cxxpn_pay_date, yc) *
                 mkt_cts->fund_leg->cpn[index_funding - 1].not ;
    }

    coupon = mkt_cts->exo_leg->cpn[i].mkt_fixed_part /
                 mkt_cts->exo_leg->cpn[i].cxxpn_not -
             fund_pv / mkt_cts->fund_leg->cpn[index_funding - 1].not ;

    coupon /= mkt_cts->exo_leg->cpn[i].cxxpn_cvg *
              swp_f_df(today, mkt_cts->exo_leg->cpn[i].cxxpn_pay_date, yc);

    /* adjust float notional */

    fund_cash_pv =
        swp_f_df(today, mkt_cts->exo_leg->cpn[i].cxxpn_start_date, yc) -
        swp_f_df(today, mkt_cts->exo_leg->cpn[i].cxxpn_pay_date, yc);
    float_notional = mkt_cts->exo_leg->cpn[i].mkt_float_part / fund_cash_pv;
    float_ratio =
        (mkt_cts->fund_leg->cpn[index_funding - 1].not -float_notional) /
        mkt_cts->fund_leg->cpn[index_funding - 1].not ;

    for (j = last_index_funding; j < index_funding; j++) {
      mdl_cts->fund_leg->cpn[j].not *= float_ratio;
    }

    /* kill the fixings of the cts */
    cts_free_coupon(cpn);

    /* realloc */

    cpn->type = 0;
    cpn->cpn_coupon = coupon;

    cpn->nfix = 1;
    cpn->used_nfix = 1;

    cpn->fix = (cts_fix *)calloc(cpn->used_nfix, sizeof(cts_fix));
    cpn->used_fixidx = (int *)calloc(cpn->used_nfix, sizeof(int));
    cpn->used_fixweights = (double *)calloc(cpn->used_nfix, sizeof(double));

    if (!cpn->fix || !cpn->used_fixidx || !cpn->used_fixweights) {
      err = "Allocation error (1) in cts_convert_to_equi_midat";
      goto FREE_RETURN;
    }

    cpn->used_fixidx[0] = 0;
    cpn->used_fixweights[0] = 1.0;

    fix = (long)(0.5 * (cpn->cpn_start_date + cpn->cpn_pay_date));

    cpn->fix[0].ref_fix_date = fix;
    cpn->fix[0].ref_fix_time = (fix - today) * YEARS_IN_DAY;
    cpn->fix[0].schedule.iNCpn = 0;

    err = cts_init_cpn_payoff_for_range(0, 0, -999, 999, 0, 9999999999, 1.0,
                                        10.0E-04, cpn);

    if (err)
      goto FREE_RETURN;
  }

  /* recalculate the funding leg */
  for (i = 0; i < mdl_cts->fund_leg->num_cpn; i++) {
    fund_cpn = &(mdl_cts->fund_leg->cpn[i]);
    fund_cpn->cpn = 0.0;
    fund_cpn->cpn_plus_ex = fund_cpn->cpn - fund_cpn->not ;

    if (i < mdl_cts->fund_leg->num_cpn - 1) {
      fund_cpn->cpn_plus_ex += (fund_cpn + 1)->not ;
    }
  }

  /* turn the call into midat */
  for (i = 0; i < mdl_cts->num_calls; i++) {
    mdl_cts->call[i].is_midat = 1;
  }

FREE_RETURN:

  return err;
}

Err cts_calibrate_equi_midat(CTS_MKT mkt, CTS_UND und, CTS cts,

                             char *exo_basis, char *fund_basis,

                             int notperiod, int *one_time_calc, double *prices,
                             double *fees,

                             double lambda,

                             double numer_tstar, double min_time) {
  Err err = NULL;
  long *fix_start_dates = NULL, *fix_start_dates_one = NULL,
       *fix_end_dates = NULL, *fix_end_dates_one = NULL,

       *fix_pay_dates = NULL, *fix_pay_dates_one = NULL,
       *float_pay_dates = NULL, *float_pay_dates_one = NULL,

       *float_start_dates = NULL, *float_start_dates_one = NULL,
       *float_end_dates = NULL, *float_end_dates_one = NULL;

  double *fix_rates = NULL, *fix_rates_one = NULL, *fix_notionals = NULL,
         *fix_notionals_one = NULL, *ex_fee = NULL, *ex_fee_one = NULL,
         *ex_prices = NULL, *ex_prices_one = NULL,

         *float_margins = NULL, *float_margins_one = NULL,
         *float_spreads = NULL, *float_spreads_one = NULL,
         *float_notionals = NULL, *float_notionals_one = NULL;

  int *ex_bool = NULL, *ex_bool_one = NULL;

  double step_calib;
  int i, nb_new_dates, nb_months;
  char exo_freq[2], fund_freq[2];
  int min_days;
  int next_i, last_i;

  /* Freq calculation for Betrand !!! */
  i = cts->exo_leg->num_cpn - 1;
  nb_months = 13;

  while (i >= 0 && nb_months != 1 && nb_months != 3 && nb_months != 6 &&
         nb_months != 12) {
    nb_months = (int)((cts->exo_leg->cpn[i].cxxpn_pay_date -
                       cts->exo_leg->cpn[i].cxxpn_start_date) /
                          30.0 +
                      0.5);
    i--;
  }

  switch (nb_months) {
  case 1:
    strcpy(exo_freq, "M");
    break;

  case 3:
    strcpy(exo_freq, "Q");
    break;

  case 6:
    strcpy(exo_freq, "S");
    break;

  case 12:
    strcpy(exo_freq, "A");
    break;
  }

  i = cts->fund_leg->num_cpn - 1;
  nb_months = 13;

  while (i >= 0 && nb_months != 1 && nb_months != 3 && nb_months != 6 &&
         nb_months != 12) {
    nb_months = (int)((cts->fund_leg->cpn[i].pay_date -
                       cts->fund_leg->cpn[i].start_date) /
                          30.0 +
                      0.5);
    i--;
  }

  switch (nb_months) {
  case 1:
    strcpy(fund_freq, "M");
    break;

  case 3:
    strcpy(fund_freq, "Q");
    break;

  case 6:
    strcpy(fund_freq, "S");
    break;

  case 12:
    strcpy(fund_freq, "A");
    break;
  }

  /* Recopy Fixed Leg */
  fix_start_dates = calloc(cts->exo_leg->num_cpn, sizeof(long));
  fix_start_dates_one = calloc(cts->exo_leg->num_cpn, sizeof(long));
  fix_end_dates = calloc(cts->exo_leg->num_cpn, sizeof(long));
  fix_end_dates_one = calloc(cts->exo_leg->num_cpn, sizeof(long));
  fix_pay_dates = calloc(cts->exo_leg->num_cpn, sizeof(long));
  fix_pay_dates_one = calloc(cts->exo_leg->num_cpn, sizeof(long));

  fix_rates = calloc(cts->exo_leg->num_cpn, sizeof(double));
  fix_rates_one = calloc(cts->exo_leg->num_cpn, sizeof(double));
  fix_notionals = calloc(cts->exo_leg->num_cpn, sizeof(double));
  fix_notionals_one = calloc(cts->exo_leg->num_cpn, sizeof(double));

  ex_fee = calloc(cts->exo_leg->num_cpn, sizeof(double));
  ex_fee_one = calloc(cts->exo_leg->num_cpn, sizeof(double));
  ex_bool = calloc(cts->exo_leg->num_cpn, sizeof(int));
  ex_bool_one = calloc(cts->exo_leg->num_cpn, sizeof(int));
  ex_prices = calloc(cts->exo_leg->num_cpn, sizeof(double));
  ex_prices_one = calloc(cts->exo_leg->num_cpn, sizeof(double));

  float_start_dates = calloc(cts->fund_leg->num_cpn, sizeof(long));
  float_start_dates_one = calloc(cts->fund_leg->num_cpn, sizeof(long));
  float_end_dates = calloc(cts->fund_leg->num_cpn, sizeof(long));
  float_end_dates_one = calloc(cts->fund_leg->num_cpn, sizeof(long));
  float_pay_dates = calloc(cts->fund_leg->num_cpn, sizeof(long));
  float_pay_dates_one = calloc(cts->fund_leg->num_cpn, sizeof(long));

  float_margins = calloc(cts->fund_leg->num_cpn, sizeof(double));
  float_margins_one = calloc(cts->fund_leg->num_cpn, sizeof(double));
  float_spreads = calloc(cts->fund_leg->num_cpn, sizeof(double));
  float_spreads_one = calloc(cts->fund_leg->num_cpn, sizeof(double));
  float_notionals = calloc(cts->fund_leg->num_cpn, sizeof(double));
  float_notionals_one = calloc(cts->fund_leg->num_cpn, sizeof(double));

  if (!fix_start_dates || !fix_end_dates || !fix_rates || !fix_notionals ||
      !ex_fee || !ex_bool || !ex_prices || !float_start_dates ||
      !float_end_dates || !float_margins || !float_spreads ||
      !float_notionals || !fix_start_dates_one || !fix_end_dates_one ||
      !fix_rates_one || !fix_notionals_one || !ex_fee_one || !ex_bool ||
      !ex_prices_one || !fix_end_dates || !fix_end_dates_one ||
      !float_start_dates_one || !float_end_dates_one || !float_margins_one ||
      !float_spreads_one || !float_notionals_one || !float_pay_dates ||
      !float_pay_dates_one) {
    err = "Memory allocation faillure (1) in cts_calibrate_equi_midat";
    goto FREE_RETURN;
  }

  for (i = 0; i < cts->exo_leg->num_cpn; i++) {
    fix_start_dates[i] = cts->exo_leg->cpn[i].cxxpn_start_date;
    fix_start_dates_one[i] = fix_start_dates[i];

    fix_end_dates[i] = cts->exo_leg->cpn[i].cxxpn_end_date;
    fix_end_dates_one[i] = fix_end_dates[i];

    fix_pay_dates[i] = cts->exo_leg->cpn[i].cxxpn_pay_date;
    fix_pay_dates_one[i] = fix_pay_dates[i];

    fix_rates[i] = cts->exo_leg->cpn[i].cxxpn_coupon;
    fix_rates_one[i] = fix_rates[i];

    fix_notionals[i] = cts->exo_leg->cpn[i].cxxpn_not;
    fix_notionals_one[i] = fix_notionals[i];

    ex_fee[i] = fees[i];
    ex_fee_one[i] = ex_fee[i];

    ex_bool[i] = one_time_calc[i];
    ex_bool_one[i] = 0;

    ex_prices[i] = prices[i];
    ex_prices_one[i] = ex_prices[i];
  }

  for (i = 0; i < cts->fund_leg->num_cpn; i++) {
    float_start_dates[i] = cts->fund_leg->cpn[i].start_date;
    float_start_dates_one[i] = float_start_dates[i];

    float_end_dates[i] = cts->fund_leg->cpn[i].end_date;
    float_end_dates_one[i] = float_end_dates[i];

    float_pay_dates[i] = float_end_dates[i];
    float_pay_dates_one[i] = float_end_dates[i];

    float_margins[i] = 0.0;
    float_margins_one[i] = float_margins[i];

    float_spreads[i] = 0.0;
    float_spreads_one[i] = float_spreads[i];

    float_notionals[i] = cts->fund_leg->cpn[i].not ;
    float_notionals_one[i] = float_notionals[i];
  }

  i = 0;
  while (i < cts->exo_leg->num_cpn && !ex_bool[i]) {
    i++;
  }

  if (cts->call[i].ex_time < 1.0E-08) {
    i++;
  }

  /* apply min_time */
  min_days = (int)(min_time * DAYS_IN_YEAR + 0.5);

  last_i = i;
  next_i = i + 1;

  while (next_i < cts->exo_leg->num_cpn) {
    if (fix_start_dates[next_i] - fix_start_dates[last_i] < min_days) {
      ex_bool[next_i] = 0;
    } else {
      last_i = next_i;
    }

    next_i++;
  }

  /* First run with only the first one-time */

  if (i < cts->exo_leg->num_cpn) {
    ex_bool_one[i] = 1;

    err = amortMidat_cpd_calib_diagonal_new(
        notperiod, mkt->yc, mkt->vc, mkt->ref, mkt->swap_basis, mkt->swap_freq,

        mkt->get_cash_vol,

        0, 0,

        ex_bool_one, ex_prices_one, ex_fee_one,

        exo_freq, exo_basis, cts->exo_leg->num_cpn, fix_start_dates_one,
        fix_end_dates_one, fix_pay_dates_one, fix_rates_one, fix_notionals_one,

        fund_freq, fund_basis, cts->fund_leg->num_cpn, float_start_dates_one,
        float_end_dates_one, float_pay_dates_one, float_margins_one,
        float_spreads_one, float_notionals_one,

        NULL, 0, 999,

        1, 0, 0, 1, 0.0001, &lambda, 1, 0, 0, 0, &(und->model.iNbPWTime),
        &(und->model.dPWTime), &(und->model.dSigma));

    if (err)
      goto FREE_RETURN;

    if (und->model.dSigma[0] > MAX_FIRST_VOL) {
      ex_bool[i] = 0;
      one_time_calc[i] = 0;
    }

    if (und->model.dPWTime)
      free(und->model.dPWTime);
    und->model.dPWTime = NULL;
    if (und->model.dSigma)
      free(und->model.dSigma);
    und->model.dSigma = NULL;
  }

  /* do the regular calibration */

  /*	add discretisation dates till the end of the deal */
  err = amortMidat_cpd_calib_diagonal_new(
      notperiod, mkt->yc, mkt->vc, mkt->ref, mkt->swap_basis, mkt->swap_freq,

      mkt->get_cash_vol,

      0, 0,

      ex_bool, ex_prices, ex_fee,

      exo_freq, exo_basis, cts->exo_leg->num_cpn, fix_start_dates,
      fix_end_dates, fix_pay_dates, fix_rates, fix_notionals,

      fund_freq, fund_basis, cts->fund_leg->num_cpn, float_start_dates,
      float_end_dates, float_pay_dates, float_margins, float_spreads,
      float_notionals,

      NULL, 0, 999,

      1, 0, 0, 1, 0.0001, &lambda, 1, 0, 0, 0, &(und->model.iNbPWTime),
      &(und->model.dPWTime), &(und->model.dSigma));

  if (err)
    goto FREE_RETURN;

  if (und->model.iNbPWTime > 1) {
    step_calib = und->model.dPWTime[1] - und->model.dPWTime[0];
  } else {
    step_calib = 1.0;
  }

  nb_new_dates =
      (int)((cts->exo_leg->cpn[cts->exo_leg->num_cpn - 1].cxxpn_pay_time -
             und->model.dPWTime[und->model.iNbPWTime - 1]) /
            step_calib) +
      1;

  und->model.iNbPWTime += nb_new_dates;

  und->model.dPWTime = (double *)realloc(und->model.dPWTime,
                                         und->model.iNbPWTime * sizeof(double));
  und->model.dSigma = (double *)realloc(und->model.dSigma,
                                        und->model.iNbPWTime * sizeof(double));

  und->model.dAlpha = (double *)calloc(und->model.iNbPWTime, sizeof(double));
  und->model.dRho = (double *)calloc(und->model.iNbPWTime, sizeof(double));
  und->model.dLambdaEps =
      (double *)calloc(und->model.iNbPWTime, sizeof(double));
  und->model.dLvlEps = (double *)calloc(und->model.iNbPWTime, sizeof(double));

  if (!und->model.dPWTime || !und->model.dSigma || !und->model.dAlpha ||
      !und->model.dRho || !und->model.dLambdaEps || !und->model.dLvlEps) {
    err = "Memory allocation error (2) in cts_calibrate_equi_midat";
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_new_dates; i++) {
    und->model.dPWTime[und->model.iNbPWTime - nb_new_dates + i] =
        und->model.dPWTime[und->model.iNbPWTime - nb_new_dates + i - 1] +
        step_calib;
    und->model.dSigma[und->model.iNbPWTime - nb_new_dates + i] =
        und->model.dSigma[und->model.iNbPWTime - nb_new_dates + i - 1];
  }

  for (i = 0; i < und->model.iNbPWTime; i++) {
    und->model.dAlpha[i] = 0.01;
    und->model.dRho[i] = 0.0;
    und->model.dLambdaEps[i] = 0.0;
    und->model.dLvlEps[i] = 0.0;
  }

  und->model.lToday = mkt->today;
  und->model.dLambdaX = lambda;
  und->model.dTau = 1.0 / lambda;
  und->model.dInitTStar = LGMSV_Tstar;
  und->model.dTStar = LGMSV_Tstar;
  und->model.dInitTStar = LGMSV_Tstar;
  und->model.lTStarDate = (long)(mkt->today + LGMSV_Tstar * DAYS_IN_YEAR + 0.5);

  ConvertTS_LGM_to_LGMSV(und->model.iNbPWTime, und->model.dPWTime,
                         und->model.dSigma, und->model.dLambdaX,
                         und->model.dTStar, 1, 0.0, 0.0, 0.0, NULL, NULL);

  Convert_Tstar_model(&(und->model), numer_tstar);

FREE_RETURN:

  if (fix_start_dates)
    free(fix_start_dates);
  if (fix_start_dates_one)
    free(fix_start_dates_one);

  if (fix_end_dates)
    free(fix_end_dates);
  if (fix_end_dates_one)
    free(fix_end_dates_one);

  if (fix_pay_dates)
    free(fix_pay_dates);
  if (fix_pay_dates_one)
    free(fix_pay_dates_one);

  if (fix_rates)
    free(fix_rates);
  if (fix_rates_one)
    free(fix_rates_one);

  if (fix_notionals)
    free(fix_notionals);
  if (fix_notionals_one)
    free(fix_notionals_one);

  if (ex_fee)
    free(ex_fee);
  if (ex_fee_one)
    free(ex_fee_one);

  if (ex_bool)
    free(ex_bool);
  if (ex_bool_one)
    free(ex_bool_one);

  if (ex_prices)
    free(ex_prices);
  if (ex_prices_one)
    free(ex_prices_one);

  if (float_start_dates)
    free(float_start_dates);
  if (float_start_dates_one)
    free(float_start_dates_one);

  if (float_end_dates)
    free(float_end_dates);
  if (float_end_dates_one)
    free(float_end_dates_one);

  if (float_pay_dates)
    free(float_pay_dates);
  if (float_pay_dates_one)
    free(float_pay_dates_one);

  if (float_margins)
    free(float_margins);
  if (float_margins_one)
    free(float_margins_one);

  if (float_spreads)
    free(float_spreads);
  if (float_spreads_one)
    free(float_spreads_one);

  if (float_notionals)
    free(float_notionals);
  if (float_notionals_one)
    free(float_notionals_one);

  return err;
}

Err cts_calc_reserve_midat(/*	Initial results */
                           CTS cts_mkt, CTS cts_mdl, CTS_UND und,
                           CTS_ADI_ARG adi_arg_ts,

                           double call, double *onetime_call,

                           /*	Calibration parameters */
                           int recalib_euro,

                           LGMSV_NUMERPARAMS NumerParams, double tstar,
                           double numer_tstar, char *cal_tenor, char *cal_ref,
                           char *cal_freq, char *cal_basis, int force_atm,
                           double max_std_long, double max_std_short,
                           double vol_shift_long,
                           DIAGCALIB_VOLTYPE vol_type_long,
                           DIAGCALIB_SHIFTTYPE vol_shift_type_long,
                           double vol_shift_short,
                           DIAGCALIB_VOLTYPE vol_type_short,
                           DIAGCALIB_SHIFTTYPE vol_shift_type_short,
                           double lambda_shift,
                           int calib_strategy, /*	-1: autocal        , 0:
                                            swaptions / cap        , 1: cap /
                                            swaptions */
                           int fix_lambda, char *short_tenor,
                           char *short_refrate, char *short_freq,
                           char *short_basis, int fix_smile,
                           int smile_calib_months, /* 0: co-terminal swaption ,
                                                otherwise underlyings with
                                                required nb months */
                           LGMSV_CalibParams *lgmsv_calib_params,
                           double min_time,
                           int skip_last,   /*	If 1        , the last option is
                                         disregarded and the forward
                                         volatility is flat from option n-1
                                       */
                           double min_fact, /*	Maximum down jump on variance */
                           double max_fact, /*	Maximum up jump on variance */
                           int use_jumps,   /*	1: we allow jumps on vol   , 0:
                                         we don't */
                           double prec, int maxiter, int keep_first,
                           int long_strike_flag,  /*	0: ATM 1: Coupon 2: Eq
                                               (PV/Lvl) */
                           int short_strike_flag, /*	0: ATM        , 1:
                                               implied digital caplet strike 2:
                                               same number of std */

                           /*	Reserve parameters */
                           double lambda_init, double lambda_reserve,
                           double onetime_vega,

                           int reserve_nstpt, /*	number of time steps for
                                           reserve calculation */
                           int reserve_nstpx, /*	number of space steps
                                           for reserve calculation */
                           double reserve_integ_mintime,

                           /*	Results */
                           double *switch_reserve, double *onetime_reserve,

                           CTS_IV fwd_iv_info, int save_extra_infos,
                           CTS_EXTRA_INFOS extra_infos) {
  CTS_MKT mkt;

  int notperiod;
  long temp_date;

  char *fund_basis, *exo_basis;

  cts_adi_arg *adi_arg_init = NULL, *adi_arg_bump = NULL,
              *adi_arg_switch = NULL;

  int free_adi_arg_init = 0, free_adi_arg_bump = 0, free_adi_arg_switch = 0;

  double call_init, switch_init;
  double call_new, switch_new;

  int nb_one_time_calc;
  int *one_time_calc = NULL;
  double *european_prices = NULL, *european_prices_vega = NULL,
         *european_fees = NULL;
  double most_expensive;

  double smilepartime, alphaeps, ldaeps, rhoeps;
  double one_time_call_temp[MAX_ONETIME];

  Err err = NULL;

  /* Check if calculation is required */
  if (fabs(onetime_vega) < 1.0E-08 &&
      fabs(lambda_reserve - lambda_init) < 1.0E-08) {
    *onetime_reserve = 0.0;
    *switch_reserve = 0.0;
    goto FREE_RETURN;
  }

  /* Setup LGM params */
  smilepartime = und->model.dPWTime[0];
  alphaeps = 0.005;
  ldaeps = 0.0;
  rhoeps = 0.0;

  /* Get the market */
  mkt = und->mkt;

  /*	Step 1: Convert the CTS into the equivalent midat */
  /*	it finds the fixed coupon which have the same PV  */
  /*	as the actual TS coupons. */
  /*	************************************************* */

  err = cts_convert_to_equi_midat(mkt, cts_mkt, cts_mdl);

  if (err)
    goto FREE_RETURN;

  err = cts_adjust_model_fwd_iv(und, cts_mkt, cts_mdl, 0, 0, 0, 0, 0, NULL);

  if (err)
    goto FREE_RETURN;

  err = translate_basis(
      &fund_basis,
      cts_mkt->fund_leg->cpn[cts_mkt->fund_leg->num_cpn - 1].basis);

  if (err)
    goto FREE_RETURN;

  err = translate_basis(
      &exo_basis,
      cts_mkt->exo_leg->cpn[cts_mkt->exo_leg->num_cpn - 1].cxxpn_basis);

  if (err)
    goto FREE_RETURN;

  notperiod = 0;
  temp_date = cts_mdl->exo_leg->cpn[cts_mdl->call[0].exo_idx].cxxpn_start_date;

  while (temp_date > cts_mdl->call[0].ex_date) {
    notperiod++;
    temp_date = add_unit(temp_date, -1, SRT_BDAY, MODIFIED_SUCCEEDING);
  }

  /* Get the one time infos */

  one_time_calc = calloc(cts_mdl->exo_leg->num_cpn, sizeof(int));
  european_prices = calloc(cts_mdl->exo_leg->num_cpn, sizeof(double));
  european_prices_vega = calloc(cts_mdl->exo_leg->num_cpn, sizeof(double));
  european_fees = calloc(cts_mdl->exo_leg->num_cpn, sizeof(double));

  adi_arg_init = calloc(1, sizeof(cts_adi_arg));
  adi_arg_bump = calloc(1, sizeof(cts_adi_arg));
  adi_arg_switch = calloc(1, sizeof(cts_adi_arg));

  if (!one_time_calc || !european_prices || !european_prices_vega ||
      !european_fees || !adi_arg_init || !adi_arg_bump || !adi_arg_switch) {
    err = "memory allocation faillure in cts_calc_reserve";
    goto FREE_RETURN;
  }

  err = cts_get_one_time_infos(mkt, adi_arg_ts, cts_mdl, onetime_vega,
                               onetime_call, one_time_calc, european_prices,
                               european_prices_vega, european_fees,
                               &nb_one_time_calc, &most_expensive, fwd_iv_info);

  if (err)
    goto FREE_RETURN;

  if (nb_one_time_calc == 0) {
    *onetime_reserve = 0.0;
    *switch_reserve = 0.0;
    goto FREE_RETURN;
  }

  /*	Step 2: Pricing of the reference bermudean */
  /*	****************************************** */

  free_LGMSV_model(&(und->model));

  if (recalib_euro) {
    err = cts_calibrate_equi_midat(
        mkt, und, cts_mdl, exo_basis, fund_basis, notperiod, one_time_calc,
        european_prices, european_fees, lambda_init, numer_tstar, min_time);
  } else {
    err = cts_calib_und(
        mkt, 0, lambda_init, 1, 0.0, 0.0, 0.0, 1, &smilepartime, &alphaeps,
        &rhoeps, &ldaeps, NULL, tstar, NumerParams, cal_tenor, cal_ref,
        cal_freq, cal_basis, force_atm, max_std_long, max_std_short,
        vol_shift_long, vol_type_long, vol_shift_type_long, vol_shift_short,
        vol_type_short, vol_shift_type_short, lambda_shift, calib_strategy, 0,
        short_tenor, short_refrate, short_freq, short_basis, 1,
        smile_calib_months, lgmsv_calib_params, min_time, skip_last, min_fact,
        max_fact, use_jumps, prec, maxiter, keep_first, 1, 0, cts_mdl, und, 0,
        NULL);
  }

  if (err)
    goto FREE_RETURN;

  /* recalibrate the adi_arg */
  err = cts_fill_algo_arg(und, cts_mdl, 0, reserve_nstpt, 1, reserve_nstpx, 1,
                          0, 0, reserve_integ_mintime, 0, 0, adi_arg_init);

  free_adi_arg_init = 1;

  if (err)
    goto FREE_RETURN;

  smessage(
      "Launching adi        , time steps requested: %d        , actual: %d",
      reserve_nstpt, adi_arg_init->nstp);

  err = cts_launch_algo(cts_mdl, und, adi_arg_init, &call_init,
                        &(one_time_call_temp[0]), NULL);

  if (err)
    goto FREE_RETURN;

  switch_init = call_init - most_expensive;

  extra_infos->call1 = call_init;

  /*	Step 3: Pricing of the bermudean with bumped one time callable */
  /*	************************************************************** */

  if (fabs(onetime_vega) > 1.0E-08) {
    free_LGMSV_model(&(und->model));

    if (recalib_euro) {
      err = cts_calibrate_equi_midat(mkt, und, cts_mdl, exo_basis, fund_basis,
                                     notperiod, one_time_calc,
                                     european_prices_vega, european_fees,
                                     lambda_init, numer_tstar, min_time);
    } else {
      err = cts_calib_und(
          mkt, 0, lambda_init, 1, 0.0, 0.0, 0.0, 1, &smilepartime, &alphaeps,
          &rhoeps, &ldaeps, NULL, tstar, NumerParams, cal_tenor, cal_ref,
          cal_freq, cal_basis, force_atm, max_std_long, max_std_short,
          vol_shift_long, vol_type_long, vol_shift_type_long, vol_shift_short,
          vol_type_short, vol_shift_type_short, lambda_shift, calib_strategy, 0,
          short_tenor, short_refrate, short_freq, short_basis, 1,
          smile_calib_months, lgmsv_calib_params, min_time, skip_last, min_fact,
          max_fact, use_jumps, prec, maxiter, keep_first, 1, 0, cts_mdl, und, 0,
          NULL);
    }

    if (err)
      goto FREE_RETURN;

    /* recalibrate the adi_arg */
    err = cts_fill_algo_arg(und, cts_mdl, 0, reserve_nstpt, 1, reserve_nstpx, 1,
                            0, 0, reserve_integ_mintime, 0, 0, adi_arg_bump);

    free_adi_arg_bump = 1;

    if (err)
      goto FREE_RETURN;

    smessage(
        "Launching adi        , time steps requested: %d        , actual: %d",
        reserve_nstpt, adi_arg_bump->nstp);

    err = cts_launch_algo(cts_mdl, und, adi_arg_bump, &call_new,
                          &(one_time_call_temp[0]), NULL);

    if (err)
      goto FREE_RETURN;

    *onetime_reserve = call_new - call_init;
  } else {
    *onetime_reserve = 0.0;
  }

  /*	Step 4: Pricing of the bermudean with the new lambda */
  /*	**************************************************** */

  if (fabs(lambda_reserve - lambda_init) > 1.0E-08 && nb_one_time_calc > 1) {
    free_LGMSV_model(&(und->model));

    if (recalib_euro) {
      err = cts_calibrate_equi_midat(mkt, und, cts_mdl, exo_basis, fund_basis,
                                     notperiod, one_time_calc, european_prices,
                                     european_fees, lambda_reserve, numer_tstar,
                                     min_time);
    } else {
      err = cts_calib_und(
          mkt, 0, lambda_reserve, 1, 0.0, 0.0, 0.0, 1, &smilepartime, &alphaeps,
          &rhoeps, &ldaeps, NULL, tstar, NumerParams, cal_tenor, cal_ref,
          cal_freq, cal_basis, force_atm, max_std_long, max_std_short,
          vol_shift_long, vol_type_long, vol_shift_type_long, vol_shift_short,
          vol_type_short, vol_shift_type_short, lambda_shift, calib_strategy, 0,
          short_tenor, short_refrate, short_freq, short_basis, 1,
          smile_calib_months, lgmsv_calib_params, min_time, skip_last, min_fact,
          max_fact, use_jumps, prec, maxiter, keep_first, 1, 0, cts_mdl, und, 0,
          NULL);
    }

    if (err)
      goto FREE_RETURN;

    /* recalibrate the adi_arg */
    err = cts_fill_algo_arg(und, cts_mdl, 0, reserve_nstpt, 1, reserve_nstpx, 1,
                            0, 0, reserve_integ_mintime, 0, 0, adi_arg_switch);

    free_adi_arg_switch = 1;

    if (err)
      goto FREE_RETURN;

    smessage(
        "Launching adi        , time steps requested: %d        , actual: %d",
        reserve_nstpt, adi_arg_switch->nstp);

    err = cts_launch_algo(cts_mdl, und, adi_arg_switch, &call_new,
                          &(one_time_call_temp[0]), NULL);

    if (err)
      goto FREE_RETURN;

    switch_new = call_new - most_expensive;

    *switch_reserve = switch_new - switch_init;
  } else {
    *switch_reserve = 0.0;
  }

FREE_RETURN:

  if (adi_arg_init) {
    if (free_adi_arg_init)
      cts_free_adi_arg(adi_arg_init);
    free(adi_arg_init);
  }

  if (adi_arg_bump) {
    if (free_adi_arg_bump)
      cts_free_adi_arg(adi_arg_bump);
    free(adi_arg_bump);
  }

  if (adi_arg_switch) {
    if (free_adi_arg_switch)
      cts_free_adi_arg(adi_arg_switch);
    free(adi_arg_switch);
  }

  if (one_time_calc)
    free(one_time_calc);
  if (european_prices)
    free(european_prices);
  if (european_prices_vega)
    free(european_prices_vega);
  if (european_fees)
    free(european_fees);

  return err;
}

Err cts_calc_reserve_timeswap(/*	Initial results */
                              CTS cts_mkt, CTS cts_mdl, CTS_UND und,

                              double call, double *onetime_call,
                              int nb_european,

                              int for_fund,

                              /*	Model parameters */
                              int nsmilepar, double *smilepartime,
                              double *alphaeps, double *rhoeps, double *ldaeps,
                              double tstar,

                              /*	Calibration parameters */
                              LGMSV_NUMERPARAMS NumerParams, double numer_tstar,
                              char *cal_tenor, char *cal_ref, char *cal_freq,
                              char *cal_basis, int force_atm,
                              double max_std_long, double max_std_short,
                              double vol_shift_long,
                              DIAGCALIB_VOLTYPE vol_type_long,
                              DIAGCALIB_SHIFTTYPE vol_shift_type_long,
                              double vol_shift_short,
                              DIAGCALIB_VOLTYPE vol_type_short,
                              DIAGCALIB_SHIFTTYPE vol_shift_type_short,
                              double lambda_shift,
                              int calib_strategy, /*	-1: autocal        , 0:
                                               swaptions / cap        , 1: cap /
                                               swaptions */
                              int fix_lambda, char *short_tenor,
                              char *short_refrate, char *short_freq,
                              char *short_basis, int fix_smile,
                              int smile_calib_months, /* 0: co-terminal swaption
                                                         , otherwise underlyings
                                                   with required nb months
                                                 */
                              LGMSV_CalibParams
                                  *lgmsv_calib_params, /*	shift on Rho SV
                                                        */
                              double min_time,
                              int skip_last, /*	If 1        , the last option is
                                          disregarded and the forward
                                          volatility is flat from option
                                          n-1 */
                              double
                                  min_fact, /*	Maximum down jump on variance */
                              double
                                  max_fact,  /*	Maximum up jump on variance */
                              int use_jumps, /*	1: we allow jumps on vol , 0:
                                          we don't */
                              double prec, int maxiter, int keep_first,
                              int long_strike_flag,  /*	0: ATM 1: Coupon 2: Eq
                                                  (PV/Lvl) */
                              int short_strike_flag, /*	0: ATM        , 1:
                                                  implied digital caplet strike
                                                  2: same number of std */

                              int calc_fwd_iv, int adj_fee,

                              /*	Reserve parameters */
                              double lambda_init, double lambda_reserve,
                              double onetime_vega,

                              int euro_nb_iter,   /*	Number of Levenberg
                                               iteration for calibration of
                                               european */
                              int one_time_index, /*	One Time callable used
                                               for reserve calculation */

                              int reserve_nstpt,   /*	number of time steps for
                                                reserve calculation */
                              int reserve_nstppsi, /*	number of psi steps for
                                                reserve calculation */
                              int reserve_nstpx,   /*	number of space steps
                                                for reserve calculation */
                              int reserve_nstpvol, /*	number of vol steps for
                                                reserve calculation */
                              double reserve_integ_mintime,

                              /*	Results */
                              double *switch_reserve, double *onetime_reserve,

                              CTS_IV fwd_iv_info, int save_extra_infos,
                              CTS_EXTRA_INFOS extra_infos) {
  CTS_MKT mkt;

  double fitting_error;

  cts_adi_arg *adi_arg = NULL;

  int free_adi_arg = 0;
  int i;

  double most_expensive, switch_init, call_new, switch_new;

  Err err = NULL;

  /* Memory allocation */
  adi_arg = calloc(1, sizeof(cts_adi_arg));

  if (!adi_arg) {
    err = "Memory allocation faillure in cts_calc_reserve_timeswap";
    goto FREE_RETURN;
  }

  /* Get the market */
  mkt = und->mkt;

  /*	Step 1: Get the most expensive PV  */
  /*	********************************** */

  most_expensive = onetime_call[0];

  for (i = 1; i < nb_european; i++) {
    if (onetime_call[i] > most_expensive) {
      most_expensive = onetime_call[i];
    }
  }

  switch_init = call - most_expensive;

  /*	Step 2: Pricing of the switch reserve by recalibrating european time
   * swap */
  /*	*************************************************************************
   */

  free_LGMSV_model(&(und->model));

  /* first recalibrate the model with the new lambda */
  err = cts_calib_und(
      mkt, 0, lambda_reserve, 1, 0.0, 0.0, 0.0, nsmilepar, smilepartime,
      alphaeps, rhoeps, ldaeps, NULL, tstar, NumerParams, cal_tenor, cal_ref,
      cal_freq, cal_basis, force_atm, max_std_long, max_std_short,
      vol_shift_long, vol_type_long, vol_shift_type_long, vol_shift_short,
      vol_type_short, vol_shift_type_short, lambda_shift, calib_strategy,
      (fabs(lambda_reserve) > 1.0E-08), short_tenor, short_refrate, short_freq,
      short_basis, fix_smile, smile_calib_months, lgmsv_calib_params, min_time,
      skip_last, min_fact, max_fact, use_jumps, prec, maxiter, keep_first,
      long_strike_flag, short_strike_flag, cts_mdl, und, 0, NULL);

  if (err)
    goto FREE_RETURN;

  /* recalibrate the adi_arg with the numer tstar */

  Convert_Tstar_model(&(und->model), numer_tstar);

  err = cts_fill_algo_arg(und, cts_mdl, 0, reserve_nstpt, reserve_nstppsi,
                          reserve_nstpx, reserve_nstpvol, 0, 0,
                          reserve_integ_mintime, 1, one_time_index, adi_arg);

  free_adi_arg = 1;

  Convert_Tstar_model(&(und->model), und->model.dInitTStar);

  if (err)
    goto FREE_RETURN;

  err = cts_calibrate_european_and_price(
      und, cts_mkt, cts_mdl, for_fund, calc_fwd_iv, adj_fee, adi_arg,
      NumerParams, numer_tstar, nb_european, onetime_call, euro_nb_iter,
      &call_new, &fitting_error);

  if (err)
    goto FREE_RETURN;

  switch_new = call_new - most_expensive;

  *onetime_reserve = 0.0;
  *switch_reserve = switch_new - switch_init;

FREE_RETURN:

  if (adi_arg) {
    if (free_adi_arg)
      cts_free_adi_arg(adi_arg);
    free(adi_arg);
  }

  return err;
}

Err cts_calc_reserve(/*	Initial results */
                     CTS cts_mkt, CTS cts_mdl, CTS_UND und, CTS_ADI_ARG adi_arg,

                     int for_fund,

                     double call, double *onetime_call,

                     /*	Model parameters */
                     int nsmilepar, double *smilepartime, double *alphaeps,
                     double *rhoeps, double *ldaeps, double *rho2eps,
                     double tstar,

                     /*	Calibration parameters */
                     LGMSV_NUMERPARAMS NumerParams, double numer_tstar,
                     char *cal_tenor, char *cal_ref, char *cal_freq,
                     char *cal_basis, int force_atm, double max_std_long,
                     double max_std_short, double vol_shift_long,
                     DIAGCALIB_VOLTYPE vol_type_long,
                     DIAGCALIB_SHIFTTYPE vol_shift_type_long,
                     double vol_shift_short, DIAGCALIB_VOLTYPE vol_type_short,
                     DIAGCALIB_SHIFTTYPE vol_shift_type_short,
                     double lambda_shift,
                     int calib_strategy, /*	-1: autocal        , 0:
                                      swaptions / cap        , 1: cap /
                                      swaptions */
                     int fix_lambda, char *short_tenor, char *short_refrate,
                     char *short_freq, char *short_basis, int fix_smile,
                     int smile_calib_months, /* 0: co-terminal swaption        ,
                                          otherwise underlyings with
                                          required nb months */
                     LGMSV_CalibParams
                         *lgmsv_calib_params, /*	shift on Rho SV */
                     double min_time,
                     int skip_last,   /*	If 1        , the last option is
                                   disregarded   and the forward volatility is flat
                                   from   option n-1 */
                     double min_fact, /*	Maximum down jump on variance */
                     double max_fact, /*	Maximum up jump on variance */
                     int use_jumps,   /*	1: we allow jumps on vol        , 0: we
                                       * don't
                                       */
                     double prec, int maxiter, int keep_first,
                     int long_strike_flag,  /*	0: ATM 1: Coupon 2: Eq (PV/Lvl)
                                             */
                     int short_strike_flag, /*	0: ATM        , 1: implied
                                         digital caplet strike 2: same number of
                                         std */

                     int calc_fwd_iv, int adj_fee,

                     /*	Reserve parameters */
                     double lambda_reserve, double lgm_alpha, double lgm_gamma,
                     double lgm_rho, double onetime_vega, int vega_shift_type,

                     int reserve_method,    /*	0: old method        , 1: new
                                             * method
                                             */
                     int lgm_reserve,       /*	1: calculate the reserve in the
                                         one factor */
                     int midat_reserve,     /*	1: reserve calculated using the
                                         midat-replication */
                     int recalib_european,  /*	If set to 1        , then
                                         european are  recalibrated when changing
                                         lambda
                                       */
                     int euro_nb_iter,      /*	Number of Levenberg iteration
                                         for calibration of european */
                     int one_time_index,    /*	One Time callable used for
                                         reserve calculation */
                     int recalc_one_factor, /*	Recalculate MCEB in one factor
                                         for reserve */

                     int reserve_nstpt,   /*	number of time steps for reserve
                                       calculation */
                     int reserve_nstpx,   /*	number of space steps for
                                       reserve calculation */
                     int reserve_nstppsi, /*	number of time steps for reserve
                                       calculation */
                     int reserve_nstpvol, /*	number of space steps for
                                       reserve calculation */
                     double reserve_integ_mintime,

                     /*	Results */
                     double *switch_reserve, double *onetime_reserve,

                     CTS_IV fwd_iv_info, int save_extra_infos,
                     CTS_EXTRA_INFOS extra_infos) {

  Err err = NULL;

  CTS_MKT mkt;

  double lambda_init;
  int nbfactor_init;
  int i;

  CTS_ADI_ARG adi_arg_init;
  int free_adi_arg_lgm = 0, free_adi_arg_mc = 0;

  cts_adi_arg *adi_arg_lgm = NULL, *adi_arg_mc = NULL;

  int init_nsmilepar;
  double *init_smilepartime = NULL, *init_alpha = NULL, *init_rho = NULL,
         *init_ldaeps = NULL, *init_rho2 = NULL;

  double call_init, call_new, fwd_iv;
  double vol_shift_tot;

  LGMSV_CalibParams *local_lgm_sv_calib_params;

  /* Get the market */
  mkt = und->mkt;

  /* Get the initial model informations */
  lambda_init = und->model.dLambdaX;
  nbfactor_init = und->model.iOne2F;

  init_nsmilepar = und->model.iNbPWTime;
  init_smilepartime = calloc(init_nsmilepar, sizeof(double));
  init_alpha = calloc(init_nsmilepar, sizeof(double));
  init_rho = calloc(init_nsmilepar, sizeof(double));
  init_ldaeps = calloc(init_nsmilepar, sizeof(double));
  init_rho2 = calloc(init_nsmilepar, sizeof(double));

  local_lgm_sv_calib_params = calloc(1, sizeof(LGMSV_CalibParams));

  if (!local_lgm_sv_calib_params || !init_smilepartime || !init_alpha ||
      !init_rho || !init_ldaeps || !init_rho2) {
    err = "Memory allocation faillure in cts_calc_reserve";
    goto FREE_RETURN;
  }

  /* Copy Params */
  LGMSV_Copy_CalibParams(lgmsv_calib_params, local_lgm_sv_calib_params);

  for (i = 0; i < init_nsmilepar; i++) {
    init_smilepartime[i] = und->model.dPWTime[i];
    init_alpha[i] = und->model.dAlpha[i] / 2.0;
    init_ldaeps[i] = und->model.dInitLambdaEps / 2.0;
    init_rho[i] = und->model.dRho[i];

    if (und->model.iOne2F == 2) {
      init_rho2[i] = und->model.dRho2[i];
    } else {
      init_rho2[i] = und->model.dRho[i];
    }

    /* tunr off the one_Factor Rho */
    if (und->model.iOne2F == 2) {
      /* Remove the one factor rho since it is already converted */
      local_lgm_sv_calib_params->onefac_rho = 0;
    }
  }

  /* Numer TStar */
  cts_get_numer_tstar(mkt->today, cts_mdl, &numer_tstar);

  switch (reserve_method) {
  /* Initial method as of 08/2003 */
  /* **************************** */
  case 0: {
    /* Get and freeze the initial lambda */
    fix_lambda = 1;

    /* Freeze the smile parameters if needed */
    if (lgm_reserve || midat_reserve) {
      /* switch to one factor model */
      for (i = 0; i < nsmilepar; i++) {
        alphaeps[i] = 0.005;
        rhoeps[i] = 0.0;
        ldaeps[i] = 0.0;
      }

      fix_smile = 1;

      /* remove the steps on vol and phi */
      reserve_nstppsi = 1;
      reserve_nstpvol = 1;
    }

    /* First part: we get the one-time callable values */

    if (one_time_index == 0) {
      /* First Calculate the index of the one time Call */
      err = cts_find_one_time_index(und, cts_mdl, &one_time_index);

      if (err)
        goto FREE_RETURN;
    }

    /*	1 Step: get the one time callable infos */
    /*	*************************************** */

    free_adi_arg_lgm = 0;

    adi_arg_lgm = calloc(1, sizeof(cts_adi_arg));

    if (!adi_arg_lgm) {
      err = "Memory allocation faillure in cts_calc_reserve";
      goto FREE_RETURN;
    }

    if (lgm_reserve) {
      /* Recalibrate the model using LGM */

      free_LGMSV_model(&(und->model));

      err = cts_calib_und(
          mkt, 0, lambda_init, 1, 0.0, 0.0, 0.0, nsmilepar, smilepartime,
          alphaeps, rhoeps, ldaeps, NULL, tstar, NumerParams, cal_tenor,
          cal_ref, cal_freq, cal_basis, force_atm, max_std_long, max_std_short,
          vol_shift_long, vol_type_long, vol_shift_type_long, vol_shift_short,
          vol_type_short, vol_shift_type_short, lambda_shift, calib_strategy,
          fix_lambda, short_tenor, short_refrate, short_freq, short_basis,
          fix_smile, smile_calib_months, local_lgm_sv_calib_params, min_time,
          skip_last, min_fact, max_fact, use_jumps, prec, maxiter, keep_first,
          long_strike_flag, short_strike_flag, cts_mdl, und, 0, NULL);

      if (err)
        goto FREE_RETURN;

      err = cts_calc_model_fwd_iv(und, cts_mdl, calc_fwd_iv, adj_fee,
                                  NumerParams);

      if (err)
        goto FREE_RETURN;

      Convert_Tstar_model(&(und->model), numer_tstar);

      err = cts_fill_algo_arg(und, cts_mdl, 0, reserve_nstpt, reserve_nstppsi,
                              reserve_nstpx, reserve_nstpvol, 0, 0,
                              reserve_integ_mintime, 1, one_time_index,
                              adi_arg_lgm);

      free_adi_arg_lgm = 1;

      if (err)
        goto FREE_RETURN;

      err = cts_adjust_model_fwd_iv(und, cts_mkt, cts_mdl, for_fund,
                                    calc_fwd_iv, adj_fee, 0, 0, NULL);

      if (err)
        goto FREE_RETURN;

      /* Launch the adi */
      smessage(
          "Launching adi        , time steps requested: %d        , actual: %d",
          reserve_nstpt, adi_arg_lgm->nstp);

      err = cts_launch_algo(cts_mdl, und, adi_arg_lgm, &call,
                            &(onetime_call[0]), NULL);

      if (err)
        goto FREE_RETURN;

      adi_arg_init = adi_arg_lgm;
    } else {
      adi_arg_init = adi_arg;
    }

    /*	Step2: reserve calculation */
    /*	************************** */

    if (midat_reserve) {
      err = cts_calc_reserve_midat(
          cts_mkt, cts_mdl, und, adi_arg_init, call, onetime_call,
          recalib_european, NumerParams, tstar, numer_tstar, cal_tenor, cal_ref,
          cal_freq, cal_basis, force_atm, max_std_long, max_std_short,
          vol_shift_long, vol_type_long, vol_shift_type_long, vol_shift_short,
          vol_type_short, vol_shift_type_short, lambda_shift, calib_strategy,
          fix_lambda, short_tenor, short_refrate, short_freq, short_basis,
          fix_smile, smile_calib_months, local_lgm_sv_calib_params, min_time,
          skip_last, min_fact, max_fact, use_jumps, prec, maxiter, keep_first,
          long_strike_flag, short_strike_flag, lambda_init, lambda_reserve,
          onetime_vega, reserve_nstpt, reserve_nstpx, reserve_integ_mintime,
          switch_reserve, onetime_reserve, fwd_iv_info, save_extra_infos,
          extra_infos);

      if (err)
        goto FREE_RETURN;

      Convert_Tstar_model(&(und->model), und->model.dInitTStar);
    } else {
      cts_calc_reserve_timeswap(
          cts_mkt, cts_mdl, und, call, onetime_call, adi_arg_init->has_one_time,
          for_fund, nsmilepar, smilepartime, alphaeps, rhoeps, ldaeps, tstar,
          NumerParams, numer_tstar, cal_tenor, cal_ref, cal_freq, cal_basis,
          force_atm, max_std_long, max_std_short, vol_shift_long, vol_type_long,
          vol_shift_type_long, vol_shift_short, vol_type_short,
          vol_shift_type_short, lambda_shift, calib_strategy, fix_lambda,
          short_tenor, short_refrate, short_freq, short_basis, fix_smile,
          smile_calib_months, local_lgm_sv_calib_params, min_time, skip_last,
          min_fact, max_fact, use_jumps, prec, maxiter, keep_first,
          long_strike_flag, short_strike_flag, calc_fwd_iv, adj_fee,
          lambda_init, lambda_reserve, onetime_vega, euro_nb_iter,
          one_time_index, reserve_nstpt, reserve_nstppsi, reserve_nstpx,
          reserve_nstpvol, reserve_integ_mintime, switch_reserve,
          onetime_reserve, fwd_iv_info, save_extra_infos, extra_infos);

      if (err)
        goto FREE_RETURN;
    }

    if (save_extra_infos) {
      /* We save the latest term structure in Extra Infos */
      ConvertTS_LGMSV_to_LGM(und->model.iNbPWTime, und->model.dPWTime,
                             und->model.dSigma, und->model.dLambdaX,
                             und->model.dTStar);

      extra_infos->ncol = 2;
      extra_infos->nrow = und->model.iNbPWTime + 1;

      extra_infos->extra_infos =
          dmatrix(0, extra_infos->nrow - 1, 0, extra_infos->ncol - 1);

      if (!extra_infos->extra_infos) {
        err = "Memory allocation faillure in cts_caller";
        goto FREE_RETURN;
      }

      extra_infos->extra_infos[0][0] = und->model.dLambdaX;
      extra_infos->extra_infos[0][1] = und->model.dTau;

      for (i = 0; i < und->model.iNbPWTime; i++) {
        extra_infos->extra_infos[i + 1][0] =
            mkt->today + und->model.dPWTime[i] * DAYS_IN_YEAR;
        extra_infos->extra_infos[i + 1][1] = und->model.dSigma[i];
      }
    }

    break;
  }

  /* New reserve methodology as of 06/2004 */
  /* ************************************* */
  case 1: {
    free_adi_arg_mc = 0;

    adi_arg_mc = calloc(1, sizeof(cts_adi_arg));

    if (!adi_arg_mc) {
      err = "Memory allocation faillure in cts_calc_reserve";
      goto FREE_RETURN;
    }

    /* Switch Reserve = Pricing in the 2F model */

    if (fix_lambda == 0 && und->model.iOne2F == 2) {
      /* Initial price made in the 2F model */
      *switch_reserve = 0.0;
    } else {
      /* Initial call price */
      if (recalc_one_factor) {
        /* Recalculate one factor price */
        err = cts_fill_algo_arg(und, cts_mdl, 1, reserve_nstpt, reserve_nstppsi,
                                reserve_nstpx, reserve_nstpvol,
                                adi_arg->mc_mintime, adi_arg->npaths,
                                reserve_integ_mintime, 0, 0, adi_arg_mc);

        free_adi_arg_mc = 1;

        if (err)
          goto FREE_RETURN;

        err = cts_launch_algo(cts_mdl, und, adi_arg_mc, &call_init, NULL, NULL);

        if (err)
          goto FREE_RETURN;

        if (free_adi_arg_mc)
          cts_free_adi_arg(adi_arg_mc);
        free_adi_arg_mc = 0;
      } else {
        call_init = call;
      }

      /* 2F Price */
      free_LGMSV_model(&(und->model));

      /* We reserve against full calibrated model */
      err = cts_calib_und(
          mkt, 0, lambda_init, 2, lgm_alpha, lgm_gamma, lgm_rho, init_nsmilepar,
          init_smilepartime, init_alpha, init_rho, init_ldaeps, init_rho2,
          tstar, NumerParams, cal_tenor, cal_ref, cal_freq, cal_basis,
          force_atm, max_std_long, max_std_short, vol_shift_long, vol_type_long,
          vol_shift_type_long, vol_shift_short, vol_type_short,
          vol_shift_type_short, lambda_shift, calib_strategy, 0, short_tenor,
          short_refrate, short_freq, short_basis, 1, smile_calib_months,
          local_lgm_sv_calib_params, min_time, skip_last, min_fact, max_fact,
          use_jumps, prec, maxiter, keep_first, long_strike_flag,
          short_strike_flag, cts_mdl, und, 0, NULL);

      if (err)
        goto FREE_RETURN;

      Convert_Tstar_model(&(und->model), numer_tstar);

      err = cts_calc_model_fwd_iv(und, cts_mdl, calc_fwd_iv, adj_fee,
                                  NumerParams);

      if (err)
        goto FREE_RETURN;

      err = cts_fill_algo_arg(und, cts_mdl, 1, adi_arg->initnstp,
                              adi_arg->nstppsi, adi_arg->nstpx, adi_arg->nstpz,
                              adi_arg->mc_mintime, adi_arg->npaths,
                              reserve_integ_mintime, 0, 0, adi_arg_mc);

      free_adi_arg_mc = 1;

      if (err)
        goto FREE_RETURN;

      err = cts_adjust_model_fwd_iv(und, cts_mkt, cts_mdl, for_fund,
                                    calc_fwd_iv, adj_fee, 1, 0, NULL);

      if (err)
        goto FREE_RETURN;

      /* Launch the MC */
      err = cts_launch_algo(cts_mdl, und, adi_arg_mc, &call_new, NULL, &fwd_iv);

      if (err)
        goto FREE_RETURN;

      *switch_reserve = call_new - call_init;

      if (save_extra_infos && extra_infos) {
        extra_infos->algo_fwd_iv2 = fwd_iv;
      }

      if (free_adi_arg_mc)
        cts_free_adi_arg(adi_arg_mc);
      free_adi_arg_mc = 0;
    }

    /* One Time Reserve */
    /* We shift only the swaptions */
    if (fabs(onetime_vega) > 1.0E-08) {
      /* 2F Price */
      free_LGMSV_model(&(und->model));

      vol_shift_tot = vol_shift_long + onetime_vega;

      err = cts_calib_und(
          mkt, 0, lambda_init, nbfactor_init, lgm_alpha, lgm_gamma, lgm_rho,
          init_nsmilepar, init_smilepartime, init_alpha, init_rho, init_ldaeps,
          init_rho2, tstar, NumerParams, cal_tenor, cal_ref, cal_freq,
          cal_basis, force_atm, max_std_long, max_std_short, vol_shift_tot,
          LOGNORMAL_VOL, ADDITIVE, vol_shift_short, vol_type_short,
          vol_shift_type_short, lambda_shift, calib_strategy, fix_lambda,
          short_tenor, short_refrate, short_freq, short_basis, 1,
          smile_calib_months, local_lgm_sv_calib_params, min_time, skip_last,
          min_fact, max_fact, use_jumps, prec, maxiter, keep_first,
          long_strike_flag, short_strike_flag, cts_mdl, und, 0, NULL);

      if (err)
        goto FREE_RETURN;

      Convert_Tstar_model(&(und->model), numer_tstar);

      err = cts_calc_model_fwd_iv(und, cts_mdl, calc_fwd_iv, adj_fee,
                                  NumerParams);

      if (err)
        goto FREE_RETURN;

      err = cts_fill_algo_arg(
          und, cts_mdl, adi_arg->pde_or_mc, adi_arg->initnstp, adi_arg->nstppsi,
          adi_arg->nstpx, adi_arg->nstpz, adi_arg->mc_mintime, adi_arg->npaths,
          reserve_integ_mintime, 0, 0, adi_arg_mc);

      free_adi_arg_mc = 1;

      if (err)
        goto FREE_RETURN;

      err =
          cts_adjust_model_fwd_iv(und, cts_mkt, cts_mdl, for_fund, calc_fwd_iv,
                                  adj_fee, adi_arg->pde_or_mc, 0, NULL);

      if (err)
        goto FREE_RETURN;

      /* Launch the MC */
      err = cts_launch_algo(cts_mdl, und, adi_arg_mc, &call_new, NULL, NULL);

      if (err)
        goto FREE_RETURN;

      *onetime_reserve = call_new - call;
    } else {
      *onetime_reserve = 0.0;
    }
  }
  }

FREE_RETURN:

  if (adi_arg_lgm) {
    if (free_adi_arg_lgm)
      cts_free_adi_arg(adi_arg_lgm);
    free(adi_arg_lgm);
  }

  if (adi_arg_mc) {
    if (free_adi_arg_mc)
      cts_free_adi_arg(adi_arg_mc);
    free(adi_arg_mc);
  }

  if (local_lgm_sv_calib_params)
    free(local_lgm_sv_calib_params);

  if (init_smilepartime)
    free(init_smilepartime);
  if (init_alpha)
    free(init_alpha);
  if (init_rho)
    free(init_rho);
  if (init_ldaeps)
    free(init_ldaeps);
  if (init_rho2)
    free(init_rho2);

  return err;
}

void shift_LGMSV_model(long today, char *yc, LGMSV_MODEL model, double shift,
                       double log_norm) {
  int i;
  double short_rate;

  /* first transform into LGM */
  ConvertTS_LGMSV_to_LGM(model->iNbPWTime, model->dPWTime, model->dSigma,
                         model->dLambdaX, model->dTStar);
  if (log_norm) {
    for (i = 0; i < model->iNbPWTime; i++) {
      model->dSigma[i] += shift;

      if (model->dSigma[i] < LGMSV_MINVOL) {
        model->dSigma[i] = LGMSV_MINVOL;
      }
    }
  } else {
    for (i = 0; i < model->iNbPWTime; i++) {
      short_rate = swp_f_zr(today + model->dPWTime[i] * DAYS_IN_YEAR,
                            today + model->dPWTime[i] * DAYS_IN_YEAR + 1, yc);

      model->dSigma[i] += shift * model->dSigma[i] / short_rate;

      if (model->dSigma[i] < LGMSV_MINVOL) {
        model->dSigma[i] = LGMSV_MINVOL;
      }
    }
  }

  ConvertTS_LGM_to_LGMSV(model->iNbPWTime, model->dPWTime, model->dSigma,
                         model->dLambdaX, model->dTStar, 1, 0.0, 0.0, 0.0, NULL,
                         NULL);
}
