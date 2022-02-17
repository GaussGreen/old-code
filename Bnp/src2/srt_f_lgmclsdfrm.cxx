/* ===================================================================================
   FILENAME:      srt_f_lgmclsdfrm.cxx

   PURPOSE:       Compute swaption        , bond option        ,caps/floors
   prices in LGM        , via a closed form analytical solution.
   ===================================================================================
 */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lgmclsdfrm.h"
#include "utconst.h"

#define MAXITER 20
#define Version 1
/* if Version == 2 then we use Hermite interpolation with n = 25 (may be
   generic) Version == 1 then we use Hermite interpolation with n = 12
   (hardcoded) Version == 0 then we use the old version
*/

static double lgm_discount_bond_option(
    double fixing_time, double period_start_time, double period_end_time,
    double period_start_df, double period_end_df, double bond_strike,
    TermStruct *ts, SrtReceiverType pay_rec, SrtMdlDim mdl_dim);

static double d_function_onefac(int n, double d, double *coupon, double *df,
                                double *s_i, double bond_strike);

static double d_function_twofac(int n, double d, double g, double *coupon,
                                double *df, double *s_i_squared, double *S_i_3,
                                double *S_i_4, double bond_strike);

static Err onefac_exer_probs_under_QTi(double *s_i, double *coupon,
                                       double bond_strike, double *df, int n,
                                       SrtReceiverType rec_pay, double *prob);

static Err twofac_exer_probs_under_QTi_for_g(double *prob, double g,
                                             double *s_i_3, double *s_i_4,
                                             double *s_i_squared,
                                             double *coupon, double bond_strike,
                                             double *df, int n,
                                             SrtReceiverType rec_pay);

/* ----------------------------------------------------------------------------
        For the pricing of a cap/floor in LGM (with spreads already in the
   coupon)
   ----------------------------------------------------------------------------
 */

double srt_f_lgm_capfloor(
    TermStruct *ts, double *fixing_time,
    double *period_time,         /* period_time[0] is the strike payment time */
    double *df, double *payment, /* This is cvg * ( cash_fwd + spread ) */
    int num_caplets, SrtReceiverType rec_pay, SrtMdlDim mdl_dim) {
  double price;
  double caplet;
  int i;

  /*	We think of a cap as a sum of zero coupon bond options */
  price = 0.0;
  for (i = 0; i < num_caplets; i++) {
    /* Skip caplets that fixed today or before */
    if (fixing_time[i] > 0) {
      caplet = lgm_discount_bond_option(
          fixing_time[i], period_time[i], period_time[i + 1], df[i], df[i + 1],
          1.0 / (1.0 + payment[i + 1]), ts, rec_pay, mdl_dim);
      caplet *= (1 + payment[i + 1]);
      price += caplet;
    }
  }

  return price;

} /* END srt_f_lgm_capfloor (...) */

/* ----------------------------------------------------------------------------
 */

/* ==================================================================

    lgm_coupon_bond_option()

        prices a bond option or a swaption (option on bond at par) in LGM

        the inputs:
        n		     dates/coupons of underlying go from [0] to [n]
        bond_strike	 strike of swaption/bond option
        *ts		     term structure of volatility
        *fixing      time    time in years from today to fixing date
        *coupon      array of n coupons of underlying swap/bond        ,
                                 starting with coupon[1]
        *pay_time	 array of n times of corresponding coupons (in years ,
                                 from today);
                                 Index of first coupon is 1.  pay_time[0]
                                 is the maturity of the option and not used.
                                 (Today is 0.0        , hence there is no
                                 variable "today")

        *df	array of n+1 dfs for n+1 times in pay_time.
                        (indexed as above)

      pay_rec        rec means call on bond        , which is equivalent
                        to put on rates; receive fixed (receiver);
                        pay means put on bond; call on rates; pay fixed

        returns price.
   ======================================================================== */

double srt_f_lgm_coupon_bond_option(int n, double bond_strike, TermStruct *ts,
                                    double fixing_time, double *coupon,
                                    double *pay_time, double *df,
                                    SrtReceiverType pay_rec,
                                    SrtMdlDim mdl_dim) {
  int i, j;
  double *s_i, *s_i_squared, *s_i_3, *s_i_4, *cond_prob, *Reflected_cond_prob,
      *prob;
  double price;
  double intrinsic;
  int rp;
  Err err;

  /* Check if option has expired */
  if (fixing_time < 0.0)
    return 0.0;
  rp = (pay_rec == SRT_RECEIVER ? 1 : -1);

  /* Compute intrisic value */
  intrinsic = 0.0;
  for (i = 1; i <= n; i++)
    intrinsic += rp * coupon[i] * df[i];

  intrinsic -= rp * bond_strike * df[0];

  if (intrinsic < 0.0)
    intrinsic = 0.0;

  if (fixing_time == 0.0)
    return intrinsic;

  /* Memory allocation */
  s_i = srt_calloc(n + 1, sizeof(double));
  s_i_squared = srt_calloc(n + 1, sizeof(double));
  s_i_3 = srt_calloc(n + 1, sizeof(double));
  s_i_4 = srt_calloc(n + 1, sizeof(double));
  cond_prob = srt_calloc(n + 1, sizeof(double));
  Reflected_cond_prob = srt_calloc(n + 1, sizeof(double));
  prob = srt_calloc(n + 1, sizeof(double));
  memset(prob, 0, (n + 1) * sizeof(double));

  /* Compute s_i_squared        , cum vols for each cashflow */
  for (i = 1; i <= n; i++) {
    if (mdl_dim == ONE_FAC) {
      s_i[i] = sqrt(
          srt_f_lgm_cum_vol(ts, fixing_time, pay_time[0], pay_time[i], NULL));
    } else if (mdl_dim == TWO_FAC) {
      err = make_2f_lgm_cum_vols(ts, fixing_time, pay_time[0], pay_time[i],
                                 &s_i_3[i], &s_i_4[i], &s_i_squared[i]);
    }
  }

  if (mdl_dim == ONE_FAC) {
    err = onefac_exer_probs_under_QTi(s_i, coupon, bond_strike, df, n, pay_rec,
                                      prob);
  } else if (mdl_dim == TWO_FAC) {
/* Modif by Alan : Use Hermit point instead of The previous horrible linear
 * integration */
/* Version 2 : 25 pts Hermite interpolation
                                        Generic code with Hermite points
   generated by NR gauss_hermite */
#if Version == 2
    int NumHermitePoints = 25;
    double *gi, *wi;
    gi = (double *)malloc((NumHermitePoints + 1) * sizeof(double));
    wi = (double *)malloc((NumHermitePoints + 1) * sizeof(double));
    gauss_hermite(gi, wi, NumHermitePoints);
    for (i = 1; i <= NumHermitePoints; i++) {
      err = twofac_exer_probs_under_QTi_for_g(cond_prob, gi[i] * SQRT_TWO,
                                              s_i_3, s_i_4, s_i_squared, coupon,
                                              bond_strike, df, n, pay_rec);
      for (j = 0; j <= n; j++)
        prob[j] += INV_SQRT_PI * wi[i] * exp(-0.5 * s_i_3[j] * s_i_3[j]) *
                   cond_prob[j] * exp(-s_i_3[j] * gi[i] * SQRT_TWO);
    }
    free(gi);
    free(wi);

    /* Version 1 : Fast 12 pts Hermite interpolation */

#elif Version == 1
    double gi[6] = {0.444403001944139, 1.340375197151620, 2.259464451000790,
                    3.223709828770100, 4.271825847932280, 5.500901704467750};
    double wi[6] = {3.216643615128420E-01, 1.469670480453520E-01,
                    2.911668791236190E-02, 2.203380687533160E-03,
                    4.837184922590710E-05, 1.499927167637000E-07};
    for (i = 0; i <= 5; i++) {
      err = twofac_exer_probs_under_QTi_for_g(cond_prob, gi[i], s_i_3, s_i_4,
                                              s_i_squared, coupon, bond_strike,
                                              df, n, pay_rec);
      err = twofac_exer_probs_under_QTi_for_g(Reflected_cond_prob, -gi[i],
                                              s_i_3, s_i_4, s_i_squared, coupon,
                                              bond_strike, df, n, pay_rec);
      for (j = 0; j <= n; j++)
        prob[j] += wi[i] * exp(-0.5 * s_i_3[j] * s_i_3[j]) *
                   (cond_prob[j] * exp(-s_i_3[j] * gi[i]) +
                    Reflected_cond_prob[j] * exp(s_i_3[j] * gi[i]));
    }

/* Version 0 : default linear interpolation  form -4; +4 and 50 pts */
#else
    double length, delta, g;
    int seg_numb;
    length = 4.0;
    seg_numb = 50;
    delta = 2 * length / (double)seg_numb;
    g = -length;
    for (i = 1; i <= seg_numb; i++) {
      err = twofac_exer_probs_under_QTi_for_g(cond_prob, g, s_i_3, s_i_4,
                                              s_i_squared, coupon, bond_strike,
                                              df, n, pay_rec);
      for (j = 0; j <= n; j++)
        prob[j] += cond_prob[j] * gauss(g + s_i_3[j]) * delta;

      g += delta;
    }
#endif
  }

  /* Option price: (compare black scholes) */

  price = 0.0;
  for (i = 1; i <= n; i++)
    price += rp * coupon[i] * df[i] * prob[i];
  price -= rp * bond_strike * df[0] * prob[0];

  srt_free(s_i_squared);
  srt_free(s_i);
  srt_free(s_i_3);
  srt_free(s_i_4);
  srt_free(prob);
  srt_free(cond_prob);
  srt_free(Reflected_cond_prob);
  return price;

} /* END lgm_coupon_bond_option(...) */

/* ------------------------------------------------------------------------------
 */
/* Price an option on a discount bond (zero coupon)        , using one factor
   gauss markov model */

static double lgm_discount_bond_option(
    double fixing_time, double period_start_time, double period_end_time,
    double period_start_df, double period_end_df, double bond_strike,
    TermStruct *ts, SrtReceiverType pay_rec, SrtMdlDim mdl_dim) {

  double s, d;
  double price;
  int rp;
  double cv3, cv4;
  Err err;

  rp = (pay_rec == SRT_RECEIVER ? 1 : -1);

  if (mdl_dim == TWO_FAC) {
    err = make_2f_lgm_cum_vols(ts, fixing_time, period_start_time,
                               period_end_time, &cv3, &cv4, &s);
  } else if (mdl_dim == ONE_FAC) {
    s = srt_f_lgm_cum_vol(ts, fixing_time, period_start_time, period_end_time,
                          NULL);
  }
  s = sqrt(s);

  if (s == 0.0)
    return 0.0;

  d = 1.0 / s * log(period_end_df / (bond_strike * period_start_df)) - s / 2.0;

  price = rp * period_end_df * norm(rp * (d + s)) -
          rp * period_start_df * norm(rp * d) * bond_strike;

  return price;
}

/* ------------------------------------------------------------------------------
 */

/*  =======================================================================
        The function of d is the following

        Sum(i=1        ,n) { c_i B(0        ,Ti)exp(-S_i^2/2 - S_i*d) }
                                                - bond_strike * B(0        ,T0)

        (see The One Factor Gauss Markov Mode and the Two Factor Model in Grfn)

   =========================================================================  */

static double d_function_onefac(int n, double d, double *coupon, double *df,
                                double *s_i, double bond_strike) {
  double sum = 0;
  int i;

  for (i = 1; i <= n; i++) {
    sum += coupon[i] * df[i] * exp(-0.5 * s_i[i] * s_i[i] - d * s_i[i]);
  }
  sum -= bond_strike * df[0];

  return sum;
}

/*  =======================================================================
        compute a function of d (itself a function of a random number g )

        sum_1^n{ c_i B(t        ,Ti)exp(-S_i^2/2 - s_i_3*g - s_i_4*d) }
                                                - bond_strike * B(t        ,T)

        (for newton rhapson iteration in lgm_coupon_bond_option())
        =========================================================================
 */

static double d_function_twofac(int n, double d, double g, double *coupon,
                                double *df, double *s_i_squared, double *S_i_3,
                                double *S_i_4, double bond_strike) {
  double sum = 0;
  int i;

  for (i = 1; i <= n; i++) {
    sum += coupon[i] * df[i] *
           exp(-0.5 * s_i_squared[i] - S_i_3[i] * g - S_i_4[i] * d);
  }
  sum -= bond_strike * df[0];

  return sum;
}

/* =========================================================================

   Compute the exercise probabilities under each numeraire        , knowing the
   value of g (that correspond to the value of one of the factor):
        prob[0] => probability under B(...        ,T) prob measure  knowing g
        prob[1] => probability under B(...        ,T1) prob measure  knowing g
         ...
        prob[n] => probability under B(...        ,Tn) prob measure  knowing g

   This requires that we  compute d        , (a function of g) which is the
   unique solution of:

        sum(1->n){ c_i B(t        ,Ti)exp(- .5 * s_i_squared - s_i_3*g - s_i_4*d
   )} = bond_strike * B(t        ,T) (this is computed using  d_function() with
   Newton-Rhapson style iteration

   ========================================================================= */
static Err onefac_exer_probs_under_QTi(double *s_i, double *coupon,
                                       double bond_strike, double *df, int n,
                                       SrtReceiverType rec_pay, double *prob) {
  Err err = NULL;
  double rp;
  double nstop, d;
  double a[3], b[3];
  int count, i;

  rp = (rec_pay == SRT_RECEIVER ? 1 : -1);

  /** Initial guesses: **/

  nstop = 0.0;
  a[0] = 1;
  b[0] = d_function_onefac(n, a[0], coupon, df, s_i, bond_strike);
  a[1] = 1.01;
  b[1] = d_function_onefac(n, a[1], coupon, df, s_i, bond_strike);
  a[2] = 1.02;
  count = 0;

  /** Newton iterations **/

  while (nstop < 1.0 && count < MAXITER) {
    b[2] = d_function_onefac(n, a[2], coupon, df, s_i, bond_strike);
    newton(0.0, 3, a, b, &nstop);
    count++;
  }
  /* Corrective term due to the fixing - period start lag */
  d = a[2];

  prob[0] = norm(rp * d);
  for (i = 1; i <= n; i++) {
    prob[i] = norm(rp * (d + s_i[i]));
  }

  return err;
}

/* -------------------------------------------------------------------------- */

static Err twofac_exer_probs_under_QTi_for_g(double *prob, double g,
                                             double *s_i_3, double *s_i_4,
                                             double *s_i_squared,
                                             double *coupon, double bond_strike,
                                             double *df, int n,
                                             SrtReceiverType rec_pay) {
  Err err = NULL;
  double rp;
  double nstop, d;
  double a[3], b[3];
  int count, i;
  rp = (rec_pay == SRT_RECEIVER ? 1 : -1);

  /** initial guesses: **/

  nstop = 0.0;
  a[0] = 1;
  b[0] = d_function_twofac(n, a[0], g, coupon, df, s_i_squared, s_i_3, s_i_4,
                           bond_strike);
  a[1] = 1.01;
  b[1] = d_function_twofac(n, a[1], g, coupon, df, s_i_squared, s_i_3, s_i_4,
                           bond_strike);
  a[2] = 1.02;
  count = 0;

  /** iterations **/

  while (nstop < 1.0 && count < MAXITER) {
    b[2] = d_function_twofac(n, a[2], g, coupon, df, s_i_squared, s_i_3, s_i_4,
                             bond_strike);
    newton(0.0, 2.0, a, b, &nstop);
    count++;
    /* Sets an arbitrary bound to prevent overflow */
    if (a[2] > fabs(30 * s_i_4[n] / s_i_3[n]))
      a[2] = fabs(30 * s_i_4[n] / s_i_3[n]);
    else if (a[2] < -fabs(30 * s_i_4[n] / s_i_3[n]))
      a[2] = -fabs(30 * s_i_4[n] / s_i_3[n]);
  }

  d = a[2];

  prob[0] = norm(rp * d);
  for (i = 1; i <= n; i++)
    prob[i] = norm(rp * (d + s_i_4[i]));

  return err;
}

/* -----------------------------------------------------------------------------
 */
/*  NEW MODEL NEW LGM
/*
------------------------------------------------------------------------------
*/
/* Price an option on a discount bond (zero coupon)        , using one factor
   gauss markov model */

static double newlgm_discount_bond_option(
    double fixing_time, double period_start_time, double period_end_time,
    double period_start_df, double period_end_df, double bond_strike,
    TermStruct *ts, SrtReceiverType pay_rec, SrtMdlDim mdl_dim) {

  double s, d;
  double price, vol;
  int rp;

  rp = (pay_rec == SRT_RECEIVER ? 1 : -1);

  s = find_struct_interp(period_end_time, H, ts);
  s -= find_struct_interp(period_start_time, H, ts);
  vol = find_struct_interp(fixing_time, G, ts);
  vol = sqrt(vol);

  s *= vol;

  if (vol == 0.0)
    return 0.0;

  d = 1.0 / s * log(period_end_df / (bond_strike * period_start_df)) - s / 2.0;

  price = rp * period_end_df * norm(rp * (d + s)) -
          rp * period_start_df * norm(rp * d) * bond_strike;

  return price;
}

double srt_f_newlgm_capfloor(
    TermStruct *ts, double *fixing_time,
    double *period_time,         /* period_time[0] is the strike payment time */
    double *df, double *payment, /* This is cvg * ( cash_fwd + spread ) */
    int num_caplets, SrtReceiverType rec_pay, SrtMdlDim mdl_dim) {
  double price;
  double caplet;
  int i;

  /*	We think of a cap as a sum of zero coupon bond options */
  price = 0.0;
  for (i = 0; i < num_caplets; i++) {
    /* Skip caplets that fixed today or before */
    if (fixing_time[i] > 0) {
      caplet = newlgm_discount_bond_option(
          fixing_time[i], period_time[i], period_time[i + 1], df[i], df[i + 1],
          1.0 / (1.0 + payment[i + 1]), ts, rec_pay, mdl_dim);
      caplet *= (1 + payment[i + 1]);
      price += caplet;
    }
  }

  return price;

} /* END srt_f_newlgm_capfloor (...) */

/* ----------------------------------------------------------------------------
 */

/* ==================================================================

    newlgm_coupon_bond_option()

        prices a bond option or a swaption (option on bond at par) in LGM

        the inputs:
        n		     dates/coupons of underlying go from [0] to [n]
        bond_strike	 strike of swaption/bond option
        *ts		     term structure of volatility
        *fixing      time    time in years from today to fixing date
        *coupon      array of n coupons of underlying swap/bond        ,
                                 starting with coupon[1]
        *pay_time	 array of n times of corresponding coupons (in years ,
                                 from today);
                                 Index of first coupon is 1.  pay_time[0]
                                 is the maturity of the option and not used.
                                 (Today is 0.0        , hence there is no
                                 variable "today")

        *df	array of n+1 dfs for n+1 times in pay_time.
                        (indexed as above)

      pay_rec        rec means call on bond        , which is equivalent
                        to put on rates; receive fixed (receiver);
                        pay means put on bond; call on rates; pay fixed

        returns price.
   ======================================================================== */

double srt_f_newlgm_coupon_bond_option(int n, double bond_strike,
                                       TermStruct *ts, double fixing_time,
                                       double *coupon, double *pay_time,
                                       double *df, SrtReceiverType pay_rec,
                                       SrtMdlDim mdl_dim) {
  int i;
  double *s_i, *prob;
  double price, s;
  double vol;
  double intrinsic;
  int rp;
  Err err;

  /* Check if option has expired */
  if (fixing_time < 0.0)
    return 0.0;
  rp = (pay_rec == SRT_RECEIVER ? 1 : -1);

  /* Compute intrisic value */
  intrinsic = 0.0;
  for (i = 1; i <= n; i++)
    intrinsic += rp * coupon[i] * df[i];

  intrinsic -= rp * bond_strike * df[0];

  if (intrinsic < 0.0)
    intrinsic = 0.0;

  if (fixing_time == 0.0)
    return intrinsic;

  /* Memory allocation */
  s_i = srt_calloc(n + 1, sizeof(double));
  prob = srt_calloc(n + 1, sizeof(double));
  memset(prob, 0, (n + 1) * sizeof(double));

  /* Compute s_i_squared        , cum vols for each cashflow */
  s = find_struct_interp(pay_time[0], H, ts);
  vol = find_struct_interp(fixing_time, G, ts);
  vol = sqrt(vol);

  for (i = 1; i <= n; i++) {
    s_i[i] = find_struct_interp(pay_time[i], H, ts);
    s_i[i] -= s;
    s_i[i] *= vol;
  }

  err = onefac_exer_probs_under_QTi(s_i, coupon, bond_strike, df, n, pay_rec,
                                    prob);

  /* Option price: (compare black scholes) */

  price = 0.0;
  for (i = 1; i <= n; i++)
    price += rp * coupon[i] * df[i] * prob[i];
  price -= rp * bond_strike * df[0] * prob[0];

  srt_free(s_i);
  srt_free(prob);

  return price;

} /* END newlgm_coupon_bond_option(...) */
