
#include "fxsabradi.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"
#include "swp_h_swap_pricing.h"

#define MAX_CPN 600
#define MAX_INST 600
#define NUM_HERMITE 6
#define NITER 10
#define BUMP_LAM_LM 5.0E-04
#define MIN_CALIB_TIME 0.03

/*	Maximum decreasing factor allowed on variance
                ex: the following structure
                                1jan2003	1.00%
                                1jan2004	0.25%
                will be refused and replaced by the following
                                1jan2003	1.00%
                                1jan2004	0.50%
                even if the second option is not perfectly matched */

#define MAX_FACT 0.25 /* 0.50 ^ 2 */

/*	Set parameter defaults */
void cpd_calib_set_default_param(CPD_DIAG_CALIB_PARAM param) {
  param->fx_vol_shift = 0.0;
  param->vol_shift = 0.0;
  param->vol_type = LOGNORMAL_VOL;
  param->shift_type = MULTIPLICATIVE;

  param->strike_type = 2;

  param->lambda_shift = 0.0;
  param->lambda_min = -0.3;
  param->lambda_max = 0.3;

  param->transform_vol = 0;
  param->max_std = 1.0;
  param->min_time = 1.0;
  param->skip_last = 0;
  param->min_calib_time = MIN_CALIB_TIME;
  param->keep_first = 1;
  param->min_fact = MAX_FACT;
  param->max_fact = MAX_FACT / 4.0;
  param->use_jumps = 0;
  param->keep_first = 0;

  param->smile_vol_shift = 0.0;
  param->smile_vol_type = LOGNORMAL_VOL;
  param->smile_shift_type = MULTIPLICATIVE;
  param->smile_strike_type = 5;

  param->nb_iter_max = 7;
  param->precision = 0.00001;
  param->vega_prec = 0;

  param->smile_use_jumps = 1;
  param->smile_nb_iter_max = 7;
  param->smile_precision = 0.0001;
}

void diag_calib_lm_params_set_default_param(DIAG_CALIB_LM_PARAMS param) {
  param->nb_iter = 10;
  param->use_moment = 1;
  param->vega_weight = 0;
  param->freq_short = 1;
  param->shift_freq = 0;
  param->nb_moment = 1;
  param->break_moment = 0;
  param->precision = 0.0001;

  param->use_new = 0;
}

Err diagcalib_interp_voltype(const char *constStr,
                             DIAGCALIB_VOLTYPE *vol_type) {
  char str[30 + 1];
  if (!constStr)
    return "Empty string in diagcalib_interp_voltype";

  strncpy(str, constStr, 30);
  str[30] = '\0';
  strupper(str);
  strip_white_space(str);

  if (!strcmp(str, "NORMAL")) {
    *vol_type = NORMAL_VOL;
    return 0;
  }
  if (!strcmp(str, "NORM")) {
    *vol_type = NORMAL_VOL;
    return 0;
  }
  if (!strcmp(str, "LOGNORMAL")) {
    *vol_type = LOGNORMAL_VOL;
    return 0;
  }
  if (!strcmp(str, "LOG")) {
    *vol_type = LOGNORMAL_VOL;
    return 0;
  }
  if (!strcmp(str, "LGM")) {
    *vol_type = LGM_VOL;
    return 0;
  }
  if (!strcmp(str, "CAL")) {
    *vol_type = LGM_VOL;
    return 0;
  }

  if (!strcmp(str, "0")) {
    *vol_type = NORMAL_VOL;
    return 0;
  }
  if (!strcmp(str, "1")) {
    *vol_type = LOGNORMAL_VOL;
    return 0;
  }
  if (!strcmp(str, "2")) {
    *vol_type = LGM_VOL;
    return 0;
  }

  return serror("unknown diagcalib_interp_voltype. %s", str);
}

Err diagcalib_interp_shifttype(const char *constStr,
                               DIAGCALIB_SHIFTTYPE *shift_type) {
  char str[30 + 1];
  if (!constStr)
    return "Empty string in diagcalib_interp_shifttype";

  strncpy(str, constStr, 30);
  str[30] = '\0';
  strupper(str);
  strip_white_space(str);

  if (!strcmp(str, "ADDITIVE")) {
    *shift_type = ADDITIVE;
    return 0;
  }
  if (!strcmp(str, "ADD")) {
    *shift_type = ADDITIVE;
    return 0;
  }
  if (!strcmp(str, "MULTIPLICATIVE")) {
    *shift_type = MULTIPLICATIVE;
    return 0;
  }
  if (!strcmp(str, "MULT")) {
    *shift_type = MULTIPLICATIVE;
    return 0;
  }

  if (!strcmp(str, "0")) {
    *shift_type = ADDITIVE;
    return 0;
  }
  if (!strcmp(str, "1")) {
    *shift_type = MULTIPLICATIVE;
    return 0;
  }

  return serror("unknown diagcalib_interp_shifttype. %s", str);
}

/*	Static functions */

/*	Just a safer Gaussian */
double static_lgmsafenorm(double z) {
  if (z < -10) {
    z = -10.0;
  } else if (z > 10) {
    z = 10.0;
  }

  return norm_accurate(z);
  //	return norm (z);
}

/*	In order to help static_lgmystar: PV of coupons when reconstruction is
                df (t        , Ti) = df (0        , t        , Ti) * exp (a[i] +
   b[i] * Gaussian) and Gaussian is y */
static double static_lgmiv(double y,     /*	The y */
                           int ncpn,     /*	Num coupons */
                           double cpn[], /*	Discounted coupons */
                           double a[],   /*	The ai */
                           double b[])   /*	The bi */
{
  int i;
  double val;

  val = 0.0;
  for (i = 0; i < ncpn; i++) {
    val += cpn[i] * exp(a[i] + b[i] * y);
  }

  return val;
}

#define YPREC 1.0e-08
/*	Solve y / static_lgmiv (y) = 0 with Newton */
static double static_solvelgmy(
    int ncpn,     /*	Num coupons */
    double cpn[], /*	Discounted coupons */
    double a[],   /*	The ai */
    double b[],   /*	The bi */
    int *dir)     /*	Output        , whether iv is increasing in y
                                          1: decreasing        , -1: increasing */
{
  double y1, y2, iv1, iv2, temp;
  int it;

  y1 = 0.0;

  iv1 = static_lgmiv(y1, ncpn, cpn, a, b);

  y2 = 1.0e-02;

  iv2 = static_lgmiv(y2, ncpn, cpn, a, b);

  if (iv2 < iv1) {
    *dir = 1;
  } else {
    *dir = -1;
  }

  if (fabs(iv1) < YPREC || fabs(iv2) < YPREC) {

    if (fabs(iv1) < fabs(iv2)) {
      return y1;
    } else {
      return y2;
    }
  }

  for (it = 1; it < 25; it++) {
    if (fabs(y2 - y1) < YPREC || fabs(iv2 - iv1) < YPREC)
      break;

    temp = y2;
    y2 -= iv2 * (y2 - y1) / (iv2 - iv1);
    y1 = temp;
    iv1 = iv2;

    iv2 = static_lgmiv(y2, ncpn, cpn, a, b);

    if (fabs(iv2) < YPREC)
      break;
  }

  if (fabs(iv2) > YPREC) {
    /*	We should return an error here */
  }

  return y2;
}

/*	Solve y* in the model: 1F case */
static double static_lgmystar1F(
    int ncpn,     /*	Num coupons */
    double cpn[], /*	Discounted coupons */
    /*	G at coupon dates and zeta and G at exercise date        , in order to
       define a and b */
    double cpn_G[], double ex_zeta, double ex_sqzeta, double ex_G,
    int *dir) /*	Output        , whether iv is increasing in y
                                      1: decreasing        , -1: increasing */
{
  int i;
  double a[MAX_CPN], b[MAX_CPN];

  /*	Set a        , b */
  for (i = 0; i < ncpn; i++) {
    a[i] = -0.5 * ex_zeta * (cpn_G[i] - ex_G) * (cpn_G[i] - ex_G);
    b[i] = -(cpn_G[i] - ex_G) * ex_sqzeta;
  }

  /*	Solver */
  return static_solvelgmy(ncpn, cpn, a, b, dir);
}

/*	Solve y* in the model: LGM Stoch Vol case */
double static_lgmystar_stochvol(
    int ncpn,     /*	Num coupons */
    double cpn[], /*	Discounted coupons */
    double beta[], double psi,
    int *dir) /*	Output        , whether iv is increasing in y
                                      1: decreasing        , -1: increasing */
{
  static int i;
  static double a[MAX_CPN], b[MAX_CPN];

  /*	Set a        , b */
  for (i = 0; i < ncpn; i++) {
    a[i] = -0.5 * psi * beta[i] * beta[i];
    b[i] = -beta[i];
  }

  /*	Solver */
  return static_solvelgmy(ncpn, cpn, a, b, dir);
}

/*	Precalculate exponential factors in the IV        , and decide what
   variable to use in integration        , and which one to use in y* solving */
static void static_lgm2Fcalcexpfact(
    int ncpn, double cpn[], double cpn_G1[], double cpn_G2[], double ex_zeta1,
    double sqz1, double ex_zeta2, double sqz2, double ex_zeta12, double c1,
    double c2, double ex_G1, double ex_G2, double cstfact[], double n1fact[],
    double n2fact[],
    int *highest) /*	Index of the variable to be used for solving
                                          - the one with the highest sensitivity
                   */
{
  int i;
  double s1 = 0.0, s2 = 0.0;
  double half_z1 = 0.5 * ex_zeta1, half_z2 = 0.5 * ex_zeta2;
  double c1sqz2 = c1 * sqz2, c2sqz2 = c2 * sqz2;

  for (i = 0; i < ncpn; i++) {
    cstfact[i] = half_z1 * (cpn_G1[i] - ex_G1) * (cpn_G1[i] - ex_G1);
    cstfact[i] += half_z2 * (cpn_G2[i] - ex_G2) * (cpn_G2[i] - ex_G2);
    cstfact[i] += ex_zeta12 * (cpn_G1[i] - ex_G1) * (cpn_G2[i] - ex_G2);

    n1fact[i] = (cpn_G1[i] - ex_G1) * sqz1 + (cpn_G2[i] - ex_G2) * c1sqz2;

    n2fact[i] = (cpn_G2[i] - ex_G2) * c2sqz2;

    s1 -= cpn[i] * n1fact[i] * exp(-cstfact[i]);
    s2 -= cpn[i] * n2fact[i] * exp(-cstfact[i]);
  }

  if (fabs(s1) > fabs(s2)) {
    *highest = 1;
  } else {
    *highest = 2;
  }
}

/*	Solve y* in the model: 2F case */
static double static_lgmystar2F(
    int ncpn,     /*	Num coupons */
    double cpn[], /*	Discounted coupons */
                  /*	As output from static_lgm2Fcalcexpfact */
    double cstfact[], double n1fact[], double n2fact[], double highest,
    /*	Condition on Intergral variable = n */
    double n,
    int *dir) /*	Output        , whether iv is increasing in y
                                      1: decreasing        , -1: increasing */
{
  int i;

  double a[MAX_CPN], b[MAX_CPN];

  for (i = 0; i < ncpn; i++) {
    a[i] = -cstfact[i];

    if (highest == 1) {
      a[i] -= n2fact[i] * n;
      b[i] = -n1fact[i];
    } else {
      a[i] -= n1fact[i] * n;
      b[i] = -n2fact[i];
    }
  }

  /*	Solver */
  return static_solvelgmy(ncpn, cpn, a, b, dir);
}

/*	Calculate zeta1 from sigma        , lambda in the 2F model */
static void
static_lgmcalczeta1(int nsig,       /*	Number of Sigmas */
                    double sig_t[], /*	Sigma Times */
                    double sig[],   /*	Sigmas */
                    /*	Lambda        , Alpha        , Beta        , Rho */
                    double lambda, double alpha, double gamma, double rho,
                    /*	Output */
                    int nzeta, double zeta_t[], /*	Zeta times */
                    double zeta[])              /*	Output */
{
  int i, j;
  double t0;
  double z, dz;
  double zeta_at_sig_t[MAX_CPN];

  t0 = 0.0;
  z = 0.0;

  for (i = 0; i < nsig; i++) {
    dz = sig[i] * sig[i] * (exp(2 * lambda * sig_t[i]) - exp(2 * lambda * t0)) /
         2 / lambda;
    t0 = sig_t[i];
    z += dz;
    zeta_at_sig_t[i] = z;
  }

  i = 0;
  for (j = 0; j < nzeta; j++) {
    while (i < nsig && sig_t[i] < zeta_t[j]) {
      i++;
    }

    if (i == 0) {
      t0 = 0.0;
      z = 0.0;
    } else {
      t0 = sig_t[i - 1];
      z = zeta_at_sig_t[i - 1];
    }

    dz = sig[i] * sig[i] *
         (exp(2 * lambda * zeta_t[j]) - exp(2 * lambda * t0)) / 2 / lambda;
    z += dz;
    zeta_at_sig_t[i] = z;
  }
}

void static_interpolate_zeta(int nb_old_zeta, double old_zeta_time[],
                             double old_zeta[],

                             int nlam, double lam_time[], double lam[],

                             int nb_new_zeta, double new_zeta_time[],
                             double new_zeta[]) {
  double t1, t2, zeta1, zeta2, expfact1, expfact2;
  int i, last_index;
  double ratio;
  double last_time, last_zeta, last_vol2;

  i = 0;
  last_index = 0;

  /* Interpolation part */
  while (i < nb_new_zeta &&
         new_zeta_time[i] <= old_zeta_time[nb_old_zeta - 1]) {
    while (last_index < nb_old_zeta &&
           new_zeta_time[i] > old_zeta_time[last_index]) {
      last_index++;
    }

    if (old_zeta_time[last_index] - new_zeta_time[i] < 1.0E-08) {
      new_zeta[i] = old_zeta[last_index];
    } else {
      if (last_index > 0) {
        zeta1 = old_zeta[last_index - 1];
        t1 = old_zeta_time[last_index - 1];
      } else {
        zeta1 = 0.0;
        t1 = 0.0;
      }

      zeta2 = old_zeta[last_index];
      t2 = old_zeta_time[last_index];

      expfact1 = static_lgmcalcexpfact_tauts(t1, t2, nlam, lam_time, lam);
      expfact2 = static_lgmcalcexpfact_tauts(t1, new_zeta_time[i], nlam,
                                             lam_time, lam);

      ratio = expfact2 / expfact1;
      new_zeta[i] = ratio * zeta2 + (1.0 - ratio) * zeta1;
    }

    i++;
  }

  /* extrapolation for the others */
  if (i < nb_new_zeta) {
    /* calculation of the last vol */
    if (nb_old_zeta == 1) {
      last_zeta = 0.0;
      last_time = 0.0;
    } else {
      last_zeta = old_zeta[nb_old_zeta - 2];
      last_time = old_zeta_time[nb_old_zeta - 2];
    }

    expfact1 = static_lgmcalcexpfact_tauts(
        last_time, old_zeta_time[nb_old_zeta - 1], nlam, lam_time, lam);
    last_vol2 = (old_zeta[nb_old_zeta - 1] - last_zeta) / expfact1;

    last_zeta = old_zeta[nb_old_zeta - 1];
    last_time = old_zeta_time[nb_old_zeta - 1];

    while (i < nb_new_zeta) {
      expfact1 = static_lgmcalcexpfact_tauts(last_time, new_zeta_time[i], nlam,
                                             lam_time, lam);
      new_zeta[i] = expfact1 * last_vol2 + last_zeta;
      i++;
    }
  }
}

/*	Calculate exponential factor */
/*	Tau-TS enabled version */
double static_lgmcalcexpfact_tauts(double T1, double T2, int nlam,
                                   double lam_time[], double lam[]) {
  int i;
  double t0, ans, li;

  /*	Intergate lambda up to T1 */
  li = 0.0;
  i = 0;
  t0 = 0.0;
  while (i < nlam && lam_time[i] <= T1) {
    li += (lam_time[i] - t0) * lam[i];
    t0 = lam_time[i];
    i++;
  }
  if (i >= nlam)
    i--;
  li += (T1 - t0) * lam[i];

  /*	Intergate exp fact from T1 to T2 */
  ans = 0.0;
  t0 = T1;
  while (i < nlam && lam_time[i] <= T2) {
    if (fabs(lam[i]) > 1.0E-10) {
      ans += exp(2.0 * li) * (exp(2.0 * lam[i] * (lam_time[i] - t0)) - 1.0) /
             2.0 / lam[i];
    } else {
      ans += exp(2.0 * li) * (lam_time[i] - t0);
    }

    li += (lam_time[i] - t0) * lam[i];
    t0 = lam_time[i];
    i++;
  }
  if (i >= nlam)
    i--;

  if (fabs(lam[i]) > 1.0E-10) {
    ans += exp(2.0 * li) * (exp(2.0 * lam[i] * (T2 - t0)) - 1.0) / 2.0 / lam[i];
  } else {
    ans += exp(2.0 * li) * (T2 - t0);
  }

  return ans;
}

double export_lgmcalcexpfact_tauts(double T1, double T2, int nlam,
                                   double lam_time[], double lam[],
                                   double gamma) {
  int i;
  double exp_fact;

  for (i = 0; i < nlam; i++) {
    lam[i] += gamma;
  }

  exp_fact = static_lgmcalcexpfact_tauts(T1, T2, nlam, lam_time, lam);

  for (i = 0; i < nlam; i++) {
    lam[i] -= gamma;
  }

  return exp_fact;
}

/*	Calculate zeta1 from sigma        , lambda in the 2F model */
/*	Tau-TS enabled version */
static void static_lgmcalczeta1_tauts(
    int nsig,          /*	Number of Sigmas */
    double sig_time[], /*	Sigma Times */
    double sig[],      /*	Sigmas */
    /*	Lambda        , Alpha        , Beta        , Rho */
    int nlam, double lam_time[], double lam[], /*	Lambdas */
    double alpha, double gamma, double rho,
    /*	Output */
    int nzeta, double zeta_t[], /*	Zeta times */
    double zeta[])              /*	Output */
{
  int i, j;
  double t0;
  double z, dz;
  double exp_fact;
  double zeta_at_ts_t[MAX_CPN];

  t0 = 0.0;
  z = 0.0;

  for (i = 0; i < nsig; i++) {
    exp_fact =
        static_lgmcalcexpfact_tauts(t0, sig_time[i], nlam, lam_time, lam);
    dz = sig[i] * sig[i] * exp_fact;
    t0 = sig_time[i];
    z += dz;
    zeta_at_ts_t[i] = z;
  }

  /* Interpolation of zeta */
  i = 0;
  for (j = 0; j < nzeta; j++) {
    while (i < nsig && sig_time[i] < zeta_t[j]) {
      i++;
    }

    if (i < nsig) {
      if (i == 0) {
        t0 = 0.0;
        z = 0.0;
      } else {
        t0 = sig_time[i - 1];
        z = zeta_at_ts_t[i - 1];
      }

      exp_fact =
          static_lgmcalcexpfact_tauts(t0, zeta_t[j], nlam, lam_time, lam);
      dz = sig[i] * sig[i] * exp_fact;
      z += dz;
    } else {
      t0 = sig_time[nsig - 1];
      z = zeta_at_ts_t[nsig - 1];

      exp_fact =
          static_lgmcalcexpfact_tauts(t0, zeta_t[j], nlam, lam_time, lam);
      dz = sig[nsig - 1] * sig[nsig - 1] * exp_fact;
      z += dz;
    }

    zeta[j] = z;
  }
}

/*	Calculate zeta2 and zeta12 from zeta1        , lambda        , alpha ,
 * gamma and rho in the 2F model */
static void static_lgmcalczeta2zeta12(
    int n,          /*	Number of dates */
    double t[],     /*	Times */
    double zeta1[], /*	Zeta1 */
    /*	Lambda        , Alpha        , Beta        , Rho */
    double lambda, double alpha, double gamma, double rho,
    /*	Output */
    double zeta2[], double zeta12[]) {
  int i;
  double t0;
  double z1, z2, z12, dz1, dz2, dz12;
  double l1 = lambda, l2 = lambda + gamma;

  t0 = 0.0;
  z1 = z2 = z12 = 0.0;

  for (i = 0; i < n; i++) {
    dz1 = zeta1[i] - z1;
    dz2 = dz1 * alpha * alpha *
          ((exp(2 * l2 * t[i]) - exp(2 * l2 * t0)) / 2 / l2) /
          ((exp(2 * l1 * t[i]) - exp(2 * l1 * t0)) / 2 / l1);
    zeta2[i] = z2 + dz2;
    dz12 = dz1 * alpha * rho *
           ((exp((l1 + l2) * t[i]) - exp((l1 + l2) * t0)) / (l1 + l2)) /
           ((exp(2 * l1 * t[i]) - exp(2 * l1 * t0)) / 2 / l1);
    zeta12[i] = z12 + dz12;

    t0 = t[i];
    z1 = zeta1[i];
    z2 = zeta2[i];
    z12 = zeta12[i];
  }
}

void export_lgmcalczeta1_tauts(
    int nsig,          /*	Number of Sigmas */
    double sig_time[], /*	Sigma Times */
    double sig[],      /*	Sigmas */
    /*	Lambda        , Alpha        , Beta        , Rho */
    int nlam, double lam_time[], double lam[], /*	Lambdas */
    double alpha, double gamma, double rho,
    /*	Output */
    int nzeta, double zeta_t[], /*	Zeta times */
    double zeta[])              /*	Output */
{
  static_lgmcalczeta1_tauts(nsig, sig_time, sig, nlam, lam_time, lam, alpha,
                            gamma, rho, nzeta, zeta_t, zeta);
}

/*	Calculate zeta2 and zeta12 from zeta1        , lambda        , alpha ,
 * gamma and rho in the 2F model */
void export_lgmcalczeta2zeta12(
    int n,          /*	Number of dates */
    double t[],     /*	Times */
    double zeta1[], /*	Zeta1 */
    /*	Lambda        , Alpha        , Beta        , Rho */
    double lambda, double alpha, double gamma, double rho,
    /*	Output */
    double zeta2[], double zeta12[]) {
  static_lgmcalczeta2zeta12(n, t, zeta1, lambda, alpha, gamma, rho, zeta2,
                            zeta12);
}

/*	Calculate zeta2 and zeta12 from zeta1        , lambda        , alpha ,
 * gamma and rho in the 2F model */
void export_lgmcalczeta2zeta12_ts(
    int n,          /*	Number of dates */
    double t[],     /*	Times */
    double zeta1[], /*	Zeta1 */
    /*	Lambda        , Alpha        , Beta        , Rho */
    int nlambda, double lambda_time[], double lambda[], double alpha,
    double gamma, double rho,
    /*	Output */
    double zeta2[], double zeta12[]) {
  int i;
  double t0;
  double z1, z2, z12, dz1, dz2, dz12;

  /// previously
  /// double		l1 = lambda        , l2 = lambda + gamma;

  double dAvgl1_ti, dAvgl1_t0, dAvgl2_ti, dAvgl2_t0, dAvgl1, dAvgl2;

  t0 = 0.0;
  z1 = z2 = z12 = 0.0;

  for (i = 0; i < n; i++) {
    dz1 = zeta1[i] - z1;

    dAvgl1_ti = _average_lambda_(t[i], nlambda, lambda_time, lambda);
    dAvgl1_t0 = _average_lambda_(t0, nlambda, lambda_time, lambda);
    dAvgl1 = (dAvgl1_ti * t[i] - dAvgl1_t0 * t0) / (t[i] - t0);

    dAvgl2_ti = dAvgl1_ti + gamma;
    dAvgl2_t0 = dAvgl1_t0 + gamma;
    dAvgl2 = dAvgl1 + gamma;

    ///// previously
    // dz2 = dz1 * alpha * alpha
    //	* ((exp (2 * l2 * t[i]) - exp (2 * l2 * t0)) / 2 / l2)
    //	/ ((exp (2 * l1 * t[i]) - exp (2 * l1 * t0)) / 2 / l1);
    dz2 = dz1 * alpha * alpha *
          ((exp(2 * dAvgl2_ti * t[i]) - exp(2 * dAvgl2_t0 * t0)) / 2 / dAvgl2) /
          ((exp(2 * dAvgl1_ti * t[i]) - exp(2 * dAvgl1_t0 * t0)) / 2 / dAvgl1);

    zeta2[i] = z2 + dz2;

    ///// previously
    // dz12 = dz1 * alpha * rho
    //	* ((exp ((l1 + l2) * t[i]) - exp ((l1 + l2) * t0)) / (l1 + l2))
    //	/ ((exp (2 * l1 * t[i]) - exp (2 * l1 * t0)) / 2 / l1);
    dz12 = dz1 * alpha * rho *
           ((exp(dAvgl1_ti * t[i] + dAvgl2_ti * t[i]) -
             exp((dAvgl1_t0 + dAvgl2_t0) * t0)) /
            (dAvgl1 + dAvgl2)) /
           ((exp(2 * dAvgl1_ti * t[i]) - exp(2 * dAvgl1_t0 * t0)) / 2 / dAvgl1);

    zeta12[i] = z12 + dz12;

    t0 = t[i];
    z1 = zeta1[i];
    z2 = zeta2[i];
    z12 = zeta12[i];
  }
}

#define MAXNLAM 256
/*	Calculate zeta2 and zeta12 from zeta1        , lambda        , alpha ,
 * gamma and rho in the 2F model */
/*	Tau-TS enabled version */
void static_lgmcalczeta2zeta12_tauts(
    int n,          /*	Number of dates */
    double t[],     /*	Times */
    double zeta1[], /*	Zeta1 */
    /*	Lambda        , Alpha        , Beta        , Rho */
    int nlam, double lam_time[], double lam[], double alpha, double gamma,
    double rho,
    /*	Output */
    double zeta2[], double zeta12[]) {
  int i;
  double t0;
  double z1, z2, z12, dz1, dz2, dz12;
  double lam1[MAXNLAM], lam2[MAXNLAM], lam12[MAXNLAM];
  double exp_fact1, exp_fact2, exp_fact12;

  /*	Fill lambda tabs */
  for (i = 0; i < nlam; i++) {
    lam1[i] = lam[i];
    lam2[i] = lam[i] + gamma;
    lam12[i] = 0.5 * (lam1[i] + lam2[i]);
  }

  t0 = 0.0;
  z1 = z2 = z12 = 0.0;

  for (i = 0; i < n; i++) {
    exp_fact1 = static_lgmcalcexpfact_tauts(t0, t[i], nlam, lam_time, lam1);
    exp_fact2 =
        static_lgmcalcexpfact_tauts(t0, t[i], nlam, lam_time, lam2) / exp_fact1;
    exp_fact12 = static_lgmcalcexpfact_tauts(t0, t[i], nlam, lam_time, lam12) /
                 exp_fact1;

    dz1 = zeta1[i] - z1;
    dz2 = dz1 * alpha * alpha * exp_fact2;
    zeta2[i] = z2 + dz2;
    dz12 = dz1 * alpha * rho * exp_fact12;
    zeta12[i] = z12 + dz12;

    t0 = t[i];
    z1 = zeta1[i];
    z2 = zeta2[i];
    z12 = zeta12[i];
  }
}

/*	Main functions */

/*	Value of European option on a general set of cash flows within LGM 1F */
double lgmopval1F(int ncpn,       /*	Number of cash-flows */
                  double cpn[],   /*	Discounted Cash-Flows */
                  double cpn_G[], /*	G at cash-flow dates */
                  double ex_zeta, /*	Zeta at exercise date */
                  double ex_G)    /*	G at exercise date */
{
  double sqz = sqrt(ex_zeta);
  double ystar;
  int i;
  int dir;
  double val = 0.0;

  /*	Check for intrinsic */
  if (ex_zeta < 1.0e-10) {
    for (i = 0; i < ncpn; i++) {
      val += cpn[i];
    }

    return max(0.0, val);
  }

  ystar = static_lgmystar1F(ncpn, cpn, cpn_G, ex_zeta, sqz, ex_G, &dir);

  for (i = 0; i < ncpn; i++) {
    val += cpn[i] * static_lgmsafenorm(dir * (ystar + sqz * (cpn_G[i] - ex_G)));
  }

  return val;
}

/*	Value of European option on a general set of cash flows within LGM Stoch
 * Vol */
double lgmopval_stochvol(
    int ncpn,      /*	Number of cash-flows */
    double cpn[],  /*	Discounted Cash-Flows */
    double beta[], /*	beta(T*        , T) at cash-flow dates */
    double psi,    /*	Psi at exercise date */
    double driftf, /*  Std of f(t        ,T*) */
    double rho2,   /*	rho * rho */
    double sqrho2) /*	sqrt(1.0 - rho * rho) */

{
  double sqpsi = sqrt(psi);
  double stdf = sqpsi * sqrho2;
  double ystar;
  int i;
  int dir;
  double val = 0.0;

  /*	Check for intrinsic */
  if (psi < 1.0e-10) {
    for (i = 0; i < ncpn; i++) {
      val += cpn[i];
    }

    return val;
  }

  ystar = static_lgmystar_stochvol(ncpn, cpn, beta, psi, &dir);
  ystar -= driftf;
  ystar /= stdf;

  for (i = 0; i < ncpn; i++) {
    val += cpn[i] * exp(-beta[i] * (driftf + 0.5 * beta[i] * psi * rho2)) *
           static_lgmsafenorm(dir * (ystar + stdf * beta[i]));
  }

  return val;
}

/*	Value of European option on a general set of cash flows within LGM 2F */
double x[NUM_HERMITE + 1], w[NUM_HERMITE + 1];
double lgmopval2F(int ncpn,         /*	Number of cash-flows */
                  double cpn[],     /*	Discounted Cash-Flows */
                  double cpn_G1[],  /*	G1 at cash-flow dates */
                  double cpn_G2[],  /*	G2 at cash-flow dates */
                  double ex_zeta1,  /*	Z1 at exercise date */
                  double ex_zeta2,  /*	Z2 at exercise date */
                  double ex_zeta12, /*	Z12 at exercise date */
                  double ex_G1,     /*	G1 at exercise date */
                  double ex_G2)     /*	G2 at exercise date */
{
  double cstfact[MAX_CPN], n1fact[MAX_CPN], n2fact[MAX_CPN];
  double sqz1 = sqrt(ex_zeta1), sqz2 = sqrt(ex_zeta2);
  double c1 = ex_zeta12 / sqz1 / sqz2, c1sqz2 = c1 * sqz2,
         c2 = sqrt(1.0 - c1 * c1), c2sqz2 = c2 * sqz2;
  double *X = &(x[1]), *W = &(w[1]);
  double ystar;
  int i, j;
  int highest, dir;
  double temp, val = 0.0;

  /*	Check for intrinsic */
  if (ex_zeta1 < 1.0e-10 && ex_zeta2 < 1.0e-10) {
    for (i = 0; i < ncpn; i++) {
      val += cpn[i];
    }

    return max(0.0, val);
  }

  /*	Precalculate exponential factors in the IV */
  static_lgm2Fcalcexpfact(ncpn, cpn, cpn_G1, cpn_G2, ex_zeta1, sqz1, ex_zeta2,
                          sqz2, ex_zeta12, c1, c2, ex_G1, ex_G2, cstfact,
                          n1fact, n2fact, &highest);

  /*	Do Hermite integration */
  for (i = 0; i < NUM_HERMITE; i++) {
    /*	Disregard points of little weight */
    if (fabs(X[i]) > 5.0)
      continue;

    temp = 0.0;
    for (j = 0; j < ncpn; j++) {
      ystar = static_lgmystar2F(ncpn, cpn, cstfact, n1fact, n2fact, highest,
                                X[i] - (highest == 1 ? n2fact[j] : n1fact[j]),
                                &dir);

      temp +=
          cpn[j] * static_lgmsafenorm(
                       dir * (ystar + (highest == 1 ? n1fact[j] : n2fact[j])));
    }

    val += W[i] * temp;
  }

  return val;
}

/*	Value of European Swap option within LGM 1F */
double lgmsopval1F(int ncpn,     /*	Number of cash-flow dates        , including
                                                   start and end date */
                   double df[],  /*	Df to cash flow dates */
                   double cvg[], /*	cvg from i-1 to i */
                   double cpn_G[], /*	G at cash-flow dates */
                   double ex_zeta, /*	Z at exercise date */
                   double ex_G,    /*	G at exercise date */
                   double strike)  /*	Strike */
{
  int i;
  static double cpn[MAX_CPN];

  cpn[0] = -df[0];

  for (i = 1; i < ncpn; i++) {
    cpn[i] = df[i] * cvg[i] * strike;
  }

  cpn[ncpn - 1] += df[ncpn - 1];

  return lgmopval1F(ncpn, cpn, cpn_G, ex_zeta, ex_G);
}

/*	Value of European Swap option within LGM 2F */
double lgmsopval2F(int ncpn,     /*	Number of cash-flow dates        , including
                                                   start and end date */
                   double df[],  /*	Df to cash flow dates */
                   double cvg[], /*	cvg from i-1 to i */
                   double cpn_G1[],  /*	G1 at cash-flow dates */
                   double cpn_G2[],  /*	G2 at cash-flow dates */
                   double ex_zeta1,  /*	Z1 at exercise date */
                   double ex_zeta2,  /*	Z2 at exercise date */
                   double ex_zeta12, /*	Z12 at exercise date */
                   double ex_G1,     /*	G1 at exercise date */
                   double ex_G2,     /*	G2 at exercise date */
                   double strike)    /*	Strike */
{
  int i;
  static double cpn[MAX_CPN];

  cpn[0] = -df[0];

  for (i = 1; i < ncpn; i++) {
    cpn[i] = df[i] * cvg[i] * strike;
  }

  cpn[ncpn - 1] += df[ncpn - 1];

  return lgmopval2F(ncpn, cpn, cpn_G1, cpn_G2, ex_zeta1, ex_zeta2, ex_zeta12,
                    ex_G1, ex_G2);
}

/*	Value of European Cap within LGM 1F */
/*	NOTE: cap = sum of swaptions
                each swaption has for underlyings all the coupons from one
   exercise date to the following */
double lgmcapval1F(int ncpn,     /*	Number of cash-flow dates        , including
                                                   start and end date */
                   double df[],  /*	Df to cash flow dates */
                   double cvg[], /*	cvg from i-1 to i */
                   double cpn_G[], /*	G at cash-flow dates */
                   int nex,        /*	Number of exercise dates */
                   int ex_cpn[],   /*	For each exercise date        , first
                                      coupon   to be exercised */
                   int ex_ncpn[],  /*	For each exercise date        , number
                                of  coupons  to be exercised */
                   double ex_weight[],
                   double ex_zeta[], /*	Z at exercise date */
                   double ex_G[],    /*	G at exercise date */
                   double strike[])  /*	Strikes */
{
  double val;
  int i;

  val = 0.0;

  for (i = 0; i < nex; i++) {
    val += ex_weight[i] * lgmsopval1F(ex_ncpn[i], &(df[ex_cpn[i]]),
                                      &(cvg[ex_cpn[i]]), &(cpn_G[ex_cpn[i]]),
                                      ex_zeta[i], ex_G[i], strike[i]);
  }

  return val;
}

/*	Value of European Cap within LGM 2F */
/*	NOTE: cap = sum of swaptions
                each swaption has for underlyings all the coupons from one
   exercise date to the following */
double lgmcapval2F(int ncpn,     /*	Number of cash-flow dates        , including
                                                   start and end date */
                   double df[],  /*	Df to cash flow dates */
                   double cvg[], /*	cvg from i-1 to i */
                   double cpn_G1[],    /*	G1 at cash-flow dates */
                   double cpn_G2[],    /*	G2 at cash-flow dates */
                   int nex,            /*	Number of exercise dates */
                   int ex_cpn[],       /*	For each exercise date        , first
                                          coupon       to be exercised */
                   int ex_ncpn[],      /*	For each exercise date        , number
                                    of      coupons      to be exercised */
                   double ex_weight[], /*	Weights on each caplet */
                   double ex_zeta1[],  /*	Z1 at exercise date */
                   double ex_zeta2[],  /*	Z2 at exercise date */
                   double ex_zeta12[], /*	Z12 at exercise date */
                   double ex_G1[],     /*	G1 at exercise date */
                   double ex_G2[],     /*	G2 at exercise date */
                   double strike[])    /*	Strikes */
{
  double val;
  int i;

  val = 0.0;

  for (i = 0; i < nex; i++) {
    val += ex_weight[i] * lgmsopval2F(ex_ncpn[i], &(df[ex_cpn[i]]),
                                      &(cvg[ex_cpn[i]]), &(cpn_G1[ex_cpn[i]]),
                                      &(cpn_G2[ex_cpn[i]]), ex_zeta1[i],
                                      ex_zeta2[i], ex_zeta12[i], ex_G1[i],
                                      ex_G2[i], strike[i]);
  }

  return val;
}

/*	Setup G function from lambda */
static void static_lgmsetupG(double lambda,
                             int ncpn,          /*	Number of cash-flows */
                             double cpn_time[], /*	Cash-Flow times */
                             double cpn_G[],    /*	Output: G at cash-flow dates
                                                                          G(T)
                                             = (1.0 - exp (- lambda * T )) */
                             int nex,           /*	Number of exercise dates */
                             double ex_time[],  /*	Exercise times */
                             double ex_G[]) /*	Output: G at exercise dates */
{
  int i;

  for (i = 0; i < ncpn; i++) {
    cpn_G[i] = (1.0 - exp(-lambda * cpn_time[i])) / lambda;
  }

  for (i = 0; i < nex; i++) {
    ex_G[i] = (1.0 - exp(-lambda * ex_time[i])) / lambda;
  }
}

void export_lgmsetupG(double lambda, int ncpn, /*	Number of cash-flows */
                      double cpn_time[],       /*	Cash-Flow times */
                      double cpn_G[],          /*	Output: G at cash-flow dates
                                                                         G(T) =
                                            (1.0 - exp (- lambda * T )) */
                      int nex,                 /*	Number of exercise dates */
                      double ex_time[],        /*	Exercise times */
                      double ex_G[])           /*	Output: G at exercise dates */
{
  static_lgmsetupG(lambda, ncpn, cpn_time, cpn_G, nex, ex_time, ex_G);
}

void export_lgmsetupG_ts(int nlambda, double lambda_time[], double lambda[],

                         int ncpn,          /*	Number of cash-flows */
                         double cpn_time[], /*	Cash-Flow times */
                         double cpn_G[],    /*	Output: G at cash-flow dates
                                                                      G(T) =
                                         (1.0 - exp (- lambda * T )) */
                         int nex,           /*	Number of exercise dates */
                         double ex_time[],  /*	Exercise times */
                         double ex_G[])     /*	Output: G at exercise dates */
{
  static_lgmsetupG_tauts(nlambda, lambda_time, lambda, ncpn, cpn_time, cpn_G,
                         nex, ex_time, ex_G);
}

//// integrate lambda function from 0 to cpn_time[i]
double _average_lambda_(double dEndTime, int nlambda, double lambda_time[],
                        double lambda[]) {
  ///// first and last lambda times and values
  double dLambdaTime_t0, dLambdaTime_tlast;
  double dLambda_t0, dLambda_tlast;
  double dSumLambda, dSumTime;
  // double dPrevTime;
  int nI = 0, nJ;

  // if(dEndTime< 0.)
  //	{
  //		return "_average_lambda_(..) : end time must be positive!";
  //	}
#if 0
		if ( lambda[0] > 0.02+1.e-15 ||  lambda[0] < 0.02-1.e-15)
		{
			lambda[0] = 0.02;
		}
#endif
  dLambdaTime_t0 = lambda_time[0];
  dLambdaTime_tlast = lambda_time[nlambda - 1];

  dLambda_t0 = lambda[0];
  dLambda_tlast = lambda[nlambda - 1];

  //// 2 degenerate cases
  if (dEndTime <= dLambdaTime_t0)
    return dLambda_t0;

  if (dLambdaTime_tlast <= 0.)
    return dLambda_tlast;

  /// other cases
  dSumTime = lambda_time[0];
  dSumLambda = lambda_time[0] * lambda[0];
  nJ = 0;
  for (nI = 1; nI < nlambda && lambda_time[nI] < dEndTime; ++nI) {
    //// NB: term structure is left open right close - ( lo        , hi ]
    dSumTime += (lambda_time[nI] - lambda_time[nI - 1]);
    dSumLambda += (lambda_time[nI] - lambda_time[nI - 1]) * lambda[nI];
    nJ = nI;
  }

  nJ = ((nJ + 1) < nlambda) ? (nJ + 1) : (nlambda - 1);

  dSumLambda += (dEndTime - dSumTime) * lambda[nJ];
#if 0
		if( (dSumLambda/dEndTime)  > 0.02+1.e-15 ||   (dSumLambda/dEndTime)  < 0.02-1.e-15)
		{
			 (dSumLambda) = 0.02 ;
		}
#endif
  return dSumLambda / dEndTime;
}

/*	Calc G(T) function from lambda */
/*	Tau-TS enabled version */
static double static_lgmsetup_1G_tauts(int nlam, double lam_time[],
                                       double lam[], double T) {
  int i;
  double t0, ans, li;

  li = 0.0;
  i = 0;
  t0 = 0.0;
  ans = 0.0;
  while (i < nlam && lam_time[i] <= T) {
    ans += exp(-li) * (1.0 - exp(-lam[i] * (lam_time[i] - t0))) / lam[i];
    li += (lam_time[i] - t0) * lam[i];
    t0 = lam_time[i];
    i++;
  }
  if (i >= nlam)
    i--;
  ans += exp(-li) * (1.0 - exp(-lam[i] * (T - t0))) / lam[i];

  return ans;
}

/*	Setup G function from lambda */
/*	Tau-TS enabled version */
void static_lgmsetupG_tauts(int nlam, double lam_time[], double lam[],
                            int ncpn,          /*	Number of cash-flows */
                            double cpn_time[], /*	Cash-Flow times */
                            double cpn_G[],    /*	Output: G at cash-flow dates
                                                                         G(T)
                                            = (1.0 - exp (- lambda * T )) */
                            int nex,           /*	Number of exercise dates */
                            double ex_time[],  /*	Exercise times */
                            double ex_G[])     /*	Output: G at exercise dates */
{
  int i;

  for (i = 0; i < ncpn; i++) {
    cpn_G[i] = static_lgmsetup_1G_tauts(nlam, lam_time, lam, cpn_time[i]);
  }

  for (i = 0; i < nex; i++) {
    ex_G[i] = static_lgmsetup_1G_tauts(nlam, lam_time, lam, ex_time[i]);
  }
}

/*	Setup G2 function (2F) */
static void
static_lgmsetupG2(double lambda, double gamma,
                  int ncpn,          /*	Number of cash-flows */
                  double cpn_time[], /*	Cash-Flow times */
                  double cpn_G2[],   /*	Output: G2 at cash-flow dates
                                                               G2(T) = (1.0
                                  - exp (- lambda2 * T )) */
                  int nex,           /*	Number of exercise dates */
                  double ex_time[],  /*	Exercise times */
                  double ex_G2[])    /*	Output: G2 at exercise dates */
{
  int i;
  double lambda2 = lambda + gamma;

  for (i = 0; i < ncpn; i++) {
    cpn_G2[i] = (1.0 - exp(-lambda2 * cpn_time[i])) / lambda2;
  }

  for (i = 0; i < nex; i++) {
    ex_G2[i] = (1.0 - exp(-lambda2 * ex_time[i])) / lambda2;
  }
}

void export_lgmsetupG2(double lambda, double gamma,
                       int ncpn,          /*	Number of cash-flows */
                       double cpn_time[], /*	Cash-Flow times */
                       double cpn_G2[],   /*	Output: G2 at cash-flow dates
                                                                    G2(T) =
                                       (1.0 - exp (- lambda2 * T )) */
                       int nex,           /*	Number of exercise dates */
                       double ex_time[],  /*	Exercise times */
                       double ex_G2[])    /*	Output: G2 at exercise dates */
{
  static_lgmsetupG2(lambda, gamma, ncpn, cpn_time, cpn_G2, nex, ex_time, ex_G2);
}

void export_lgmsetupG2_ts(int nlambda, double lambda_time[], double lambda[],

                          double gamma, int ncpn, /*	Number of cash-flows */
                          double cpn_time[],      /*	Cash-Flow times */
                          double cpn_G2[],  /*	Output: G2 at cash-flow dates
                                                                      G2(T)
                                         = (1.0 - exp (- lambda2 * T )) */
                          int nex,          /*	Number of exercise dates */
                          double ex_time[], /*	Exercise times */
                          double ex_G2[])   /*	Output: G2 at exercise dates */
{
  // int			i;
  // double		dAvglambda2        , dAvglambda;// + gamma;

  static_lgmsetupG2_tauts(nlambda, lambda_time, lambda, gamma, ncpn, cpn_time,
                          cpn_G2, nex, ex_time, ex_G2);
}

/*	Setup G2 function (2F) */
/*	Tau-TS enabled version */
void static_lgmsetupG2_tauts(int nlam, double lam_time[], double lam[],
                             double gamma, int ncpn, /*	Number of cash-flows */
                             double cpn_time[],      /*	Cash-Flow times */
                             double cpn_G2[],        /*	Output: G2 at cash-flow
                                                  dates        G2(T) = (1.0 - exp (-
                                                  lambda2 * T )) */
                             int nex,          /*	Number of exercise dates */
                             double ex_time[], /*	Exercise times */
                             double ex_G2[]) /*	Output: G2 at exercise dates */
{
  int i;
  static double lam2[MAX_CPN];

  for (i = 0; i < nlam; i++) {
    lam2[i] = lam[i] + gamma;
  }

  for (i = 0; i < ncpn; i++) {
    cpn_G2[i] = static_lgmsetup_1G_tauts(nlam, lam_time, lam2, cpn_time[i]);
  }

  for (i = 0; i < nex; i++) {
    ex_G2[i] = static_lgmsetup_1G_tauts(nlam, lam_time, lam2, ex_time[i]);
  }
}

/*	Calibrate zeta to diagonal given G: 1F case */
Err lgmcalibzeta1F(
    int ncpn,           /*	Total number of cash-flow dates */
    double cpn_time[],  /*	Cash-Flow times */
    double cpn_df[],    /*	Df to cash-flow dates */
    double cpn_cvg[],   /*	cvg from i-1 to i */
    double cpn_G[],     /*	G at cash-flow dates */
    int nex,            /*	Total number of exercise dates */
    double ex_time[],   /*	Exercise times */
    int ex_cpn[],       /*	Index of the first cash-flow to be exercised */
    double ex_G[],      /*	G at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas */
    double lambda,
    int skip_last,   /*	If 1        , the last option is disregarded and the
                  forward   volatility is flat from option n-1 */
    double prec,     /*	Precision on primary instruments */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps)   /*	Allow vol term structure to jump */
{
  double s1, s2, s_last, ds;
  double t, zeta;
  int i, j, it;
  double dz1, dz2, z1, z2;
  double pr1, pr2;
  double exp_fact;
  double temp;
  double minvar, maxvar, d;
  double quad_var;
  int niter;

  /*	Initial guess: total vol of 85 - 115 bps */
  /*	s = local variance */
  s1 = 0.0085;
  s1 *= s1;
  s2 = 0.0115;
  s2 *= s2;

  /*	Initialisation */
  quad_var = 0.0;
  s_last = s1;
  ds = s2 / s1;
  zeta = 0.0;
  t = 0.0;

  /*	If only one option        , no skipping last */
  if (nex == 1) {
    skip_last = 0;
  }

  for (i = 0; i < nex - (skip_last > 0); i++) {
    j = ex_cpn[i];
    exp_fact =
        ((exp(2 * lambda * ex_time[i]) - exp(2 * lambda * t)) / 2 / lambda);

    /*	Initial guess: last vol */
    s1 = s_last;
    /*	Initial guess 2: last vol * ds */
    s2 = ds * s1;

    /*	First price */
    dz1 = s1 * exp_fact;
    z1 = zeta + dz1;

    pr1 = lgmsopval1F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G + j, z1, ex_G[i],
                      strike[i]);

    if (fabs(mkt_price[i] - pr1) < prec) {
      s_last = s1;
      quad_var += s1 * (ex_time[i] - t);
      ex_zeta[i] = zeta = z1;
      t = ex_time[i];
      continue;
    }

    /*	Second price */
    dz2 = s2 * exp_fact;
    z2 = zeta + dz2;

    pr2 = lgmsopval1F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G + j, z2, ex_G[i],
                      strike[i]);

    if (fabs(mkt_price[i] - pr2) < prec) {
      s_last = s2;
      quad_var += s2 * (ex_time[i] - t);
      ex_zeta[i] = zeta = z2;
      t = ex_time[i];
      continue;
    }

    /*	First option: no limits        , 15 iterations */
    if (i == 0) {
      niter = 3 * NITER;
      minvar = 1.0e-16;
      maxvar = 1.0;
    } else
    /*	Next ones: limited vol variation and only 5 iterations */
    {
      niter = NITER;
      minvar = quad_var / t * min_fact;
      maxvar = quad_var / t / max_fact;
    }

    d = 0.0;

    for (it = 1; it < niter; it++) {
      temp = s2;

      /*	Newton iteration */
      s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);

      /*	Out of lower bound */
      if (s2 < minvar) {
        /*	Calibrate to market price + 1bp */
        d = prec;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 < minvar) {
          s2 = minvar;
        }
      }

      /*	Out of upper bound */
      if (s2 > maxvar) {
        /*	Calibrate to market price - 1bp */
        d = -prec;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 > maxvar) {
          s2 = maxvar;
        }
      }

      if (fabs(s2 - temp) < 1.0e-16)
        break;

      s1 = temp;
      pr1 = pr2;

      /*	Reprice with new vol */
      dz2 = s2 * exp_fact;
      z2 = zeta + dz2;

      pr2 = lgmsopval1F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G + j, z2,
                        ex_G[i], strike[i]);

      if (fabs((mkt_price[i] + d) - pr2) < prec)
        break;
    }

    /*	If failed        , keep last vol */
    if (fabs((mkt_price[i] + d) - pr2) > prec) {
      smessage("Failed to calibrate for exercise date %d", i + 1);

      if (!use_jumps) {
        s2 = s_last;
      }

      dz2 = s2 * exp_fact;
      z2 = zeta + dz2;
    }

    s_last = s2;
    quad_var += s2 * (ex_time[i] - t);
    ex_zeta[i] = zeta = z2;
    t = ex_time[i];
  }

  /*	Replace last vol by last-1 if relevant */
  if (skip_last) {
    dz2 = s_last * ((exp(2 * lambda * ex_time[nex - 1]) - exp(2 * lambda * t)) /
                    2 / lambda);
    ex_zeta[nex - 1] = zeta + dz2;
  }

  return NULL;
}

/*	Calibrate zeta to diagonal given G: 1F case */
/*	Tau-TS enabled version */

Err lgmcalibzeta1F_tauts(
    int ncpn,           /*	Total number of cash-flow dates */
    double cpn_time[],  /*	Cash-Flow times */
    double cpn_df[],    /*	Df to cash-flow dates */
    double cpn_cvg[],   /*	cvg from i-1 to i */
    double cpn_G[],     /*	G at cash-flow dates */
    int nex,            /*	Total number of exercise dates */
    double ex_time[],   /*	Exercise times */
    int ex_cpn[],       /*	Index of the first cash-flow to be exercised */
    double ex_G[],      /*	G at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas */
    int nlam, double lam_time[], double lam[],
    int skip_last,   /*	If 1        , the last option is disregarded and the
                  forward   volatility is flat from option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps)   /*	Allow vol term structure to jump */
{
  double s1, s2, s_last, ds;
  double t, zeta;
  int i, j, it;
  double dz1, dz2, z1, z2;
  double pr1, pr2;
  double exp_fact;
  double temp;
  double minvar, maxvar, d;
  double quad_var;
  int niter;

  /*	Initial guess: total vol of 85 - 115 bps */
  /*	s = local variance */
  s1 = 0.0085;
  s1 *= s1;
  s2 = 0.0115;
  s2 *= s2;

  /*	Initialisation */
  quad_var = 0.0;
  s_last = s1;
  ds = s2 / s1;
  zeta = 0.0;
  t = 0.0;

  /*	If only one option        , no skipping last */
  if (nex == 1) {
    skip_last = 0;
  }

  for (i = 0; i < nex - (skip_last > 0); i++) {
    j = ex_cpn[i];
    exp_fact = static_lgmcalcexpfact_tauts(t, ex_time[i], nlam, lam_time, lam);

    /*	Initial guess: last vol */
    s1 = s_last;
    /*	Initial guess 2: last vol * ds */
    s2 = ds * s1;

    /*	First price */
    dz1 = s1 * exp_fact;
    z1 = zeta + dz1;

    pr1 = lgmsopval1F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G + j, z1, ex_G[i],
                      strike[i]);

    if (fabs(mkt_price[i] - pr1) < CALPRES) {
      s_last = s1;
      quad_var += s1 * (ex_time[i] - t);
      ex_zeta[i] = zeta = z1;
      t = ex_time[i];
      continue;
    }

    /*	Second price */
    dz2 = s2 * exp_fact;
    z2 = zeta + dz2;

    pr2 = lgmsopval1F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G + j, z2, ex_G[i],
                      strike[i]);

    if (fabs(mkt_price[i] - pr2) < CALPRES) {
      s_last = s2;
      quad_var += s2 * (ex_time[i] - t);
      ex_zeta[i] = zeta = z2;
      t = ex_time[i];
      continue;
    }

    /*	First option: no limits        , 15 iterations */
    if (i == 0) {
      niter = 3 * NITER;
      minvar = 1.0e-16;
      maxvar = 1.0;
    } else
    /*	Next ones: limited vol variation and only 5 iterations */
    {
      niter = NITER;
      minvar = quad_var / t * min_fact;
      maxvar = quad_var / t / max_fact;
    }

    d = 0.0;

    for (it = 1; it < niter; it++) {
      temp = s2;

      /*	Newton iteration */
      s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);

      /*	Out of lower bound */
      if (s2 < minvar) {
        /*	Calibrate to market price + 1bp */
        d = CALPRES;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 < minvar) {
          s2 = minvar;
        }
      }

      /*	Out of upper bound */
      if (s2 > maxvar) {
        /*	Calibrate to market price - 1bp */
        d = -CALPRES;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 > maxvar) {
          s2 = maxvar;
        }
      }

      if (fabs(s2 - temp) < 1.0e-16)
        break;

      s1 = temp;
      pr1 = pr2;

      /*	Reprice with new vol */
      dz2 = s2 * exp_fact;
      z2 = zeta + dz2;

      pr2 = lgmsopval1F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G + j, z2,
                        ex_G[i], strike[i]);

      if (fabs((mkt_price[i] + d) - pr2) < CALPRES)
        break;
    }

    /*	If failed        , keep last vol */
    if (fabs((mkt_price[i] + d) - pr2) > CALPRES) {
      smessage("Failed to calibrate for exercise date %d", i + 1);

      if (!use_jumps) {
        s2 = s_last;
      }

      dz2 = s2 * exp_fact;
      z2 = zeta + dz2;
    }

    s_last = s2;
    quad_var += s2 * (ex_time[i] - t);
    ex_zeta[i] = zeta = z2;
    t = ex_time[i];
  }

  /*	Replace last vol by last-1 if relevant */
  if (skip_last) {
    exp_fact =
        static_lgmcalcexpfact_tauts(t, ex_time[nex - 1], nlam, lam_time, lam);
    dz2 = s_last * exp_fact;
    ex_zeta[nex - 1] = zeta + dz2;
  }

  return NULL;
}

Err lgmcalibzeta1F_tauts2_old(
    int ncpn,          /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G[],    /*	G at cash-flow dates */
    int nex,           /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int ex_cpn[],      /*	Index of the first cash-flow to be exercised */
    /*	The extra argument */
    int ex_endcpn[], /*	Index of the last cash-flow to be exercised */
    /* */
    double ex_G[],      /*	G at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas */
    int nlam, double lam_time[], double lam[],
    int skip_last,   /*	If 1        , the last option is disregarded and the
                  forward   volatility is flat from option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact) /*	Maximum up jump on variance */
{
  double s1, s2, s_last, ds;
  double t, zeta;
  int i, j, it;
  double dz1, dz2, z1, z2;
  double pr1, pr2;
  double exp_fact;
  double temp;
  double minvar, maxvar, d;
  double quad_var;
  int niter;

  /*	Initial guess: total vol of 85 - 115 bps */
  /*	s = local variance */
  s1 = 0.0085;
  s1 *= s1;
  s2 = 0.0115;
  s2 *= s2;

  /*	Initialisation */
  quad_var = 0.0;
  s_last = s1;
  ds = s2 / s1;
  zeta = 0.0;
  t = 0.0;

  /*	If only one option        , no skipping last */
  if (nex == 1) {
    skip_last = 0;
  }

  for (i = 0; i < nex - (skip_last > 0); i++) {
    j = ex_cpn[i];
    exp_fact = static_lgmcalcexpfact_tauts(t, ex_time[i], nlam, lam_time, lam);

    /*	Initial guess: last vol */
    s1 = s_last;
    /*	Initial guess 2: last vol * ds */
    s2 = ds * s1;

    /*	First price */
    dz1 = s1 * exp_fact;
    z1 = zeta + dz1;

    pr1 = lgmsopval1F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G + j,
                      z1, ex_G[i], strike[i]);

    if (fabs(mkt_price[i] - pr1) < CALPRES) {
      s_last = s1;
      quad_var += s1 * (ex_time[i] - t);
      ex_zeta[i] = zeta = z1;
      t = ex_time[i];
      continue;
    }

    /*	Second price */
    dz2 = s2 * exp_fact;
    z2 = zeta + dz2;

    pr2 = lgmsopval1F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G + j,
                      z2, ex_G[i], strike[i]);

    if (fabs(mkt_price[i] - pr2) < CALPRES) {
      s_last = s2;
      quad_var += s2 * (ex_time[i] - t);
      ex_zeta[i] = zeta = z2;
      t = ex_time[i];
      continue;
    }

    /*	First option: no limits        , 15 iterations */
    if (i == 0) {
      niter = 3 * NITER;
      minvar = 1.0e-16;
      maxvar = 1.0;
    } else
    /*	Next ones: limited vol variation and only 5 iterations */
    {
      niter = NITER;
      minvar = quad_var / t * min_fact;
      maxvar = quad_var / t / max_fact;
    }

    d = 0.0;

    for (it = 1; it < niter; it++) {
      temp = s2;

      /*	Newton iteration */
      s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);

      /*	Out of lower bound */
      if (s2 < minvar) {
        /*	Calibrate to market price + 1bp */
        d = CALPRES;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 < minvar) {
          s2 = minvar;
        }
      }

      /*	Out of upper bound */
      if (s2 > maxvar) {
        /*	Calibrate to market price - 1bp */
        d = -CALPRES;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 > maxvar) {
          s2 = maxvar;
        }
      }

      if (fabs(s2 - temp) < 1.0e-16)
        break;

      s1 = temp;
      pr1 = pr2;

      /*	Reprice with new vol */
      dz2 = s2 * exp_fact;
      z2 = zeta + dz2;

      pr2 = lgmsopval1F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j,
                        cpn_G + j, z2, ex_G[i], strike[i]);

      if (fabs((mkt_price[i] + d) - pr2) < CALPRES)
        break;
    }

    /*	If failed        , keep last vol */
    if (fabs((mkt_price[i] + d) - pr2) > CALPRES) {
      smessage("Failed to calibrate for exercise date %d", i + 1);
      s2 = s_last;
      dz2 = s2 * exp_fact;
      z2 = zeta + dz2;
    }

    s_last = s2;
    quad_var += s2 * (ex_time[i] - t);
    ex_zeta[i] = zeta = z2;
    t = ex_time[i];
  }

  /*	Replace last vol by last-1 if relevant */
  if (skip_last) {
    exp_fact =
        static_lgmcalcexpfact_tauts(t, ex_time[nex - 1], nlam, lam_time, lam);
    dz2 = s_last * exp_fact;
    ex_zeta[nex - 1] = zeta + dz2;
  }

  return NULL;
}

/*	Calibrate zeta to diagonal given G: 1F case */
/*	Tau-TS enabled version */

Err lgmcalibzeta1F_tauts2_new(
    int ncpn,          /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G[],    /*	G at cash-flow dates */
    int nex,           /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int ex_cpn[],      /*	Index of the first cash-flow to be exercised */
    /*	The extra argument */
    int ex_endcpn[], /*	Index of the last cash-flow to be exercised */
    /* */
    double ex_G[],                       /*	G at exercise date */
    double strike[],                     /*	Strikes */
    double mkt_price[],                  /*	Market prices */
    double mkt_vega[], double ex_zeta[], /*	Output: zetas */
    int nlam, double lam_time[], double lam[],
    int skip_last, /*	If 1        , the last option is disregarded and the
                forward volatility is flat from option n-1 */
    double prec, int vega_prec,
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps)   /*	Use Jumps */
{
  double s1, s2, s_last, ds;
  double t, zeta;
  int i, j, it, l, k;
  double dz1, dz2, z1, z2;
  double pr1, pr2;
  double exp_fact;
  double temp;
  double minvar, maxvar, d;
  double quad_var;
  int niter, nIter = 10, iFirstMult = 5;
  double vega;
  /* the matrix for storing previous results */
  double **res_iter;

  /* allocate memory for the previous result matrix */
  res_iter = dmatrix(0, iFirstMult * nIter, 0, 1);

  /*	Initial guess: total vol of 85 - 115 bps */
  /*	s = local variance */
  s1 = 0.0085 * 0.0085;
  s2 = 0.0115 * 0.0115;

  /*	Initialisation */
  quad_var = 0.0;
  s_last = s1;
  ds = s2 / s1;
  zeta = 0.0;
  t = 0.0;

  /*	If only one option        , no skipping last */
  if (nex == 1) {
    skip_last = 0;
  }

  for (i = 0; i < nex - (skip_last > 0); i++) {
    /* calculate the vega */
    if (vega_prec && mkt_vega) {
      vega = mkt_vega[i];
    } else {
      vega = 1.0;
    }

    j = ex_cpn[i];
    exp_fact = static_lgmcalcexpfact_tauts(t, ex_time[i], nlam, lam_time, lam);

    /*	Initial guess: last vol */
    s1 = s_last;
    dz1 = s1 * exp_fact;
    z1 = zeta + dz1;
    pr1 = lgmsopval1F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G + j,
                      z1, ex_G[i], strike[i]);

    /*	Initial guess 2: last vol * ds */
    s2 = ds * s1;
    dz2 = s2 * exp_fact;
    z2 = zeta + dz2;
    pr2 = lgmsopval1F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G + j,
                      z2, ex_G[i], strike[i]);

    /* Check to see if guess 2 is okay */
    if (fabs(mkt_price[i] - pr2) / vega < prec) {
      s_last = s2;
      quad_var += s2 * (ex_time[i] - t);
      ex_zeta[i] = zeta = z2;
      t = ex_time[i];
      continue;
    }

    /*	First option: no limits on var        , 15 iterations */
    if (i == 0) {
      niter = iFirstMult * nIter;
      minvar = 1.0e-16;
      maxvar = 1.0;
    } else
    /*	Next ones: limited vol variation and only 5 iterations */
    {
      niter = nIter;
      minvar = quad_var / t * min_fact;
      maxvar = quad_var / t / max_fact;
    }

    d = 0.0;

    res_iter[0][0] = s1;
    res_iter[0][1] = pr1;
    res_iter[1][0] = s2;
    res_iter[1][1] = pr2;
    for (it = 1; it < niter; it++) {
      temp = s2;

      /*	Newton iteration */
      s2 = solve_for_next_coef(res_iter, it + 1, mkt_price[i], 0);

      /*	Out of lower bound */
      if (s2 < minvar) {
        /*	Calibrate to market price + 1bp
        d = prec;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                Still out of bounds
        if (s2 < minvar)
        {
                s2 = minvar;
        }*/
        s2 = minvar;
      }

      /*	Out of upper bound */
      if (s2 > maxvar) {
        /*	Calibrate to market price - 1bp
        d = -prec;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                Still out of bounds
        if (s2 > maxvar)
        {
                s2 = maxvar;
        } */
        s2 = maxvar;
      }

      /* Check to see if the var parameter has moved */
      if (fabs(s2 - temp) < 1.0e-16)
        break;

      s1 = temp;
      pr1 = pr2;

      /*	Reprice with new vol */
      dz2 = s2 * exp_fact;
      z2 = zeta + dz2;

      pr2 = lgmsopval1F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j,
                        cpn_G + j, z2, ex_G[i], strike[i]);

      /* if the price has converged then stop */
      if (fabs((mkt_price[i] + d) - pr2) / vega < prec)
        break;

      /* save the result */
      l = 0;
      while (l < it && res_iter[l][0] < s2)
        l++;

      if (l < it + 1)
        for (k = it; k >= l; k--) {
          res_iter[k + 1][0] = res_iter[k][0];
          res_iter[k + 1][1] = res_iter[k][1];
        }

      res_iter[l][0] = s2;
      res_iter[l][1] = pr2;
    }

    /*	If failed        , keep last vol */
    if (fabs((mkt_price[i] + d) - pr2) / vega > prec) {
      smessage("Failed to calibrate for exercise date %d", i + 1);

      if (!use_jumps) {
        s2 = s_last;
        dz2 = s2 * exp_fact;
        z2 = zeta + dz2;
      }
    }

    s_last = s2;
    quad_var += s2 * (ex_time[i] - t);
    ex_zeta[i] = zeta = z2;
    t = ex_time[i];
  }

  /*	Replace last vol by last-1 if relevant */
  if (skip_last) {
    exp_fact =
        static_lgmcalcexpfact_tauts(t, ex_time[nex - 1], nlam, lam_time, lam);
    dz2 = s_last * exp_fact;
    ex_zeta[nex - 1] = zeta + dz2;
  }

  /* free memory for the result matrix */
  free_dmatrix(res_iter, 0, iFirstMult * NITER, 0, 1);

  return NULL;
}

/*	Calibrate zeta to diagonal given G: 1F case */
/*	Tau-TS enabled version */
Err lgmcalibzeta1F_tauts2(
    int ncpn,          /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G[],    /*	G at cash-flow dates */
    int nex,           /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int ex_cpn[],      /*	Index of the first cash-flow to be exercised */
    /*	The extra argument */
    int ex_endcpn[], /*	Index of the last cash-flow to be exercised */
    /* */
    double ex_G[],                       /*	G at exercise date */
    double strike[],                     /*	Strikes */
    double mkt_price[],                  /*	Market prices */
    double mkt_vega[], double ex_zeta[], /*	Output: zetas */
    int nlam, double lam_time[], double lam[],
    int skip_last, /*	If 1        , the last option is disregarded and the
                forward volatility is flat from option n-1 */
    double prec, int vega_prec,
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps)   /*	Use Jumps */
{
  return lgmcalibzeta1F_tauts2_new(
      ncpn,     /*	Total number of cash-flow dates */
      cpn_time, /*	Cash-Flow times */
      cpn_df,   /*	Df to cash-flow dates */
      cpn_cvg,  /*	cvg from i-1 to i */
      cpn_G,    /*	G at cash-flow dates */
      nex,      /*	Total number of exercise dates */
      ex_time,  /*	Exercise times */
      ex_cpn,   /*	Index of the first cash-flow to be exercised */
      /*	The extra argument */
      ex_endcpn, /*	Index of the last cash-flow to be exercised */
      /* */
      ex_G,              /*	G at exercise date */
      strike,            /*	Strikes */
      mkt_price,         /*	Market prices */
      mkt_vega, ex_zeta, /*	Output: zetas */
      nlam, lam_time, lam,
      skip_last, /*	If 1        , the last option is disregarded and the
              forward volatility is flat from option n-1 */
      prec, vega_prec, min_fact, /*	Maximum down jump on variance */
      max_fact, use_jumps);      /*	Maximum up jump on variance */
}

/*	Calibrate zeta to swaptions given G: 1F case */
/*	New version: not necessarily the diagonal */
Err lgmcalibzeta1F_2(
    int ncpn,          /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G[],    /*	G at cash-flow dates */
    int nex,           /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int ex_cpn[],      /*	Index of the first cash-flow to be exercised */
    /*	The extra argument */
    int ex_endcpn[], /*	Index of the last cash-flow to be exercised */
    /* */
    double ex_G[],      /*	G at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas */
    double lambda,
    int skip_last) /*	If 1        , the last option is disregarded
                                     and the forward
                volatility is flat from option n-1 */
{
  double s1, s2, s_last, ds;
  double t, zeta;
  int i, j, it;
  double dz1, dz2, z1, z2;
  double pr1, pr2;
  double exp_fact;
  double temp;
  double minvar, maxvar, d;
  double quad_var;
  int niter;

  /*	Initial guess: total vol of 85 - 115 bps */
  /*	s = local variance */
  s1 = 0.0085;
  s1 *= s1;
  s2 = 0.0115;
  s2 *= s2;

  /*	Initialisation */
  quad_var = 0.0;
  s_last = s1;
  ds = s2 / s1;
  zeta = 0.0;
  t = 0.0;

  /*	If only one option        , no skipping last */
  if (nex == 1) {
    skip_last = 0;
  }

  for (i = 0; i < nex - (skip_last > 0); i++) {
    j = ex_cpn[i];
    exp_fact =
        ((exp(2 * lambda * ex_time[i]) - exp(2 * lambda * t)) / 2 / lambda);

    /*	Initial guess: last vol */
    s1 = s_last;
    /*	Initial guess 2: last vol * ds */
    s2 = ds * s1;

    /*	First price */
    dz1 = s1 * exp_fact;
    z1 = zeta + dz1;

    pr1 = lgmsopval1F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G + j,
                      z1, ex_G[i], strike[i]);

    if (fabs(mkt_price[i] - pr1) < CALPRES) {
      s_last = s1;
      quad_var += s1 * (ex_time[i] - t);
      ex_zeta[i] = zeta = z1;
      t = ex_time[i];
      continue;
    }

    /*	Second price */
    dz2 = s2 * exp_fact;
    z2 = zeta + dz2;

    pr2 = lgmsopval1F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G + j,
                      z2, ex_G[i], strike[i]);

    if (fabs(mkt_price[i] - pr2) < CALPRES) {
      s_last = s2;
      quad_var += s2 * (ex_time[i] - t);
      ex_zeta[i] = zeta = z2;
      t = ex_time[i];
      continue;
    }

    /*	First option: no limits        , 15 iterations */
    if (i == 0) {
      niter = 3 * NITER;
      minvar = 1.0e-16;
      maxvar = 1.0;
    } else
    /*	Next ones: limited vol variation and only 5 iterations */
    {
      niter = NITER;
      minvar = quad_var / t * MAX_FACT;
      maxvar = 4.0 * quad_var / t / MAX_FACT;
    }

    d = 0.0;

    for (it = 1; it < niter; it++) {
      temp = s2;

      /*	Newton iteration */
      s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);

      /*	Out of lower bound */
      if (s2 < minvar) {
        /*	Calibrate to market price + 1bp */
        d = CALPRES;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 < minvar) {
          s2 = minvar;
        }
      }

      /*	Out of upper bound */
      if (s2 > maxvar) {
        /*	Calibrate to market price - 1bp */
        d = -CALPRES;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 > maxvar) {
          s2 = maxvar;
        }
      }

      if (fabs(s2 - temp) < 1.0e-16)
        break;

      s1 = temp;
      pr1 = pr2;

      /*	Reprice with new vol */
      dz2 = s2 * exp_fact;
      z2 = zeta + dz2;

      pr2 = lgmsopval1F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j,
                        cpn_G + j, z2, ex_G[i], strike[i]);

      if (fabs((mkt_price[i] + d) - pr2) < CALPRES)
        break;
    }

    /*	If failed        , keep last vol */
    if (fabs((mkt_price[i] + d) - pr2) > CALPRES) {
      smessage("Failed to calibrate for exercise date %d", i + 1);
      s2 = s_last;
      dz2 = s2 * exp_fact;
      z2 = zeta + dz2;
    }

    s_last = s2;
    quad_var += s2 * (ex_time[i] - t);
    ex_zeta[i] = zeta = z2;
    t = ex_time[i];
  }

  /*	Replace last vol by last-1 if relevant */
  if (skip_last) {
    dz2 = s_last * ((exp(2 * lambda * ex_time[nex - 1]) - exp(2 * lambda * t)) /
                    2 / lambda);
    ex_zeta[nex - 1] = zeta + dz2;
  }

  return NULL;
}

/*	Calibrate zeta to diagonal given G: 2F case */
Err lgmcalibzeta2F(
    int ncpn,           /*	Total number of cash-flow dates */
    double cpn_time[],  /*	Cash-Flow times */
    double cpn_df[],    /*	Df to cash-flow dates */
    double cpn_cvg[],   /*	cvg from i-1 to i */
    double cpn_G1[],    /*	G1 at cash-flow dates */
    double cpn_G2[],    /*	G2 at cash-flow dates */
    int nex,            /*	Total number of exercise dates */
    double ex_time[],   /*	Exercise times */
    int ex_cpn[],       /*	Index of the first cash-flow to be exercised */
    double ex_G1[],     /*	G1 at exercise date */
    double ex_G2[],     /*	G2 at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas (1) */
    /*	Lambda        , Alpha        , gamma        , rho */
    double lambda, double alpha, double gamma, double rho,
    int skip_last,   /*	If 1        , the last option is disregarded and the
                  forward   volatility is flat from option n-1 */
    double prec,     /*	Precision on primary instruments */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps)   /*	Allow vol term structure to jump */
{
  double s1, s2, s_last, ds;
  int i, j, it;
  double l1 = lambda, l2 = lambda + gamma;
  double z1, z2, z12;
  double q1, q2, q12;
  double t, zeta1, zeta2, zeta12, dz1, dz2, dz12;
  double pr1, pr2;
  double exp_fact1, exp_fact2, exp_fact12;
  double temp;
  double minvar, maxvar, d;
  double quad_var;
  int niter;

  /*	Initial guess: total vol of 85 - 115 bps */
  /*	s = local variance of first factor */
  s1 = 0.0085 / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
  s1 *= s1;
  s2 = 0.0115 / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
  s2 *= s2;

  /*	Initialisation */
  quad_var = 0.0;
  s_last = s1;
  ds = s2 / s1;
  zeta1 = zeta2 = zeta12 = 0.0;
  t = 0.0;

  /*	If only one option        , no skipping last */
  if (nex == 1) {
    skip_last = 0;
  }

  for (i = 0; i < nex - (skip_last > 0); i++) {
    j = ex_cpn[i];

    if (fabs(l1) > 1.0E-10) {
      exp_fact1 = ((exp(2 * l1 * ex_time[i]) - exp(2 * l1 * t)) / 2 / l1);
    } else {
      exp_fact1 = ex_time[i] - t;
    }

    if (fabs(l2) > 1.0E-10) {
      exp_fact2 =
          ((exp(2 * l2 * ex_time[i]) - exp(2 * l2 * t)) / 2 / l2) / exp_fact1;
    } else {
      exp_fact2 = (ex_time[i] - t) / exp_fact1;
    }

    if (fabs(l1 + l2) > 1.0E-10) {
      exp_fact12 =
          ((exp((l1 + l2) * ex_time[i]) - exp((l1 + l2) * t)) / (l1 + l2)) /
          exp_fact1;
    } else {
      exp_fact12 = (ex_time[i] - t) / exp_fact1;
    }

    /*	Initial guess: last vol */
    s1 = s_last;
    /*	Initial guess 2: last vol * ds */
    s2 = ds * s1;

    /*	First price */
    dz1 = s1 * exp_fact1;
    z1 = zeta1 + dz1;

    dz2 = dz1 * alpha * alpha * exp_fact2;
    z2 = zeta2 + dz2;

    dz12 = dz1 * alpha * rho * exp_fact12;
    z12 = zeta12 + dz12;

    pr1 = lgmsopval2F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j, cpn_G2 + j,
                      z1, z2, z12, ex_G1[i], ex_G2[i], strike[i]);

    if (fabs(mkt_price[i] - pr1) < prec) {
      s_last = s1;
      quad_var += s1 * (ex_time[i] - t);
      ex_zeta[i] = zeta1 = z1;
      zeta2 = z2;
      zeta12 = z12;
      t = ex_time[i];
      continue;
    }

    /*	Second price */
    dz1 = s2 * exp_fact1;
    q1 = zeta1 + dz1;

    dz2 = dz1 * alpha * alpha * exp_fact2;
    q2 = zeta2 + dz2;

    dz12 = dz1 * alpha * rho * exp_fact12;
    q12 = zeta12 + dz12;

    pr2 = lgmsopval2F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j, cpn_G2 + j,
                      q1, q2, q12, ex_G1[i], ex_G2[i], strike[i]);

    if (fabs(mkt_price[i] - pr2) < prec) {
      s_last = s2;
      quad_var += s2 * (ex_time[i] - t);
      ex_zeta[i] = zeta1 = q1;
      zeta2 = q2;
      zeta12 = q12;
      t = ex_time[i];
      continue;
    }

    /*	First option: no limits        , 15 iterations */
    if (i == 0) {
      niter = 3 * NITER;
      minvar = 1.0e-16;
      maxvar = 1.0;
    } else
    /*	Next ones: limited vol variation and only 5 iterations */
    {
      niter = NITER;
      minvar = quad_var / t * min_fact;
      maxvar = quad_var / t / max_fact;
    }

    d = 0.0;

    for (it = 1; it < niter; it++) {
      temp = s2;

      /*	Newton iteration */
      s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);

      /*	Out of lower bound */
      if (s2 < minvar) {
        /*	Calibrate to market price + 1bp */
        d = prec;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 < minvar) {
          s2 = minvar;
        }
      }

      /*	Out of upper bound */
      if (s2 > maxvar) {
        /*	Calibrate to market price - 1bp */
        d = -prec;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 > maxvar) {
          s2 = maxvar;
        }
      }

      if (fabs(s2 - temp) < 1.0e-16)
        break;

      s1 = temp;
      pr1 = pr2;

      /*	Reprice with new vol */
      dz1 = s2 * exp_fact1;
      q1 = zeta1 + dz1;

      dz2 = dz1 * alpha * alpha * exp_fact2;
      q2 = zeta2 + dz2;

      dz12 = dz1 * alpha * rho * exp_fact12;
      q12 = zeta12 + dz12;

      pr2 = lgmsopval2F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j,
                        cpn_G2 + j, q1, q2, q12, ex_G1[i], ex_G2[i], strike[i]);

      if (fabs((mkt_price[i] + d) - pr2) < prec)
        break;
    }

    /*	If failed        , keep last vol */
    if (fabs((mkt_price[i] + d) - pr2) > prec) {
      smessage("Failed to calibrate for exercise date %d", i + 1);

      if (!use_jumps) {
        s2 = s_last;
      }

      dz1 = s2 * exp_fact1;
      dz2 = dz1 * alpha * alpha * exp_fact2;
      dz12 = dz1 * alpha * rho * exp_fact12;

      q1 = zeta1 + dz1;
      q2 = zeta2 + dz2;
      q12 = zeta12 + dz12;
    }

    s_last = s2;
    quad_var += s2 * (ex_time[i] - t);
    ex_zeta[i] = zeta1 = q1;
    zeta2 = q2;
    zeta12 = q12;
    t = ex_time[i];
  }

  /*	Replace last vol by last-1 if relevant */
  if (skip_last) {
    dz1 =
        s_last * ((exp(2 * l1 * ex_time[nex - 1]) - exp(2 * l1 * t)) / 2 / l1);
    ex_zeta[nex - 1] = zeta1 + dz1;
  }

  return NULL;
}

/*	Calibrate zeta to diagonal given G: 2F case */
/*	Tau-TS enabled version */

Err lgmcalibzeta2F_tauts(
    int ncpn,           /*	Total number of cash-flow dates */
    double cpn_time[],  /*	Cash-Flow times */
    double cpn_df[],    /*	Df to cash-flow dates */
    double cpn_cvg[],   /*	cvg from i-1 to i */
    double cpn_G1[],    /*	G1 at cash-flow dates */
    double cpn_G2[],    /*	G2 at cash-flow dates */
    int nex,            /*	Total number of exercise dates */
    double ex_time[],   /*	Exercise times */
    int ex_cpn[],       /*	Index of the first cash-flow to be exercised */
    double ex_G1[],     /*	G1 at exercise date */
    double ex_G2[],     /*	G2 at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas (1) */
    /*	Lambda        , Alpha        , gamma        , rho */
    int nlam, double lam_time[], double lam[], double alpha, double gamma,
    double rho,
    int skip_last,   /*	If 1        , the last option is disregarded and the
                  forward   volatility is flat from option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps)   /*	Allow vol term structure to jump */
{
  double s1, s2, s_last, ds;
  int i, j, it;
  double lam1[MAXNLAM], lam2[MAXNLAM], lam12[MAXNLAM];
  double z1, z2, z12;
  double q1, q2, q12;
  double t, zeta1, zeta2, zeta12, dz1, dz2, dz12;
  double pr1, pr2;
  double exp_fact1, exp_fact2, exp_fact12;
  double temp;
  double minvar, maxvar, d;
  double quad_var;
  int niter;

  /*	Fill lambda tabs */
  for (i = 0; i < nlam; i++) {
    lam1[i] = lam[i];
    lam2[i] = lam[i] + gamma;
    lam12[i] = 0.5 * (lam1[i] + lam2[i]);
  }

  /*	Initial guess: total vol of 85 - 115 bps */
  /*	s = local variance of first factor */
  s1 = 0.0085 / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
  s1 *= s1;
  s2 = 0.0115 / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
  s2 *= s2;

  /*	Initialisation */
  quad_var = 0.0;
  s_last = s1;
  ds = s2 / s1;
  zeta1 = zeta2 = zeta12 = 0.0;
  t = 0.0;

  /*	If only one option        , no skipping last */
  if (nex == 1) {
    skip_last = 0;
  }

  for (i = 0; i < nex - (skip_last > 0); i++) {
    j = ex_cpn[i];
    exp_fact1 =
        static_lgmcalcexpfact_tauts(t, ex_time[i], nlam, lam_time, lam1);
    exp_fact2 =
        static_lgmcalcexpfact_tauts(t, ex_time[i], nlam, lam_time, lam2) /
        exp_fact1;
    exp_fact12 =
        static_lgmcalcexpfact_tauts(t, ex_time[i], nlam, lam_time, lam12) /
        exp_fact1;

    /*	Initial guess: last vol */
    s1 = s_last;
    /*	Initial guess 2: last vol * ds */
    s2 = ds * s1;

    /*	First price */
    dz1 = s1 * exp_fact1;
    z1 = zeta1 + dz1;

    dz2 = dz1 * alpha * alpha * exp_fact2;
    z2 = zeta2 + dz2;

    dz12 = dz1 * alpha * rho * exp_fact12;
    z12 = zeta12 + dz12;

    pr1 = lgmsopval2F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j, cpn_G2 + j,
                      z1, z2, z12, ex_G1[i], ex_G2[i], strike[i]);

    if (fabs(mkt_price[i] - pr1) < CALPRES) {
      s_last = s1;
      quad_var += s1 * (ex_time[i] - t);
      ex_zeta[i] = zeta1 = z1;
      zeta2 = z2;
      zeta12 = z12;
      t = ex_time[i];
      continue;
    }

    /*	Second price */
    dz1 = s2 * exp_fact1;
    q1 = zeta1 + dz1;

    dz2 = dz1 * alpha * alpha * exp_fact2;
    q2 = zeta2 + dz2;

    dz12 = dz1 * alpha * rho * exp_fact12;
    q12 = zeta12 + dz12;

    pr2 = lgmsopval2F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j, cpn_G2 + j,
                      q1, q2, q12, ex_G1[i], ex_G2[i], strike[i]);

    if (fabs(mkt_price[i] - pr2) < CALPRES) {
      s_last = s2;
      quad_var += s2 * (ex_time[i] - t);
      ex_zeta[i] = zeta1 = q1;
      zeta2 = q2;
      zeta12 = q12;
      t = ex_time[i];
      continue;
    }

    /*	First option: no limits        , 15 iterations */
    if (i == 0) {
      niter = 3 * NITER;
      minvar = 1.0e-16;
      maxvar = 1.0;
    } else
    /*	Next ones: limited vol variation and only 5 iterations */
    {
      niter = NITER;
      minvar = quad_var / t * min_fact;
      maxvar = quad_var / t / max_fact;
    }

    d = 0.0;

    for (it = 1; it < niter; it++) {
      temp = s2;

      /*	Newton iteration */
      s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);

      /*	Out of lower bound */
      if (s2 < minvar) {
        /*	Calibrate to market price + 1bp */
        d = CALPRES;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 < minvar) {
          s2 = minvar;
        }
      }

      /*	Out of upper bound */
      if (s2 > maxvar) {
        /*	Calibrate to market price - 1bp */
        d = -CALPRES;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 > maxvar) {
          s2 = maxvar;
        }
      }

      if (fabs(s2 - temp) < 1.0e-16)
        break;

      s1 = temp;
      pr1 = pr2;

      /*	Reprice with new vol */
      dz1 = s2 * exp_fact1;
      q1 = zeta1 + dz1;

      dz2 = dz1 * alpha * alpha * exp_fact2;
      q2 = zeta2 + dz2;

      dz12 = dz1 * alpha * rho * exp_fact12;
      q12 = zeta12 + dz12;

      pr2 = lgmsopval2F(ncpn - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j,
                        cpn_G2 + j, q1, q2, q12, ex_G1[i], ex_G2[i], strike[i]);

      if (fabs((mkt_price[i] + d) - pr2) < CALPRES)
        break;
    }

    /*	If failed        , keep last vol */
    if (fabs((mkt_price[i] + d) - pr2) > CALPRES) {
      smessage("Failed to calibrate for exercise date %d", i + 1);

      if (!use_jumps) {
        s2 = s_last;
      }

      dz1 = s2 * exp_fact1;
      dz2 = dz1 * alpha * alpha * exp_fact2;
      dz12 = dz1 * alpha * rho * exp_fact12;

      q1 = zeta1 + dz1;
      q2 = zeta2 + dz2;
      q12 = zeta12 + dz12;
    }

    s_last = s2;
    quad_var += s2 * (ex_time[i] - t);
    ex_zeta[i] = zeta1 = q1;
    zeta2 = q2;
    zeta12 = q12;
    t = ex_time[i];
  }

  /*	Replace last vol by last-1 if relevant */
  if (skip_last) {
    exp_fact1 =
        static_lgmcalcexpfact_tauts(t, ex_time[nex - 1], nlam, lam_time, lam1);
    dz1 = s_last * exp_fact1;
    ex_zeta[nex - 1] = zeta1 + dz1;
  }

  return NULL;
}

/*	Calibrate zeta to diagonal given G: 2F case */
/*	Tau-TS enabled version */

Err lgmcalibzeta2F_tauts2(
    int ncpn,          /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G1[],   /*	G1 at cash-flow dates */
    double cpn_G2[],   /*	G2 at cash-flow dates */
    int nex,           /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int ex_cpn[],      /*	Index of the first cash-flow to be exercised */
    /*	The extra argument */
    int ex_endcpn[], /*	Index of the last cash-flow to be exercised */
    /* */
    double ex_G1[],     /*	G1 at exercise date */
    double ex_G2[],     /*	G2 at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double mkt_vega[],  /*	Market vega */
    double ex_zeta[],   /*	Output: zetas (1) */
    /*	Lambda        , Alpha        , gamma        , rho */
    int nlam, double lam_time[], double lam[], double alpha, double gamma,
    double rho,
    int skip_last, /*	If 1        , the last option is disregarded and the
                forward volatility is flat from option n-1 */
    double prec, int vega_prec,
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps)   /*	Use Jumps */
{
  double s1, s2, s_last, ds;
  int i, j, it;
  double lam1[MAXNLAM], lam2[MAXNLAM], lam12[MAXNLAM];
  double z1, z2, z12;
  double q1, q2, q12;
  double t, zeta1, zeta2, zeta12, dz1, dz2, dz12;
  double pr1, pr2;
  double exp_fact1, exp_fact2, exp_fact12;
  double temp;
  double minvar, maxvar, d;
  double quad_var;
  int niter;
  double vega;

  /*	Fill lambda tabs */
  for (i = 0; i < nlam; i++) {
    lam1[i] = lam[i];
    lam2[i] = lam[i] + gamma;
    lam12[i] = 0.5 * (lam1[i] + lam2[i]);
  }

  /*	Initial guess: total vol of 85 - 115 bps */
  /*	s = local variance of first factor */
  s1 = 0.0085 / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
  s1 *= s1;
  s2 = 0.0115 / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
  s2 *= s2;

  /*	Initialisation */
  quad_var = 0.0;
  s_last = s1;
  ds = s2 / s1;
  zeta1 = zeta2 = zeta12 = 0.0;
  t = 0.0;

  /*	If only one option        , no skipping last */
  if (nex == 1) {
    skip_last = 0;
  }

  for (i = 0; i < nex - (skip_last > 0); i++) {
    /* calculate the vega */
    if (vega_prec && mkt_vega) {
      vega = mkt_vega[i];
    } else {
      vega = 1.0;
    }

    j = ex_cpn[i];
    exp_fact1 =
        static_lgmcalcexpfact_tauts(t, ex_time[i], nlam, lam_time, lam1);
    exp_fact2 =
        static_lgmcalcexpfact_tauts(t, ex_time[i], nlam, lam_time, lam2) /
        exp_fact1;
    exp_fact12 =
        static_lgmcalcexpfact_tauts(t, ex_time[i], nlam, lam_time, lam12) /
        exp_fact1;

    /*	Initial guess: last vol */
    s1 = s_last;
    /*	Initial guess 2: last vol * ds */
    s2 = ds * s1;

    /*	First price */
    dz1 = s1 * exp_fact1;
    z1 = zeta1 + dz1;

    dz2 = dz1 * alpha * alpha * exp_fact2;
    z2 = zeta2 + dz2;

    dz12 = dz1 * alpha * rho * exp_fact12;
    z12 = zeta12 + dz12;

    pr1 = lgmsopval2F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j,
                      cpn_G2 + j, z1, z2, z12, ex_G1[i], ex_G2[i], strike[i]);

    if (fabs(mkt_price[i] - pr1) / vega < prec) {
      s_last = s1;
      quad_var += s1 * (ex_time[i] - t);
      ex_zeta[i] = zeta1 = z1;
      zeta2 = z2;
      zeta12 = z12;
      t = ex_time[i];
      continue;
    }

    /*	Second price */
    dz1 = s2 * exp_fact1;
    q1 = zeta1 + dz1;

    dz2 = dz1 * alpha * alpha * exp_fact2;
    q2 = zeta2 + dz2;

    dz12 = dz1 * alpha * rho * exp_fact12;
    q12 = zeta12 + dz12;

    pr2 = lgmsopval2F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j,
                      cpn_G2 + j, q1, q2, q12, ex_G1[i], ex_G2[i], strike[i]);

    if (fabs(mkt_price[i] - pr2) / vega < prec) {
      s_last = s2;
      quad_var += s2 * (ex_time[i] - t);
      ex_zeta[i] = zeta1 = q1;
      zeta2 = q2;
      zeta12 = q12;
      t = ex_time[i];
      continue;
    }

    /*	First option: no limits        , 15 iterations */
    if (i == 0) {
      niter = 3 * NITER;
      minvar = 1.0e-16;
      maxvar = 1.0;
    } else
    /*	Next ones: limited vol variation and only 5 iterations */
    {
      niter = NITER;
      minvar = quad_var / t * min_fact;
      maxvar = quad_var / t / max_fact;
    }

    d = 0.0;

    for (it = 1; it < niter; it++) {
      temp = s2;

      /*	Newton iteration */
      s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);

      /*	Out of lower bound */
      if (s2 < minvar) {
        /*	Calibrate to market price + 1bp */
        d = prec;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 < minvar) {
          s2 = minvar;
        }
      }

      /*	Out of upper bound */
      if (s2 > maxvar) {
        /*	Calibrate to market price - 1bp */
        d = -prec;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 > maxvar) {
          s2 = maxvar;
        }
      }

      if (fabs(s2 - temp) < 1.0e-16)
        break;

      s1 = temp;
      pr1 = pr2;

      /*	Reprice with new vol */
      dz1 = s2 * exp_fact1;
      q1 = zeta1 + dz1;

      dz2 = dz1 * alpha * alpha * exp_fact2;
      q2 = zeta2 + dz2;

      dz12 = dz1 * alpha * rho * exp_fact12;
      q12 = zeta12 + dz12;

      pr2 =
          lgmsopval2F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j,
                      cpn_G2 + j, q1, q2, q12, ex_G1[i], ex_G2[i], strike[i]);

      if (fabs((mkt_price[i] + d) - pr2) / vega < prec)
        break;
    }

    /*	If failed        , keep last vol */
    if (fabs((mkt_price[i] + d) - pr2) / vega > prec) {
      smessage("Failed to calibrate for exercise date %d", i + 1);

      if (!use_jumps) {
        s2 = s_last;

        dz1 = s2 * exp_fact1;
        dz2 = dz1 * alpha * alpha * exp_fact2;
        dz12 = dz1 * alpha * rho * exp_fact12;

        q1 = zeta1 + dz1;
        q2 = zeta2 + dz2;
        q12 = zeta12 + dz12;
      }
    }

    s_last = s2;
    quad_var += s2 * (ex_time[i] - t);
    ex_zeta[i] = zeta1 = q1;
    zeta2 = q2;
    zeta12 = q12;
    t = ex_time[i];
  }

  /*	Replace last vol by last-1 if relevant */
  if (skip_last) {
    exp_fact1 =
        static_lgmcalcexpfact_tauts(t, ex_time[nex - 1], nlam, lam_time, lam1);
    dz1 = s_last * exp_fact1;
    ex_zeta[nex - 1] = zeta1 + dz1;
  }

  return NULL;
}

/*	Calibrate zeta to diagonal given G: 2F case */
/*	New version: not necessarily the diagonal */
Err lgmcalibzeta2F_2(
    int ncpn,          /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G1[],   /*	G1 at cash-flow dates */
    double cpn_G2[],   /*	G2 at cash-flow dates */
    int nex,           /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int ex_cpn[],      /*	Index of the first cash-flow to be exercised */
    /*	The extra argument */
    int ex_endcpn[], /*	Index of the last cash-flow to be exercised */
    /* */
    double ex_G1[],     /*	G1 at exercise date */
    double ex_G2[],     /*	G2 at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas (1) */
    /*	Lambda        , Alpha        , gamma        , rho */
    double lambda, double alpha, double gamma, double rho,
    int skip_last) /*	If 1        , the last option is disregarded
                                           and the forward volatility is flat
                      from option n-1 */
{
  double s1, s2, s_last, ds;
  int i, j, it;
  double l1 = lambda, l2 = lambda + gamma;
  double z1, z2, z12;
  double q1, q2, q12;
  double t, zeta1, zeta2, zeta12, dz1, dz2, dz12;
  double pr1, pr2;
  double exp_fact1, exp_fact2, exp_fact12;
  double temp;
  double minvar, maxvar, d;
  double quad_var;
  int niter;

  /*	Initial guess: total vol of 85 - 115 bps */
  /*	s = local variance of first factor */
  s1 = 0.0085 / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
  s1 *= s1;
  s2 = 0.0115 / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
  s2 *= s2;

  /*	Initialisation */
  quad_var = 0.0;
  s_last = s1;
  ds = s2 / s1;
  zeta1 = zeta2 = zeta12 = 0.0;
  t = 0.0;

  /*	If only one option        , no skipping last */
  if (nex == 1) {
    skip_last = 0;
  }

  for (i = 0; i < nex - (skip_last > 0); i++) {
    j = ex_cpn[i];
    exp_fact1 = ((exp(2 * l1 * ex_time[i]) - exp(2 * l1 * t)) / 2 / l1);
    exp_fact2 = ((exp(2 * l2 * ex_time[i]) - exp(2 * l2 * t)) / 2 / l2) /
                ((exp(2 * l1 * ex_time[i]) - exp(2 * l1 * t)) / 2 / l1);
    exp_fact12 =
        ((exp((l1 + l2) * ex_time[i]) - exp((l1 + l2) * t)) / (l1 + l2)) /
        ((exp(2 * l1 * ex_time[i]) - exp(2 * l1 * t)) / 2 / l1);

    /*	Initial guess: last vol */
    s1 = s_last;
    /*	Initial guess 2: last vol * ds */
    s2 = ds * s1;

    /*	First price */
    dz1 = s1 * exp_fact1;
    z1 = zeta1 + dz1;

    dz2 = dz1 * alpha * alpha * exp_fact2;
    z2 = zeta2 + dz2;

    dz12 = dz1 * alpha * rho * exp_fact12;
    z12 = zeta12 + dz12;

    pr1 = lgmsopval2F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j,
                      cpn_G2 + j, z1, z2, z12, ex_G1[i], ex_G2[i], strike[i]);

    if (fabs(mkt_price[i] - pr1) < CALPRES) {
      s_last = s1;
      quad_var += s1 * (ex_time[i] - t);
      ex_zeta[i] = zeta1 = z1;
      zeta2 = z2;
      zeta12 = z12;
      t = ex_time[i];
      continue;
    }

    /*	Second price */
    dz1 = s2 * exp_fact1;
    q1 = zeta1 + dz1;

    dz2 = dz1 * alpha * alpha * exp_fact2;
    q2 = zeta2 + dz2;

    dz12 = dz1 * alpha * rho * exp_fact12;
    q12 = zeta12 + dz12;

    pr2 = lgmsopval2F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j,
                      cpn_G2 + j, q1, q2, q12, ex_G1[i], ex_G2[i], strike[i]);

    if (fabs(mkt_price[i] - pr2) < CALPRES) {
      s_last = s2;
      quad_var += s2 * (ex_time[i] - t);
      ex_zeta[i] = zeta1 = q1;
      zeta2 = q2;
      zeta12 = q12;
      t = ex_time[i];
      continue;
    }

    /*	First option: no limits        , 15 iterations */
    if (i == 0) {
      niter = 3 * NITER;
      minvar = 1.0e-16;
      maxvar = 1.0;
    } else
    /*	Next ones: limited vol variation and only 5 iterations */
    {
      niter = NITER;
      minvar = quad_var / t * MAX_FACT;
      maxvar = 4.0 * quad_var / t / MAX_FACT;
    }

    d = 0.0;

    for (it = 1; it < niter; it++) {
      temp = s2;

      /*	Newton iteration */
      s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);

      /*	Out of lower bound */
      if (s2 < minvar) {
        /*	Calibrate to market price + 1bp */
        d = CALPRES;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 < minvar) {
          s2 = minvar;
        }
      }

      /*	Out of upper bound */
      if (s2 > maxvar) {
        /*	Calibrate to market price - 1bp */
        d = -CALPRES;
        s2 = temp;
        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
        /*	Still out of bounds */
        if (s2 > maxvar) {
          s2 = maxvar;
        }
      }

      if (fabs(s2 - temp) < 1.0e-16)
        break;

      s1 = temp;
      pr1 = pr2;

      /*	Reprice with new vol */
      dz1 = s2 * exp_fact1;
      q1 = zeta1 + dz1;

      dz2 = dz1 * alpha * alpha * exp_fact2;
      q2 = zeta2 + dz2;

      dz12 = dz1 * alpha * rho * exp_fact12;
      q12 = zeta12 + dz12;

      pr2 =
          lgmsopval2F(ex_endcpn[i] + 1 - j, cpn_df + j, cpn_cvg + j, cpn_G1 + j,
                      cpn_G2 + j, q1, q2, q12, ex_G1[i], ex_G2[i], strike[i]);

      if (fabs((mkt_price[i] + d) - pr2) < CALPRES)
        break;
    }

    /*	If failed        , keep last vol */
    if (fabs((mkt_price[i] + d) - pr2) > CALPRES) {
      smessage("Failed to calibrate for exercise date %d", i + 1);
      s2 = s_last;

      dz1 = s2 * exp_fact1;
      dz2 = dz1 * alpha * alpha * exp_fact2;
      dz12 = dz1 * alpha * rho * exp_fact12;

      q1 = zeta1 + dz1;
      q2 = zeta2 + dz2;
      q12 = zeta12 + dz12;
    }

    s_last = s2;
    quad_var += s2 * (ex_time[i] - t);
    ex_zeta[i] = zeta1 = q1;
    zeta2 = q2;
    zeta12 = q12;
    t = ex_time[i];
  }

  /*	Replace last vol by last-1 if relevant */
  if (skip_last) {
    dz1 =
        s_last * ((exp(2 * l1 * ex_time[nex - 1]) - exp(2 * l1 * t)) / 2 / l1);
    ex_zeta[nex - 1] = zeta1 + dz1;
  }

  return NULL;
}

/*	Calibrate and price cap given lambda */
Err lgmprcapgivenlambda(
    int ncpn,            /*	Total number of cash-flow dates */
    double cpn_time[],   /*	Cash-Flow times */
    double cpn_df[],     /*	Df to cash-flow dates */
    double cpn_cvg[],    /*	cvg from i-1 to i */
    int nex,             /*	Total number of exercise dates */
    double ex_time[],    /*	Exercise times */
    int ex_cpn[],        /*	Index of the first cash-flow to be exercised */
    int ex_sncpn[],      /*	Number of coupons in each caplet */
    double ex_lstrike[], /*	Strikes for diagonal */
    double ex_lprice[],  /*	Market prices for diagonal */
    double ex_sstrike[], /*	Strikes for cap */
    double ex_sweight[], /*	Weights on each caplet */
    double ex_zeta[],    /*	Output: zetas */
    double lambda,       /*	Lambda */
    int one2F,           /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    double alpha, double gamma, double rho,
    int skip_last,     /*	If 1        , the last option is disregarded and the
                    forward     volatility is flat from option n-1 */
    double prec,       /*	Precision on primary instruments */
    double min_fact,   /*	Maximum down jump on variance */
    double max_fact,   /*	Maximum up jump on variance */
    int use_jumps,     /*	Allow vol term structure to jump */
    int price_cap,     /*	0: just calibrate */
    double *ex_sprice) /*	Cap price as output */
{
  static double cpn_G[MAX_CPN], cpn_G2[MAX_CPN], ex_G[MAX_CPN], ex_G2[MAX_CPN],
      zeta2[MAX_CPN], zeta12[MAX_CPN];
  Err err;

  /*	Setup G function */
  static_lgmsetupG(lambda, ncpn, cpn_time, cpn_G, nex, ex_time, ex_G);

  if (one2F == 2) {
    static_lgmsetupG2(lambda, gamma, ncpn, cpn_time, cpn_G2, nex, ex_time,
                      ex_G2);
  }

  if (one2F == 1) {
    err = lgmcalibzeta1F(ncpn, cpn_time, cpn_df, cpn_cvg, cpn_G, nex, ex_time,
                         ex_cpn, ex_G, ex_lstrike, ex_lprice, ex_zeta, lambda,
                         skip_last, prec, min_fact, max_fact, use_jumps);

    if (err) {
      return err;
    }

    if (price_cap) {
      *ex_sprice = lgmcapval1F(ncpn, cpn_df, cpn_cvg, cpn_G, nex, ex_cpn,
                               ex_sncpn, ex_sweight, ex_zeta, ex_G, ex_sstrike);
    }
  } else {
    err = lgmcalibzeta2F(ncpn, cpn_time, cpn_df, cpn_cvg, cpn_G, cpn_G2, nex,
                         ex_time, ex_cpn, ex_G, ex_G2, ex_lstrike, ex_lprice,
                         ex_zeta, lambda, alpha, gamma, rho, skip_last, prec,
                         min_fact, max_fact, use_jumps);

    if (err) {
      return err;
    }

    static_lgmcalczeta2zeta12(nex, ex_time, ex_zeta, lambda, alpha, gamma, rho,
                              zeta2, zeta12);

    if (price_cap) {
      *ex_sprice = lgmcapval2F(ncpn, cpn_df, cpn_cvg, cpn_G, cpn_G2, nex,
                               ex_cpn, ex_sncpn, ex_sweight, ex_zeta, zeta2,
                               zeta12, ex_G, ex_G2, ex_sstrike);
    }
  }

  return NULL;
}

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
/*	Tau-TS enabled version
 *** no lambda calibration so far for this version *** */
Err lgmprcapgivenlambda_tauts(
    int ncpn,            /*	Total number of cash-flow dates */
    double cpn_time[],   /*	Cash-Flow times */
    double cpn_df[],     /*	Df to cash-flow dates */
    double cpn_cvg[],    /*	cvg from i-1 to i */
    int nex,             /*	Total number of exercise dates */
    double ex_time[],    /*	Exercise times */
    int ex_cpn[],        /*	Index of the first cash-flow to be exercised */
    double ex_lstrike[], /*	Strikes */
    double ex_lprice[],  /*	Market prices */
    double ex_zeta[],    /*	Output: zetas */
    int one2F,           /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    int nlam, /*	Lambda TS: may NOT be changed in the process */
    double lam_time[], double lam[], double alpha, double gamma, double rho,
    int skip_last,   /*	If 1        , the last option is disregarded and the
                  forward   volatility is flat from option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps)   /*	Allow vol term structure to jump */

{
  static double cpn_G[MAX_CPN], cpn_G2[MAX_CPN], ex_G[MAX_CPN], ex_G2[MAX_CPN],
      zeta2[MAX_CPN], zeta12[MAX_CPN];
  Err err;

  /*	Setup G function */
  static_lgmsetupG_tauts(nlam, lam_time, lam, ncpn, cpn_time, cpn_G, nex,
                         ex_time, ex_G);

  if (one2F == 2) {
    static_lgmsetupG2_tauts(nlam, lam_time, lam, gamma, ncpn, cpn_time, cpn_G2,
                            nex, ex_time, ex_G2);
  }

  if (one2F == 1) {
    err = lgmcalibzeta1F_tauts(ncpn, cpn_time, cpn_df, cpn_cvg, cpn_G, nex,
                               ex_time, ex_cpn, ex_G, ex_lstrike, ex_lprice,
                               ex_zeta, nlam, lam_time, lam, skip_last,
                               min_fact, max_fact, use_jumps);

    if (err) {
      return err;
    }
  } else {
    err = lgmcalibzeta2F_tauts(
        ncpn, cpn_time, cpn_df, cpn_cvg, cpn_G, cpn_G2, nex, ex_time, ex_cpn,
        ex_G, ex_G2, ex_lstrike, ex_lprice, ex_zeta, nlam, lam_time, lam, alpha,
        gamma, rho, skip_last, min_fact, max_fact, use_jumps);

    if (err) {
      return err;
    }

    static_lgmcalczeta2zeta12_tauts(nex, ex_time, ex_zeta, nlam, lam_time, lam,
                                    alpha, gamma, rho, zeta2, zeta12);
  }

  return NULL;
}

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
/*	New version: calibrates not necessarily to digonal
 *** no lambda calibration so far for this version *** */
Err lgmprcapgivenlambda_2(
    int ncpn,            /*	Total number of cash-flow dates */
    double cpn_time[],   /*	Cash-Flow times */
    double cpn_df[],     /*	Df to cash-flow dates */
    double cpn_cvg[],    /*	cvg from i-1 to i */
    int nex,             /*	Total number of exercise dates */
    double ex_time[],    /*	Exercise times */
    int ex_cpn[],        /*	Index of the first cash-flow to be exercised */
    int ex_endcpn[],     /*	Index of the last cash-flow to be exercised */
    double ex_lstrike[], /*	Strikes */
    double ex_lprice[],  /*	Market prices */
    double ex_zeta[],    /*	Output: zetas */
    double lambda,       /*	Lambda: may NOT be changed in the process */
    int one2F,           /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    double alpha, double gamma, double rho,
    int skip_last) /*	If 1        , the last option is disregarded
                                           and the forward volatility is flat
                      from option n-1 */
{
  static double cpn_G[MAX_CPN], cpn_G2[MAX_CPN], ex_G[MAX_CPN], ex_G2[MAX_CPN],
      zeta2[MAX_CPN], zeta12[MAX_CPN];
  Err err;

  /*	Setup G function */
  static_lgmsetupG(lambda, ncpn, cpn_time, cpn_G, nex, ex_time, ex_G);

  if (one2F == 2) {
    static_lgmsetupG2(lambda, gamma, ncpn, cpn_time, cpn_G2, nex, ex_time,
                      ex_G2);
  }

  if (one2F == 1) {
    err = lgmcalibzeta1F_2(ncpn, cpn_time, cpn_df, cpn_cvg, cpn_G, nex, ex_time,
                           ex_cpn, ex_endcpn, ex_G, ex_lstrike, ex_lprice,
                           ex_zeta, lambda, skip_last);

    if (err) {
      return err;
    }
  } else {
    err = lgmcalibzeta2F_2(ncpn, cpn_time, cpn_df, cpn_cvg, cpn_G, cpn_G2, nex,
                           ex_time, ex_cpn, ex_endcpn, ex_G, ex_G2, ex_lstrike,
                           ex_lprice, ex_zeta, lambda, alpha, gamma, rho,
                           skip_last);

    if (err) {
      return err;
    }

    static_lgmcalczeta2zeta12(nex, ex_time, ex_zeta, lambda, alpha, gamma, rho,
                              zeta2, zeta12);
  }

  return NULL;
}

static int NCPN_LM, NEX_LM, *EX_CPN_LM, *EX_LENDCPN_LM, *EX_SENDCPN_LM,
    ONE2F_LM, SKIP_LAST_LM, NB_MOM_LM, FREQ_SHORT_LM, SHIFT_FREQ_LM, USE_JUMPS;

static double *CPN_TIME_LM, *CPN_DF_LM, *CPN_CVG_LM, *EX_TIME_LM,
    *EX_LSTRIKE_LM, *EX_LPRICE_LM, *EX_SSTRIKE_LM, *EX_SPRICE_LM, *EX_ZETA_LM,
    *LAM_TIME_LM, LAM_LM[MAX_CPN], ALPHA_LM, GAMMA_LM, RHO_LM,
    CPN_G_LM[MAX_CPN], CPN_G2_LM[MAX_CPN], EX_G_LM[MAX_CPN], EX_G2_LM[MAX_CPN],
    ZETA2_LM[MAX_CPN], ZETA12_LM[MAX_CPN], PRICE_LM[MAX_CPN],
    SENSI_LM[MAX_CPN][MAX_CPN];

void init_static_lgmprcapgiventauts(
    int ncpn,          /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    int nex,           /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int ex_cpn[],      /*	Index of the first cash-flow to be exercised */
    int ex_lendcpn[], double ex_lstrike[], /*	Strikes */
    double ex_lprice[],                    /*	Market prices */
    int ex_sendcpn[], double ex_sstrike[], /*	Strikes */
    double ex_sprice[],                    /*	Market prices */
    double ex_zeta[],                      /*	Output: zetas */
    int one2F,                             /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    int nlam, /*	Lambda TS: may NOT be changed in the process */
    double lam_time[], double lam[], double alpha, double gamma, double rho,
    int skip_last,
    DIAG_CALIB_LM_PARAMS ln_params) /*	If 1        , the last option is
                                       disregarded and the forward
                                       volatility is flat from option n-1 */
{
  NCPN_LM = ncpn;
  NEX_LM = nex;
  EX_CPN_LM = ex_cpn;
  EX_LENDCPN_LM = ex_lendcpn;
  EX_SENDCPN_LM = ex_sendcpn;
  ONE2F_LM = one2F;
  SKIP_LAST_LM = skip_last;

  CPN_TIME_LM = cpn_time;
  CPN_DF_LM = cpn_df;
  CPN_CVG_LM = cpn_cvg;
  EX_TIME_LM = ex_time;
  EX_LSTRIKE_LM = ex_lstrike;
  EX_LPRICE_LM = ex_lprice;
  EX_SSTRIKE_LM = ex_sstrike;
  EX_SPRICE_LM = ex_sprice;
  EX_ZETA_LM = ex_zeta;
  LAM_TIME_LM = lam_time;
  ALPHA_LM = alpha;
  GAMMA_LM = gamma;
  RHO_LM = rho;
  NB_MOM_LM = ln_params->nb_moment;
  FREQ_SHORT_LM = ln_params->freq_short;
  SHIFT_FREQ_LM = ln_params->shift_freq;
}

Err static_lgmprcapgiventauts(double index, double lam[], double *price,
                              double *gradient, int nlam) {
  static int i, j, ind;
  Err err;
  double *lambda = &(lam[1]);

  ind = (int)(index);

  if (ind - SHIFT_FREQ_LM == 0) {
    /* We calculate everything */
    if (ONE2F_LM == 1) {
      /* Calculate Sensi */
      for (j = 0; j < nlam; j++) {
        lambda[j] += BUMP_LAM_LM;
        static_lgmsetupG_tauts(nlam, LAM_TIME_LM, lambda, NCPN_LM, CPN_TIME_LM,
                               CPN_G_LM, NEX_LM, EX_TIME_LM, EX_G_LM);
        err = lgmcalibzeta1F_tauts2(
            NCPN_LM, CPN_TIME_LM, CPN_DF_LM, CPN_CVG_LM, CPN_G_LM, NEX_LM,
            EX_TIME_LM, EX_CPN_LM, EX_LENDCPN_LM, EX_G_LM, EX_LSTRIKE_LM,
            EX_LPRICE_LM, NULL, EX_ZETA_LM, nlam, LAM_TIME_LM, lambda,
            SKIP_LAST_LM, CALPRES, 0, MAX_FACT, MAX_FACT, 0);

        if (err) {
          return err;
        }

        for (i = 0; i < NEX_LM; i++) {
          if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM)) {
            SENSI_LM[i][j] = lgmsopval1F(
                EX_SENDCPN_LM[i] + 1 - EX_CPN_LM[i], &(CPN_DF_LM[EX_CPN_LM[i]]),
                &(CPN_CVG_LM[EX_CPN_LM[i]]), &(CPN_G_LM[EX_CPN_LM[i]]),
                EX_ZETA_LM[i], EX_G_LM[i], EX_SSTRIKE_LM[i]);
          }
        }

        lambda[j] -= BUMP_LAM_LM;
      }

      static_lgmsetupG_tauts(nlam, LAM_TIME_LM, lambda, NCPN_LM, CPN_TIME_LM,
                             CPN_G_LM, NEX_LM, EX_TIME_LM, EX_G_LM);

      err = lgmcalibzeta1F_tauts2(
          NCPN_LM, CPN_TIME_LM, CPN_DF_LM, CPN_CVG_LM, CPN_G_LM, NEX_LM,
          EX_TIME_LM, EX_CPN_LM, EX_LENDCPN_LM, EX_G_LM, EX_LSTRIKE_LM,
          EX_LPRICE_LM, NULL, EX_ZETA_LM, nlam, LAM_TIME_LM, lambda,
          SKIP_LAST_LM, CALPRES, 0, MAX_FACT, MAX_FACT, 0);

      if (err) {
        return err;
      }

      /* Calculate Prices */
      for (i = 0; i < NEX_LM; i++) {
        if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM)) {
          PRICE_LM[i] = lgmsopval1F(
              EX_SENDCPN_LM[i] + 1 - EX_CPN_LM[i], &(CPN_DF_LM[EX_CPN_LM[i]]),
              &(CPN_CVG_LM[EX_CPN_LM[i]]), &(CPN_G_LM[EX_CPN_LM[i]]),
              EX_ZETA_LM[i], EX_G_LM[i], EX_SSTRIKE_LM[i]);

          for (j = 0; j < nlam; j++) {
            SENSI_LM[i][j] = (SENSI_LM[i][j] - PRICE_LM[i]) / BUMP_LAM_LM;
          }
        }
      }
    } else {
      /* Calculate Sensi */
      for (j = 0; j < nlam; j++) {
        lambda[j] += BUMP_LAM_LM;
        static_lgmsetupG_tauts(nlam, LAM_TIME_LM, lambda, NCPN_LM, CPN_TIME_LM,
                               CPN_G_LM, NEX_LM, EX_TIME_LM, EX_G_LM);
        static_lgmsetupG2_tauts(nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPN_LM,
                                CPN_TIME_LM, CPN_G2_LM, NEX_LM, EX_TIME_LM,
                                EX_G2_LM);

        err = lgmcalibzeta2F_tauts2(
            NCPN_LM, CPN_TIME_LM, CPN_DF_LM, CPN_CVG_LM, CPN_G_LM, CPN_G2_LM,
            NEX_LM, EX_TIME_LM, EX_CPN_LM, EX_LENDCPN_LM, EX_G_LM, EX_G2_LM,
            EX_LSTRIKE_LM, EX_LPRICE_LM, NULL, EX_ZETA_LM, nlam, LAM_TIME_LM,
            lambda, ALPHA_LM, GAMMA_LM, RHO_LM, SKIP_LAST_LM, CALPRES, 0,
            MAX_FACT, MAX_FACT / 4.0, 0);

        if (err) {
          return err;
        }

        static_lgmcalczeta2zeta12_tauts(NEX_LM, EX_TIME_LM, EX_ZETA_LM, nlam,
                                        LAM_TIME_LM, lambda, ALPHA_LM, GAMMA_LM,
                                        RHO_LM, ZETA2_LM, ZETA12_LM);

        for (i = 0; i < NEX_LM; i++) {
          if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM)) {
            SENSI_LM[i][j] = lgmsopval2F(
                EX_SENDCPN_LM[i] + 1 - EX_CPN_LM[i], &(CPN_DF_LM[EX_CPN_LM[i]]),
                &(CPN_CVG_LM[EX_CPN_LM[i]]), &(CPN_G_LM[EX_CPN_LM[i]]),
                &(CPN_G2_LM[EX_CPN_LM[i]]), EX_ZETA_LM[i], ZETA2_LM[i],
                ZETA12_LM[i], EX_G_LM[i], EX_G2_LM[i], EX_SSTRIKE_LM[i]);
          }
        }

        lambda[j] -= BUMP_LAM_LM;
      }

      static_lgmsetupG_tauts(nlam, LAM_TIME_LM, lambda, NCPN_LM, CPN_TIME_LM,
                             CPN_G_LM, NEX_LM, EX_TIME_LM, EX_G_LM);
      static_lgmsetupG2_tauts(nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPN_LM,
                              CPN_TIME_LM, CPN_G2_LM, NEX_LM, EX_TIME_LM,
                              EX_G2_LM);

      err = lgmcalibzeta2F_tauts2(
          NCPN_LM, CPN_TIME_LM, CPN_DF_LM, CPN_CVG_LM, CPN_G_LM, CPN_G2_LM,
          NEX_LM, EX_TIME_LM, EX_CPN_LM, EX_LENDCPN_LM, EX_G_LM, EX_G2_LM,
          EX_LSTRIKE_LM, EX_LPRICE_LM, NULL, EX_ZETA_LM, nlam, LAM_TIME_LM,
          lambda, ALPHA_LM, GAMMA_LM, RHO_LM, SKIP_LAST_LM, CALPRES, 0,
          MAX_FACT, MAX_FACT / 4.0, 0);

      if (err) {
        return err;
      }

      static_lgmcalczeta2zeta12_tauts(NEX_LM, EX_TIME_LM, EX_ZETA_LM, nlam,
                                      LAM_TIME_LM, lambda, ALPHA_LM, GAMMA_LM,
                                      RHO_LM, ZETA2_LM, ZETA12_LM);

      for (i = 0; i < NEX_LM; i++) {
        if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM)) {
          PRICE_LM[i] = lgmsopval2F(
              EX_SENDCPN_LM[i] + 1 - EX_CPN_LM[i], &(CPN_DF_LM[EX_CPN_LM[i]]),
              &(CPN_CVG_LM[EX_CPN_LM[i]]), &(CPN_G_LM[EX_CPN_LM[i]]),
              &(CPN_G2_LM[EX_CPN_LM[i]]), EX_ZETA_LM[i], ZETA2_LM[i],
              ZETA12_LM[i], EX_G_LM[i], EX_G2_LM[i], EX_SSTRIKE_LM[i]);

          for (j = 0; j < nlam; j++) {
            SENSI_LM[i][j] = (SENSI_LM[i][j] - PRICE_LM[i]) / BUMP_LAM_LM;
          }
        }
      }
    }
  }

  /* Return the asked result */

  (*price) = PRICE_LM[ind];
  for (j = 0; j < nlam; j++) {
    gradient[j + 1] = SENSI_LM[ind][j];
  }

  return NULL;
}

/* This method try to match the momentum */
Err static_lgmprcapgiventauts_momentum(double index, double lam[],
                                       double *price, double *gradient,
                                       int nlam) {
  static int i, j, k, ind;
  Err err;
  double *lambda = &(lam[1]);
  static double temp, temp_pow[MAX_CPN];

  ind = (int)(index);

  if (ind == 0) {
    /* We calculate everything */
    if (ONE2F_LM == 1) {
      /* Calculate Sensi */
      for (j = 0; j < nlam; j++) {
        lambda[j] += BUMP_LAM_LM;
        static_lgmsetupG_tauts(nlam, LAM_TIME_LM, lambda, NCPN_LM, CPN_TIME_LM,
                               CPN_G_LM, NEX_LM, EX_TIME_LM, EX_G_LM);
        err = lgmcalibzeta1F_tauts2(
            NCPN_LM, CPN_TIME_LM, CPN_DF_LM, CPN_CVG_LM, CPN_G_LM, NEX_LM,
            EX_TIME_LM, EX_CPN_LM, EX_LENDCPN_LM, EX_G_LM, EX_LSTRIKE_LM,
            EX_LPRICE_LM, NULL, EX_ZETA_LM, nlam, LAM_TIME_LM, lambda,
            SKIP_LAST_LM, CALPRES, 0, MAX_FACT, MAX_FACT, 0);

        if (err) {
          return err;
        }

        for (i = 0; i < NEX_LM; i++) {
          if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM)) {
            SENSI_LM[i][j] =
                (lgmsopval1F(EX_SENDCPN_LM[i] + 1 - EX_CPN_LM[i],
                             &(CPN_DF_LM[EX_CPN_LM[i]]),
                             &(CPN_CVG_LM[EX_CPN_LM[i]]),
                             &(CPN_G_LM[EX_CPN_LM[i]]), EX_ZETA_LM[i],
                             EX_G_LM[i], EX_SSTRIKE_LM[i]) -
                 EX_SPRICE_LM[i]);
          }
        }

        lambda[j] -= BUMP_LAM_LM;
      }

      static_lgmsetupG_tauts(nlam, LAM_TIME_LM, lambda, NCPN_LM, CPN_TIME_LM,
                             CPN_G_LM, NEX_LM, EX_TIME_LM, EX_G_LM);

      err = lgmcalibzeta1F_tauts2(
          NCPN_LM, CPN_TIME_LM, CPN_DF_LM, CPN_CVG_LM, CPN_G_LM, NEX_LM,
          EX_TIME_LM, EX_CPN_LM, EX_LENDCPN_LM, EX_G_LM, EX_LSTRIKE_LM,
          EX_LPRICE_LM, NULL, EX_ZETA_LM, nlam, LAM_TIME_LM, lambda,
          SKIP_LAST_LM, CALPRES, 0, MAX_FACT, MAX_FACT, 0);

      if (err) {
        return err;
      }

      /* Calculate Prices */
      for (i = 0; i < NEX_LM; i++) {
        if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM)) {
          PRICE_LM[i] = (lgmsopval1F(EX_SENDCPN_LM[i] + 1 - EX_CPN_LM[i],
                                     &(CPN_DF_LM[EX_CPN_LM[i]]),
                                     &(CPN_CVG_LM[EX_CPN_LM[i]]),
                                     &(CPN_G_LM[EX_CPN_LM[i]]), EX_ZETA_LM[i],
                                     EX_G_LM[i], EX_SSTRIKE_LM[i]) -
                         EX_SPRICE_LM[i]);
        }
      }
    } else {
      /* Calculate Sensi */
      for (j = 0; j < nlam; j++) {
        lambda[j] += BUMP_LAM_LM;
        static_lgmsetupG_tauts(nlam, LAM_TIME_LM, lambda, NCPN_LM, CPN_TIME_LM,
                               CPN_G_LM, NEX_LM, EX_TIME_LM, EX_G_LM);
        static_lgmsetupG2_tauts(nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPN_LM,
                                CPN_TIME_LM, CPN_G2_LM, NEX_LM, EX_TIME_LM,
                                EX_G2_LM);

        err = lgmcalibzeta2F_tauts2(
            NCPN_LM, CPN_TIME_LM, CPN_DF_LM, CPN_CVG_LM, CPN_G_LM, CPN_G2_LM,
            NEX_LM, EX_TIME_LM, EX_CPN_LM, EX_LENDCPN_LM, EX_G_LM, EX_G2_LM,
            EX_LSTRIKE_LM, EX_LPRICE_LM, NULL, EX_ZETA_LM, nlam, LAM_TIME_LM,
            lambda, ALPHA_LM, GAMMA_LM, RHO_LM, SKIP_LAST_LM, CALPRES, 0,
            MAX_FACT, MAX_FACT / 4.0, 0);

        if (err) {
          return err;
        }

        static_lgmcalczeta2zeta12_tauts(NEX_LM, EX_TIME_LM, EX_ZETA_LM, nlam,
                                        LAM_TIME_LM, lambda, ALPHA_LM, GAMMA_LM,
                                        RHO_LM, ZETA2_LM, ZETA12_LM);

        for (i = 0; i < NEX_LM; i++) {
          if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM)) {
            SENSI_LM[i][j] =
                (lgmsopval2F(
                     EX_SENDCPN_LM[i] + 1 - EX_CPN_LM[i],
                     &(CPN_DF_LM[EX_CPN_LM[i]]), &(CPN_CVG_LM[EX_CPN_LM[i]]),
                     &(CPN_G_LM[EX_CPN_LM[i]]), &(CPN_G2_LM[EX_CPN_LM[i]]),
                     EX_ZETA_LM[i], ZETA2_LM[i], ZETA12_LM[i], EX_G_LM[i],
                     EX_G2_LM[i], EX_SSTRIKE_LM[i]) -
                 EX_SPRICE_LM[i]);
          }
        }

        lambda[j] -= BUMP_LAM_LM;
      }

      static_lgmsetupG_tauts(nlam, LAM_TIME_LM, lambda, NCPN_LM, CPN_TIME_LM,
                             CPN_G_LM, NEX_LM, EX_TIME_LM, EX_G_LM);
      static_lgmsetupG2_tauts(nlam, LAM_TIME_LM, lambda, GAMMA_LM, NCPN_LM,
                              CPN_TIME_LM, CPN_G2_LM, NEX_LM, EX_TIME_LM,
                              EX_G2_LM);

      err = lgmcalibzeta2F_tauts2(
          NCPN_LM, CPN_TIME_LM, CPN_DF_LM, CPN_CVG_LM, CPN_G_LM, CPN_G2_LM,
          NEX_LM, EX_TIME_LM, EX_CPN_LM, EX_LENDCPN_LM, EX_G_LM, EX_G2_LM,
          EX_LSTRIKE_LM, EX_LPRICE_LM, NULL, EX_ZETA_LM, nlam, LAM_TIME_LM,
          lambda, ALPHA_LM, GAMMA_LM, RHO_LM, SKIP_LAST_LM, CALPRES, 0,
          MAX_FACT, MAX_FACT / 4.0, 0);

      if (err) {
        return err;
      }

      static_lgmcalczeta2zeta12_tauts(NEX_LM, EX_TIME_LM, EX_ZETA_LM, nlam,
                                      LAM_TIME_LM, lambda, ALPHA_LM, GAMMA_LM,
                                      RHO_LM, ZETA2_LM, ZETA12_LM);

      for (i = 0; i < NEX_LM; i++) {
        if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM)) {
          PRICE_LM[i] = (lgmsopval2F(EX_SENDCPN_LM[i] + 1 - EX_CPN_LM[i],
                                     &(CPN_DF_LM[EX_CPN_LM[i]]),
                                     &(CPN_CVG_LM[EX_CPN_LM[i]]),
                                     &(CPN_G_LM[EX_CPN_LM[i]]),
                                     &(CPN_G2_LM[EX_CPN_LM[i]]), EX_ZETA_LM[i],
                                     ZETA2_LM[i], ZETA12_LM[i], EX_G_LM[i],
                                     EX_G2_LM[i], EX_SSTRIKE_LM[i]) -
                         EX_SPRICE_LM[i]);
        }
      }
    }

    /* Computes the momentum */
    for (k = 0; k < NB_MOM_LM; k++) {
      temp = 0.0;

      for (i = 0; i < NEX_LM; i++) {
        if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM)) {
          temp += pow(PRICE_LM[i], k + 1);
        }
      }

      PRICE_LM[NEX_LM + k] = temp / fabs(temp) * exp(log(fabs(temp)) / (k + 1));
    }

    for (j = 0; j < nlam; j++) {
      for (k = 0; k < NB_MOM_LM; k++) {
        temp = 0.0;

        for (i = 0; i < NEX_LM; i++) {
          if (!fmod(i + SHIFT_FREQ_LM, FREQ_SHORT_LM)) {
            temp += pow(SENSI_LM[i][j], k + 1);
          }
        }

        SENSI_LM[NEX_LM + k][j] =
            (temp / fabs(temp) * exp(log(fabs(temp)) / (k + 1)) -
             PRICE_LM[NEX_LM + k]) /
            BUMP_LAM_LM;
      }
    }
  }

  /* Return the asked result */

  (*price) = PRICE_LM[NEX_LM + ind];
  for (j = 0; j < nlam; j++) {
    gradient[j + 1] = SENSI_LM[NEX_LM + ind][j];
  }

  return NULL;
}

Err lgmprcapgivenlambda_tauts2(
    int ncpn,          /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    int nex,           /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int ex_cpn[],      /*	Index of the first cash-flow to be exercised */
    int ex_lendcpn[], double ex_lstrike[], /*	Strikes */
    double ex_lprice[],                    /*	Market prices */
    double ex_lvega[], int ex_sendcpn[],
    double ex_sstrike[], /*	Strikes */
    double ex_zeta[],    /*	Output: zetas */
    int one2F,           /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    int nlam, /*	Lambda TS: may NOT be changed in the process */
    double lam_time[], double lam[], double alpha, double gamma, double rho,
    int skip_last, /*	If 1        , the last option is disregarded and the
                forward volatility is flat from option n-1 */
    double prec, int vega_prec,
    double min_fact,       /*	Maximum down jump on variance */
    double max_fact,       /*	Maximum up jump on variance */
    int use_jumps,         /*	Use Jumps */
    int price_cap,         /*	0: just calibrate */
    double ex_spriceres[]) /*	Cap price as output */
{
  static double cpn_G[MAX_CPN], cpn_G2[MAX_CPN], ex_G[MAX_CPN], ex_G2[MAX_CPN],
      zeta2[MAX_CPN], zeta12[MAX_CPN];
  int i;
  Err err;

  /*	Setup G function */
  static_lgmsetupG_tauts(nlam, lam_time, lam, ncpn, cpn_time, cpn_G, nex,
                         ex_time, ex_G);

  if (one2F == 2) {
    static_lgmsetupG2_tauts(nlam, lam_time, lam, gamma, ncpn, cpn_time, cpn_G2,
                            nex, ex_time, ex_G2);
  }

  if (one2F == 1) {
    err = lgmcalibzeta1F_tauts2(ncpn, cpn_time, cpn_df, cpn_cvg, cpn_G, nex,
                                ex_time, ex_cpn, ex_lendcpn, ex_G, ex_lstrike,
                                ex_lprice, ex_lvega, ex_zeta, nlam, lam_time,
                                lam, skip_last, prec, vega_prec, min_fact,
                                max_fact, use_jumps);

    if (err) {
      return err;
    }
  } else {
    err = lgmcalibzeta2F_tauts2(ncpn, cpn_time, cpn_df, cpn_cvg, cpn_G, cpn_G2,
                                nex, ex_time, ex_cpn, ex_lendcpn, ex_G, ex_G2,
                                ex_lstrike, ex_lprice, ex_lvega, ex_zeta, nlam,
                                lam_time, lam, alpha, gamma, rho, skip_last,
                                prec, vega_prec, min_fact, max_fact, use_jumps);

    if (err) {
      return err;
    }

    static_lgmcalczeta2zeta12_tauts(nex, ex_time, ex_zeta, nlam, lam_time, lam,
                                    alpha, gamma, rho, zeta2, zeta12);
  }

  if (price_cap) {
    if (one2F == 1) {
      for (i = 0; i < nex; i++) {
        ex_spriceres[i] =
            lgmsopval1F(ex_sendcpn[i] + 1 - ex_cpn[i], &(cpn_df[ex_cpn[i]]),
                        &(cpn_cvg[ex_cpn[i]]), &(cpn_G[ex_cpn[i]]), ex_zeta[i],
                        ex_G[i], ex_sstrike[i]);
      }
    } else {
      for (i = 0; i < nex; i++) {
        ex_spriceres[i] = lgmsopval2F(
            ex_sendcpn[i] + 1 - ex_cpn[i], &(cpn_df[ex_cpn[i]]),
            &(cpn_cvg[ex_cpn[i]]), &(cpn_G[ex_cpn[i]]), &(cpn_G2[ex_cpn[i]]),
            ex_zeta[i], zeta2[i], zeta12[i], ex_G[i], ex_G2[i], ex_sstrike[i]);
      }
    }
  }

  return NULL;
}

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
Err lgmcalibzetalambda(
    int ncpn,            /*	Total number of cash-flow dates */
    double cpn_time[],   /*	Cash-Flow times */
    double cpn_df[],     /*	Df to cash-flow dates */
    double cpn_cvg[],    /*	cvg from i-1 to i */
    int nex,             /*	Total number of exercise dates */
    double ex_time[],    /*	Exercise times */
    int ex_cpn[],        /*	Index of the first cash-flow to be exercised */
    int ex_sncpn[],      /*	Number of coupons in each caplet */
    double ex_lstrike[], /*	Strikes for diagonal */
    double ex_lprice[],  /*	Market prices for diagonal */
    double ex_sstrike[], /*	Strikes for cap */
    double ex_sweight[], /*	Weights on each caplet */
    double ex_sprice,    /*	Market price for cap */
    double ex_zeta[],    /*	Output: zetas */
    int fix_lambda,      /*	0: calib lambda to cap        , 1: fix lambda calib
                                           to diagonal */
    double *lambda,      /*	Lambda: may be changed in the process */
    int one2F,           /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    double alpha, double gamma, double rho,
    int skip_last,     /*	If 1        , the last option is disregarded and the
                    forward     volatility is flat from option n-1 */
    double long_prec,  /*	Precision on primary instruments */
    double short_prec, /*	Precision on secondary instruments */
    double min_fact,   /*	Maximum down jump on variance */
    double max_fact,   /*	Maximum up jump on variance */
    int use_jumps)     /*	Allow vol term structure to jump */
{
  int it;
  double lam1, lam2;
  double pr1, pr2;
  double fact, temp, temppr;
  int ifact;
  Err err;

  if (!fix_lambda && nex < 2) {
    return serror("Cannot calibrate lambda and zeta with less than 2 exercises "
                  "- choose fix lambda");
  }

  err = lgmprcapgivenlambda(ncpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time,
                            ex_cpn, ex_sncpn, ex_lstrike, ex_lprice, ex_sstrike,
                            ex_sweight, ex_zeta, *lambda, one2F, alpha, gamma,
                            rho, skip_last, long_prec, min_fact, max_fact,
                            use_jumps, (fix_lambda ? 0 : 1), &pr1);

  if (err || fix_lambda) {
    return err;
  }

  if (fabs(ex_sprice - pr1) < short_prec) {
    return NULL;
  }

  if (pr1 > ex_sprice) {
    ifact = -1;
  } else {
    ifact = 1;
  }

  lam1 = *lambda;
  fact = 0.25 * fabs(lam1);
  if (fact < 0.01) {
    fact = 0.01;
  }

  do {
    pr2 = pr1;

    fact *= 1.5;
    lam1 += ifact * fact;

    if (err = lgmprcapgivenlambda(ncpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time,
                                  ex_cpn, ex_sncpn, ex_lstrike, ex_lprice,
                                  ex_sstrike, ex_sweight, ex_zeta, lam1, one2F,
                                  alpha, gamma, rho, skip_last, long_prec,
                                  min_fact, max_fact, use_jumps, 1, &pr1)) {
      return err;
    }

    if (fabs(ex_sprice - pr1) < short_prec) {
      *lambda = lam1;
      return NULL;
    }
  } while (ifact * pr1 < ifact * ex_sprice);

  lam2 = lam1 - ifact * fact;

  if (lam1 > lam2) {
    temp = lam1;
    temppr = pr1;
    lam1 = lam2;
    pr1 = pr2;
    lam2 = temp;
    pr2 = temppr;
  }

  for (it = 1; it < 20; it++) {
    temp = lam2;
    temppr = pr2;
    //		lam2 = 0.5 * (lam1 + lam2);
    lam2 = lam1 + (ex_sprice - pr1) * (lam2 - lam1) / (pr2 - pr1);

    if (err = lgmprcapgivenlambda(ncpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time,
                                  ex_cpn, ex_sncpn, ex_lstrike, ex_lprice,
                                  ex_sstrike, ex_sweight, ex_zeta, lam2, one2F,
                                  alpha, gamma, rho, skip_last, long_prec,
                                  min_fact, max_fact, use_jumps, 1, &pr2)) {
      return err;
    }

    if (fabs(ex_sprice - pr2) < short_prec) {
      *lambda = lam2;
      return NULL;
    }

    if (pr2 < ex_sprice) {
      pr1 = pr2;
      lam1 = lam2;
      lam2 = temp;
      pr2 = temppr;
    }

    if (fabs(lam2 - lam1) < short_prec || fabs(pr2 - pr1) < short_prec)
      break;

    while (pr1 > pr2 && it < 10) {
      it++;

      temp = lam1;
      lam1 = lam2;
      pr1 = pr2;
      lam2 += lam2 - temp;

      if (err = lgmprcapgivenlambda(
              ncpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time, ex_cpn, ex_sncpn,
              ex_lstrike, ex_lprice, ex_sstrike, ex_sweight, ex_zeta, lam2,
              one2F, alpha, gamma, rho, skip_last, long_prec, min_fact,
              max_fact, use_jumps, 1, &pr2)) {
        return err;
      }

      if (fabs(ex_sprice - pr2) < short_prec) {
        *lambda = lam2;
        return NULL;
      }
    }
  }

  smessage("Calibration of lambda        , error in bp: %.00f",
           10000 * min(fabs(ex_sprice - pr1), fabs(ex_sprice - pr2)));

  *lambda = lam2;

  return NULL;
}

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
Err lgmcalibzetalambda_tauts(
    int ncpn,            /*	Total number of cash-flow dates */
    double cpn_time[],   /*	Cash-Flow times */
    double cpn_df[],     /*	Df to cash-flow dates */
    double cpn_cvg[],    /*	cvg from i-1 to i */
    int nex,             /*	Total number of exercise dates */
    double ex_time[],    /*	Exercise times */
    int ex_cpn[],        /*	Index of the first cash-flow to be exercised */
    int ex_lncpn[],      /*	Number of coupons in each caplet */
    double ex_lstrike[], /*	Strikes for diagonal */
    double ex_lprice[],  /*	Market prices for diagonal */
    int ex_sncpn[],      /*	Number of coupons in each caplet */
    double ex_sstrike[], /*	Strikes for cap */
    double ex_sprice[],  /*	Market price for cap */
    double ex_svega[],   /*	Market vega for cap */
    double ex_zeta[],    /*	Output: zetas */
    int fix_lambda,      /*	0: calib lambda to cap        , 1: fix lambda calib
                                           to diagonal */
    int nlam,            /*	Lambda TS: may NOT be changed in the process */
    double lam_time[], double lam[], int one2F, /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    double alpha, double gamma, double rho, int skip_last,
    DIAG_CALIB_LM_PARAMS lm_params) /*	If 1        , the last option is
                                       disregarded and the forward
                                       volatility is flat from option n-1 */
{
  int i;
  double ex_spriceres[MAX_CPN];
  Err err = NULL;
  double data[MAX_CPN], weight[MAX_CPN], target[MAX_CPN], lambda[MAX_CPN];

  double fitting_error;

  if (!fix_lambda && nex < 2) {
    return serror("Cannot calibrate lambda and zeta with less than 2 exercises "
                  "- choose fix lambda");
  }

  if (fix_lambda) {
    err = lgmprcapgivenlambda_tauts2(
        ncpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time, ex_cpn, ex_lncpn,
        ex_lstrike, ex_lprice, NULL, ex_sncpn, ex_sstrike, ex_zeta, one2F, nlam,
        lam_time, lam, alpha, gamma, rho, skip_last, CALPRES, 0, MAX_FACT,
        MAX_FACT, 0, 0, ex_spriceres);

    return err;
  } else {
    for (i = 1; i <= nlam; i++) {
      lambda[i] = lam[i - 1];
    }

    lm_params->shift_freq = (int)(fmod(lm_params->shift_freq, nex) + 0.5);

    if (lm_params->use_moment) {
      /* Momentum method */
      lambda[0] = lm_params->nb_moment;

      for (i = 1; i <= lm_params->nb_moment; i++) {
        data[i] = i - 1;
        if (lm_params->vega_weight) {
          weight[i] = (i * 1.0);
        } else {
          weight[i] = 1.0;
        }

        target[i] = 0.0;
      }

      init_static_lgmprcapgiventauts(ncpn, cpn_time, cpn_df, cpn_cvg, nex,
                                     ex_time, ex_cpn, ex_lncpn, ex_lstrike,
                                     ex_lprice, ex_sncpn, ex_sstrike, ex_sprice,
                                     ex_zeta, one2F, nlam, lam_time, lambda,
                                     alpha, gamma, rho, skip_last, lm_params);

      err = levenberg_marquardt(data, target, weight, lm_params->nb_moment,
                                lambda, nlam, lm_params->nb_iter,
                                static_lgmprcapgiventauts_momentum,
                                &fitting_error);

      if (err) {
        return err;
      }
    } else {
      for (i = 1; i <= nex; i++) {
        data[i] = (int)(lm_params->shift_freq +
                        (i - 1) * lm_params->freq_short + 0.5);
        if (lm_params->vega_weight) {
          weight[i] = ex_svega[(int)data[i]];
        } else {
          weight[i] = 1.0;
        }
        target[i] = ex_sprice[(int)data[i]];
      }

      init_static_lgmprcapgiventauts(
          ncpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time, ex_cpn, ex_lncpn,
          ex_lstrike, ex_lprice, ex_sncpn, ex_sstrike, ex_sprice, ex_zeta,
          one2F, nlam, lam_time, lam, alpha, gamma, rho, skip_last, lm_params);

      err = levenberg_marquardt(
          data, target, weight, nex / lm_params->freq_short, lambda, nlam,
          lm_params->nb_iter, static_lgmprcapgiventauts, &fitting_error);
    }

    for (i = 1; i <= nlam; i++) {
      lam[i - 1] = lambda[i];
      ex_zeta[i] = EX_ZETA_LM[i];
    }
  }

  return NULL;
}

void cpd_init_hermite_for_calib(int one2F) {
  if (one2F == 2) {
    HermiteStandard(x, w, NUM_HERMITE);
  }
}

/*	Calibrate lgm: main function */
Err cpd_calib_diagonal(
    char *yc_name,        /*	Name of the yield curve */
    char *vol_curve_name, /*	Name of the market vol curve */
    char *ref_rate_name,  /*	Name of the reference rate */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    double vol_shift, int shift_type, /*	0:	Additive
                                                1:	Multiplicative */
                                      /*	If ex_date is NULL        ,
                                      exercise dates will be generated 2bd before start */
    int num_ex_dates,                 /*	Exercise dates        ,
                                                                all supposed to be on or
                                   after today */
    long *ex_date,                    /*	Supposed to be sorted
                                                                NULL = 2bd before each coupon
                                 */
    long end_date,                    /*	End date for diagonal */
    double *long_strike,              /*	Diagonal swaption strikes
                                                                        NULL = ATM */
    double *short_strike,             /*	Short swaption strikes
                                                                NULL = ATM */
    int strike_type,                  /*	0: ATM
                                                        1: CASH
                                                        2: SWAP
                                                        3: STD */
    double max_std_long, double max_std_short,
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis, int fix_lambda, /*	0: calib lambda to cap        ,
                                 1: fix lambda calib to diagonal */
    int one_f_equi,                       /*	1F equivalent flag:
                                                            if set to 1        , then 2F lambda will
                                       calibrate                       to the cap priced within
                                       calibrated 1F                       with the given lambda */
    int skip_last,     /*	If 1        , the last option is disregarded
                                         and the forward volatility is flat
                    from option                 n-1 */
    double long_prec,  /*	Precision on primary instruments */
    double short_prec, /*	Precision on secondary instruments */
    double min_fact,   /*	Maximum down jump on variance */
    double max_fact,   /*	Maximum up jump on variance */
    int use_jumps,     /*	Allow vol term structure to jump */
    int proba_weight,  /*	Proba weighting for caplet */
    double *proba, double *lambda, /*	Lambda: may be changed in the process */
    int one2F,                     /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    double alpha, double gamma, double rho, int *num_sig, /*	Answer */
    double **sig_time, double **sig,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA
        inst_data) /*	NULL = don't save calibration instrument data */
{
  int i, j, k, l, nex, ncpn;
  SrtCompounding ifreq;
  SrtBasisCode ibasis;
  double ex_time[MAX_CPN], ex_lfwd[MAX_CPN], ex_llvl[MAX_CPN],
      ex_lstrike[MAX_CPN], ex_lvol[MAX_CPN], ex_lprice[MAX_CPN],
      ex_sfwd[MAX_CPN], ex_slvl[MAX_CPN], ex_sstrike[MAX_CPN], ex_svol[MAX_CPN],
      ex_sprice[MAX_CPN], ex_sweight[MAX_CPN], ex_zeta[MAX_CPN], ex_G[MAX_CPN];
  int ex_cpn[MAX_CPN], ex_sncpn[MAX_CPN];
  long cpn_date[MAX_CPN];
  long cpn_theo_date[MAX_CPN];
  double cpn_time[MAX_CPN], cpn_cvg[MAX_CPN], cpn_df[MAX_CPN], cpn_G[MAX_CPN];
  long tmplarr[MAX_CPN];
  long theo_date, act_date, temp_date;
  long today;
  double lvl, dfi, dff;
  double power;
  double swp_rte, spr;
  double cap_price;
  double std;
  SrtCurvePtr yc_ptr;
  Err err = NULL;

  *sig_time = NULL;
  *sig = NULL;

  yc_ptr = lookup_curve(yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }
  today = get_today_from_curve(yc_ptr);

  if (one2F == 2) {
    HermiteStandard(x, w, NUM_HERMITE);
  }

  /* default for min_fact and max_fact */
  if (fabs(min_fact) < 1.0E-08) {
    min_fact = MAX_FACT;
  }

  if (fabs(max_fact) < 1.0E-08) {
    max_fact = MAX_FACT / 4.0;
  }

  /*	1.)	Setup the bond schedule and its coupons */

  /*	Coupons */

  err = interp_compounding(swaption_freq, &ifreq);
  if (err) {
    goto FREE_RETURN;
  }

  err = interp_basis(swaption_basis, &ibasis);
  if (err) {
    goto FREE_RETURN;
  }

  theo_date = end_date;
  act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
  ncpn = 1;

  while (act_date > today) {
    theo_date =
        add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
    act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    ncpn++;
  }
  ncpn--;

  if (ncpn < 2) {
    err = "Not enough coupons in cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  theo_date = end_date;
  act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
  i = ncpn - 1;

  while (i >= 0) {
    cpn_time[i] = (act_date - today) * YEARS_IN_DAY;
    cpn_date[i] = act_date;
    cpn_theo_date[i] = theo_date;

    theo_date =
        add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

    temp_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    cpn_cvg[i] = coverage(temp_date, act_date, ibasis);
    cpn_df[i] = swp_f_df(today, act_date, yc_name);
    act_date = temp_date;

    i--;
  }
  cpn_cvg[0] = 0.0;

  /*	Exercise */

  if (!ex_date) {
    num_ex_dates = ncpn - 1;

    ex_date = tmplarr;

    for (i = 0; i < num_ex_dates; i++) {
      ex_date[i] = add_unit(cpn_date[i], -2, SRT_BDAY, MODIFIED_SUCCEEDING);
    }
  } else {
    memcpy(tmplarr, ex_date, num_ex_dates * sizeof(long));
    ex_date = tmplarr;
  }

  /*	Remove past dates */
  while (num_ex_dates && ex_date[0] <= today) {
    ex_date++;
    if (long_strike)
      long_strike++;
    if (short_strike)
      short_strike++;
    if (proba)
      proba++;
    num_ex_dates--;
  }

  /*	Remove redundant dates */
  j = ncpn - 1;
  l = ncpn + 1;
  for (i = num_ex_dates - 1; i >= 0; i--) {
    while (j > 0 && cpn_date[j] > ex_date[i]) {
      j--;
    }
    if (cpn_date[j] < ex_date[i]) {
      j++;
    }

    if (j >= ncpn - 1 || j == l) {
      for (k = i - 1; k >= 0; k--) {
        ex_date[k + 1] = ex_date[k];
        if (long_strike)
          long_strike[k + 1] = long_strike[k];
        if (short_strike)
          short_strike[k + 1] = short_strike[k];
        if (proba)
          proba[k + 1] = proba[k];
      }

      ex_date++;
      if (long_strike)
        long_strike++;
      if (short_strike)
        short_strike++;
      if (proba)
        proba++;
      num_ex_dates--;
    } else {
      l = j;
    }
  }

  if (num_ex_dates < 1) {
    err = "All exercise dates are past in cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  for (i = 0; i < num_ex_dates; i++) {
    if (proba_weight && proba) {
      ex_sweight[i] = proba[i];
    } else {
      ex_sweight[i] = 1.0;
    }
  }

  nex = num_ex_dates;
  j = 0;
  for (i = 0; i < nex; i++) {
    while (cpn_date[j] < ex_date[i]) {
      j++;
    }

    ex_cpn[i] = j;
    ex_time[i] = (ex_date[i] - today) * YEARS_IN_DAY;
  }

  /*	Underlyings */

  /*	Long */

  dff = swp_f_df(today, cpn_date[ncpn - 1], yc_name);
  for (i = 0; i < nex; i++) {
    j = ex_cpn[i];

    lvl = 0.0;
    for (k = j + 1; k < ncpn; k++) {
      lvl += cpn_cvg[k] * cpn_df[k];
    }
    dfi = swp_f_df(today, cpn_date[j], yc_name);

    ex_llvl[i] = lvl;
    ex_lfwd[i] = (dfi - dff) / lvl;

    /*	ATM std */
    err = get_cash_vol(vol_curve_name,
                       add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                       end_date, ex_lfwd[i], 0, ref_rate_name, &std, &power);
    if (err) {
      goto FREE_RETURN;
    }
    std += (shift_type == 1 ? std * vol_shift : vol_shift);
    if (power > 0.5) {
      power = srt_f_optblksch(ex_lfwd[i], ex_lfwd[i], std, ex_time[i], 1.0,
                              SRT_CALL, PREMIUM);
      err = srt_f_optimpvol(power, ex_lfwd[i], ex_lfwd[i], ex_time[i], 1.0,
                            SRT_CALL, SRT_NORMAL, &std);
    }
    std *= sqrt(ex_time[i]);

    /*	Strike */
    if ((!long_strike) || (!strike_type)) {
      ex_lstrike[i] = ex_lfwd[i];
    } else if (strike_type == 1) {
      ex_lstrike[i] = long_strike[i];
    } else if (strike_type == 2) {
      if (err = swp_f_ForwardRate(cpn_date[j], end_date, swaption_freq,
                                  swaption_basis, yc_name, ref_rate_name,
                                  &swp_rte)) {
        goto FREE_RETURN;
      }

      spr = swp_rte - ex_lfwd[i];

      ex_lstrike[i] = long_strike[i] - spr;
    } else if (strike_type == 3) {
      ex_lstrike[i] = ex_lfwd[i] + long_strike[i] * std;
    }

    /*	Apply max std */
    if (ex_lstrike[i] > ex_lfwd[i] + max_std_long * std) {
      ex_lstrike[i] = ex_lfwd[i] + max_std_long * std;
    } else if (ex_lstrike[i] < ex_lfwd[i] - max_std_long * std) {
      ex_lstrike[i] = ex_lfwd[i] - max_std_long * std;
    }

    /*	Make sure strikes are positive (actually more than 1bp)
                    otherwise use ATM	*/
    if (ex_lstrike[i] < 1.0e-04) {
      ex_lstrike[i] = ex_lfwd[i];
    }

    err = get_cash_vol(
        vol_curve_name, add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
        end_date, ex_lstrike[i], 0, ref_rate_name, &(ex_lvol[i]), &power);
    if (err) {
      goto FREE_RETURN;
    }
    ex_lvol[i] += (shift_type == 1 ? ex_lvol[i] * vol_shift : vol_shift);

    if (power > 0.5) {
      ex_lprice[i] = srt_f_optblksch(ex_lfwd[i], ex_lstrike[i], ex_lvol[i],
                                     ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
    } else {
      ex_lprice[i] = srt_f_optblknrm(ex_lfwd[i], ex_lstrike[i], ex_lvol[i],
                                     ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
    }
  }

  /*	Short */

  if (!fix_lambda) {
    cap_price = 0.0;
    for (i = 0; i < nex; i++) {
      if (i < nex - 1) {
        ex_sncpn[i] = ex_cpn[i + 1] - ex_cpn[i] + 1;
      } else {
        ex_sncpn[i] = ncpn - ex_cpn[i];
      }

      if (ex_sncpn[i] < 2) {
        err = "One exercise date controls less than 2 coupons in "
              "cpd_calib_diagonal";
        goto FREE_RETURN;
      }

      lvl = 0.0;
      for (k = ex_cpn[i] + 1; k < ex_cpn[i] + ex_sncpn[i]; k++) {
        lvl += cpn_cvg[k] * cpn_df[k];
      }
      dfi = swp_f_df(today, cpn_date[ex_cpn[i]], yc_name);
      dff = swp_f_df(today, cpn_date[ex_cpn[i] + ex_sncpn[i] - 1], yc_name);

      ex_slvl[i] = lvl;
      ex_sfwd[i] = (dfi - dff) / lvl;

      /*	ATM std */
      err = get_cash_vol(vol_curve_name,
                         add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                         cpn_theo_date[ex_cpn[i] + ex_sncpn[i] - 1], ex_sfwd[i],
                         0, ref_rate_name, &std, &power);
      if (err) {
        goto FREE_RETURN;
      }
      std += (shift_type == 1 ? std * vol_shift : vol_shift);
      if (power > 0.5) {
        power = srt_f_optblksch(ex_sfwd[i], ex_sfwd[i], std, ex_time[i], 1.0,
                                SRT_CALL, PREMIUM);
        err = srt_f_optimpvol(power, ex_sfwd[i], ex_sfwd[i], ex_time[i], 1.0,
                              SRT_CALL, SRT_NORMAL, &std);
      }
      std *= sqrt(ex_time[i]);

      /*	Strike */
      if ((!short_strike) || (!strike_type)) {
        ex_sstrike[i] = ex_sfwd[i];
      } else if (strike_type == 1) {
        ex_sstrike[i] = short_strike[i];
      } else if (strike_type == 2) {
        if (err = swp_f_ForwardRate(cpn_date[ex_cpn[i]],
                                    cpn_theo_date[ex_cpn[i] + ex_sncpn[i] - 1],
                                    swaption_freq, swaption_basis, yc_name,
                                    ref_rate_name, &swp_rte)) {
          goto FREE_RETURN;
        }

        spr = swp_rte - ex_sfwd[i];

        ex_sstrike[i] = short_strike[i] - spr;
      } else if (strike_type == 3) {
        ex_sstrike[i] = ex_sfwd[i] + short_strike[i] * std;
      }

      /*	Apply max std */
      if (ex_sstrike[i] > ex_sfwd[i] + max_std_short * std) {
        ex_sstrike[i] = ex_sfwd[i] + max_std_short * std;
      } else if (ex_sstrike[i] < ex_sfwd[i] - max_std_short * std) {
        ex_sstrike[i] = ex_sfwd[i] - max_std_short * std;
      }

      /*	Make sure strikes are positive (actually more than 1bp)
                      otherwise use ATM	*/
      if (ex_sstrike[i] < 1.0e-04) {
        ex_sstrike[i] = ex_sfwd[i];
      }

      err =
          get_cash_vol(vol_curve_name,
                       add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                       cpn_theo_date[ex_cpn[i] + ex_sncpn[i] - 1],
                       ex_sstrike[i], 0, ref_rate_name, &(ex_svol[i]), &power);
      if (err) {
        goto FREE_RETURN;
      }
      ex_svol[i] += (shift_type == 1 ? ex_svol[i] * vol_shift : vol_shift);

      if (power > 0.5) {
        ex_sprice[i] =
            srt_f_optblksch(ex_sfwd[i], ex_sstrike[i], ex_svol[i], ex_time[i],
                            ex_slvl[i], SRT_PUT, PREMIUM);
      } else {
        ex_sprice[i] =
            srt_f_optblknrm(ex_sfwd[i], ex_sstrike[i], ex_svol[i], ex_time[i],
                            ex_slvl[i], SRT_PUT, PREMIUM);
      }

      cap_price += ex_sweight[i] * ex_sprice[i];
    }
  }

  /*	The 1F equivalent case */

  if (one2F == 2 && fix_lambda && one_f_equi) {
    cap_price = 0.0;
    for (i = 0; i < nex; i++) {
      if (i < nex - 1) {
        ex_sncpn[i] = ex_cpn[i + 1] - ex_cpn[i] + 1;
      } else {
        ex_sncpn[i] = ncpn - ex_cpn[i];
      }

      if (ex_sncpn[i] < 2) {
        err = "One exercise date controls less than 2 coupons in "
              "cpd_calib_diagonal";
        goto FREE_RETURN;
      }

      lvl = 0.0;
      for (k = ex_cpn[i] + 1; k < ex_cpn[i] + ex_sncpn[i]; k++) {
        lvl += cpn_cvg[k] * cpn_df[k];
      }
      dfi = swp_f_df(today, cpn_date[ex_cpn[i]], yc_name);
      dff = swp_f_df(today, cpn_date[ex_cpn[i] + ex_sncpn[i] - 1], yc_name);

      ex_slvl[i] = lvl;
      ex_sfwd[i] = (dfi - dff) / lvl;
      ex_sstrike[i] = ex_sfwd[i];
    }

    err = lgmcalibzetalambda(ncpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time,
                             ex_cpn, ex_sncpn, ex_lstrike, ex_lprice,
                             ex_sstrike, ex_sweight, 0.0, ex_zeta, 1, lambda, 1,
                             0.0, 0.0, 0.0, skip_last, long_prec, short_prec,
                             min_fact, max_fact, use_jumps);

    if (err) {
      goto FREE_RETURN;
    }

    static_lgmsetupG(*lambda, ncpn, cpn_time, cpn_G, nex, ex_time, ex_G);

    cap_price = lgmcapval1F(ncpn, cpn_df, cpn_cvg, cpn_G, nex, ex_cpn, ex_sncpn,
                            ex_sweight, ex_zeta, ex_G, ex_sstrike);

    fix_lambda = 0;
  }

  /*	2.)	Calibrate lambda and zeta */

  err = lgmcalibzetalambda(ncpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time,
                           ex_cpn, ex_sncpn, ex_lstrike, ex_lprice, ex_sstrike,
                           ex_sweight, cap_price, ex_zeta, fix_lambda, lambda,
                           one2F, alpha, gamma, rho, skip_last, long_prec,
                           short_prec, min_fact, max_fact, use_jumps);

  if (err) {
    goto FREE_RETURN;
  }

  /*	3.)	Transform into sigma */

  *num_sig = nex;
  *sig_time = (double *)calloc(nex, sizeof(double));
  *sig = (double *)calloc(nex, sizeof(double));

  if (!sig_time || !sig) {
    err = "Allocation error (3) in cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  (*sig_time)[0] = ex_time[0];
  (*sig)[0] = sqrt(ex_zeta[0] * 2 * (*lambda) /
                   (exp(2 * (*lambda) * ex_time[0]) - 1.0));

  for (i = 1; i < nex; i++) {
    (*sig_time)[i] = ex_time[i];
    if (ex_zeta[i] > ex_zeta[i - 1]) {
      (*sig)[i] = sqrt((ex_zeta[i] - ex_zeta[i - 1]) * 2 * (*lambda) /
                       (exp(2 * (*lambda) * ex_time[i]) -
                        exp(2 * (*lambda) * ex_time[i - 1])));
    } else {
      smessage("Diagonal calibration failed at exercise year %.2f - "
               "Calibration stopped",
               ex_time[i]);
      for (j = i; j < nex; j++) {
        (*sig)[j] = (*sig)[i - 1];
      }
      i = nex;
    }
  }

  /*	4.)	Save instrument data if required */
  if (inst_data) {
    inst_data->num_inst = nex;
    inst_data->num_insts = nex;
    inst_data->start_dates = (long *)calloc(nex, sizeof(long));
    inst_data->end_dates = (long *)calloc(nex, sizeof(long));
    inst_data->short_strikes = (double *)calloc(nex, sizeof(double));
    inst_data->long_strikes = (double *)calloc(nex, sizeof(double));
    inst_data->short_weights = (double *)calloc(nex, sizeof(double));

    if (!inst_data->start_dates || !inst_data->end_dates ||
        !inst_data->short_strikes || !inst_data->short_weights ||
        !inst_data->long_strikes) {
      err = "Allocation error (4) in cpd_calib_diagonal";
      goto FREE_RETURN;
    }

    for (i = 0; i < nex; i++) {
      inst_data->start_dates[i] = cpn_date[ex_cpn[i]];
      inst_data->end_dates[i] = cpn_date[ncpn - 1];
      if (!fix_lambda) {
        inst_data->short_strikes[i] = ex_sstrike[i];
        inst_data->short_weights[i] = ex_sweight[i];
      } else {
        inst_data->short_strikes[i] = 0.0;
        inst_data->short_weights[i] = 1.0;
      }
      inst_data->long_strikes[i] = ex_lstrike[i];
    }
  }

FREE_RETURN:

  if (err) {
    if (*sig_time)
      free(*sig_time);
    *sig_time = NULL;

    if (*sig)
      free(*sig);
    *sig = NULL;

    /*
    if (inst_data)
    {
            cpd_free_calib_inst_data (inst_data);
    }
    */
  }

  return err;
}

/*	Calibrate lgm: main function */
/*	Tau-TS enabled version
 *** no lambda calibration so far for this version *** */
Err cpd_calib_diagonal_tauts(
    char *yc_name,        /*	Name of the yield curve */
    char *vol_curve_name, /*	Name of the market vol curve */
    char *ref_rate_name,  /*	Name of the reference rate */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    double vol_shift, int shift_type, /*	0:	Additive
                                                1:	Multiplicative */
                                      /*	If ex_date is NULL        ,
                                      exercise dates will be generated 2bd before start */
    int num_ex_dates,                 /*	Exercise dates        ,
                                                                all supposed to be on or
                                   after today */
    long *ex_date,                    /*	Supposed to be sorted
                                                                NULL = 2bd before each coupon
                                 */
    long end_date,                    /*	End date for diagonal */
    double *long_strike,              /*	Diagonal swaption strikes
                                                                        NULL = ATM */
    int strike_type,                  /*	0: ATM
                                                        1: CASH
                                                        2: SWAP
                                                        3: STD */
    double max_std_long,
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,
    int skip_last,   /*	If 1        , the last option is disregarded and the
                  forward   volatility is flat from option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps,   /*	Allow vol term structure to jump */
    int nlam,        /*	Lambda TS: may NOT be changed in the process */
    double lam_time[], double lam[], int one2F, /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    double alpha, double gamma, double rho, int *num_sig, /*	Answer */
    double **sig_time, double **sig,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA
        inst_data) /*	NULL = don't save calibration instrument data */
{
  int i, j, k, l, nex, ncpn;
  SrtCompounding ifreq;
  SrtBasisCode ibasis;
  double ex_time[MAX_CPN], ex_lfwd[MAX_CPN], ex_llvl[MAX_CPN],
      ex_lstrike[MAX_CPN], ex_lvol[MAX_CPN], ex_lprice[MAX_CPN],
      ex_zeta[MAX_CPN];
  int ex_cpn[MAX_CPN];
  long cpn_date[MAX_CPN];
  double cpn_time[MAX_CPN], cpn_cvg[MAX_CPN], cpn_df[MAX_CPN];
  long tmplarr[MAX_CPN];
  long theo_date, act_date, temp_date;
  long today;
  double lvl, dfi, dff;
  double power;
  double swp_rte, spr;
  double std;
  double exp_fact;
  SrtCurvePtr yc_ptr;
  Err err = NULL;

  *sig_time = NULL;
  *sig = NULL;

  yc_ptr = lookup_curve(yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }
  today = get_today_from_curve(yc_ptr);

  if (one2F == 2) {
    HermiteStandard(x, w, NUM_HERMITE);
  }

  /* default for min_fact and max_fact */
  if (fabs(min_fact) < 1.0E-08) {
    min_fact = MAX_FACT;
  }

  if (fabs(max_fact) < 1.0E-08) {
    max_fact = MAX_FACT / 4.0;
  }

  /*	1.)	Setup the bond schedule and its coupons */

  /*	Coupons */

  err = interp_compounding(swaption_freq, &ifreq);
  if (err) {
    goto FREE_RETURN;
  }

  err = interp_basis(swaption_basis, &ibasis);
  if (err) {
    goto FREE_RETURN;
  }

  theo_date = end_date;
  act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
  ncpn = 1;

  while (act_date > today) {
    theo_date =
        add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
    act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    ncpn++;
  }
  ncpn--;

  if (ncpn < 2) {
    err = "Not enough coupons in cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  theo_date = end_date;
  act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
  i = ncpn - 1;

  while (i >= 0) {
    cpn_time[i] = (act_date - today) * YEARS_IN_DAY;
    cpn_date[i] = act_date;

    theo_date =
        add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

    temp_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    cpn_cvg[i] = coverage(temp_date, act_date, ibasis);
    cpn_df[i] = swp_f_df(today, act_date, yc_name);
    act_date = temp_date;

    i--;
  }
  cpn_cvg[0] = 0.0;

  /*	Exercise */

  if (!ex_date) {
    num_ex_dates = ncpn - 1;

    ex_date = tmplarr;

    for (i = 0; i < num_ex_dates; i++) {
      ex_date[i] = add_unit(cpn_date[i], -2, SRT_BDAY, MODIFIED_SUCCEEDING);
    }
  } else {
    memcpy(tmplarr, ex_date, num_ex_dates * sizeof(long));
    ex_date = tmplarr;
  }

  /*	Remove past dates */
  while (num_ex_dates && ex_date[0] <= today) {
    ex_date++;
    if (long_strike)
      long_strike++;
    num_ex_dates--;
  }

  /*	Remove redundant dates */
  j = ncpn - 1;
  l = ncpn + 1;
  for (i = num_ex_dates - 1; i >= 0; i--) {
    while (j > 0 && cpn_date[j] > ex_date[i]) {
      j--;
    }
    if (cpn_date[j] < ex_date[i]) {
      j++;
    }

    if (j >= ncpn - 1 || j == l) {
      for (k = i - 1; k >= 0; k--) {
        ex_date[k + 1] = ex_date[k];
        if (long_strike)
          long_strike[k + 1] = long_strike[k];
      }

      ex_date++;
      if (long_strike)
        long_strike++;
      num_ex_dates--;
    } else {
      l = j;
    }
  }

  if (num_ex_dates < 1) {
    err = "All exercise dates are past in cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  nex = num_ex_dates;
  j = 0;
  for (i = 0; i < nex; i++) {
    while (cpn_date[j] < ex_date[i]) {
      j++;
    }

    ex_cpn[i] = j;
    ex_time[i] = (ex_date[i] - today) * YEARS_IN_DAY;
  }

  /*	Underlyings */

  /*	Long */

  dff = swp_f_df(today, cpn_date[ncpn - 1], yc_name);
  for (i = 0; i < nex; i++) {
    j = ex_cpn[i];

    lvl = 0.0;
    for (k = j + 1; k < ncpn; k++) {
      lvl += cpn_cvg[k] * cpn_df[k];
    }
    dfi = swp_f_df(today, cpn_date[j], yc_name);

    ex_llvl[i] = lvl;
    ex_lfwd[i] = (dfi - dff) / lvl;

    /*	ATM std */
    err = get_cash_vol(vol_curve_name,
                       add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                       end_date, ex_lfwd[i], 0, ref_rate_name, &std, &power);
    if (err) {
      goto FREE_RETURN;
    }
    std += (shift_type == 1 ? std * vol_shift : vol_shift);
    if (power > 0.5) {
      power = srt_f_optblksch(ex_lfwd[i], ex_lfwd[i], std, ex_time[i], 1.0,
                              SRT_CALL, PREMIUM);
      err = srt_f_optimpvol(power, ex_lfwd[i], ex_lfwd[i], ex_time[i], 1.0,
                            SRT_CALL, SRT_NORMAL, &std);
    }
    std *= sqrt(ex_time[i]);

    /*	Strike */
    if ((!long_strike) || (!strike_type)) {
      ex_lstrike[i] = ex_lfwd[i];
    } else if (strike_type == 1) {
      ex_lstrike[i] = long_strike[i];
    } else if (strike_type == 2) {
      if (err = swp_f_ForwardRate(cpn_date[j], end_date, swaption_freq,
                                  swaption_basis, yc_name, ref_rate_name,
                                  &swp_rte)) {
        goto FREE_RETURN;
      }

      spr = swp_rte - ex_lfwd[i];

      ex_lstrike[i] = long_strike[i] - spr;
    } else if (strike_type == 3) {
      ex_lstrike[i] = ex_lfwd[i] + long_strike[i] * std;
    }

    /*	Apply max std */
    if (ex_lstrike[i] > ex_lfwd[i] + max_std_long * std) {
      ex_lstrike[i] = ex_lfwd[i] + max_std_long * std;
    } else if (ex_lstrike[i] < ex_lfwd[i] - max_std_long * std) {
      ex_lstrike[i] = ex_lfwd[i] - max_std_long * std;
    }

    /*	Make sure strikes are positive (actually more than 1bp)
                    otherwise use ATM	*/
    if (ex_lstrike[i] < 1.0e-04) {
      ex_lstrike[i] = ex_lfwd[i];
    }

    err = get_cash_vol(
        vol_curve_name, add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
        end_date, ex_lstrike[i], 0, ref_rate_name, &(ex_lvol[i]), &power);
    if (err) {
      goto FREE_RETURN;
    }
    ex_lvol[i] += (shift_type == 1 ? ex_lvol[i] * vol_shift : vol_shift);

    if (power > 0.5) {
      ex_lprice[i] = srt_f_optblksch(ex_lfwd[i], ex_lstrike[i], ex_lvol[i],
                                     ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
    } else {
      ex_lprice[i] = srt_f_optblknrm(ex_lfwd[i], ex_lstrike[i], ex_lvol[i],
                                     ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
    }
  }

  /*	2.)	Calibrate zeta */

  err = lgmprcapgivenlambda_tauts(ncpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time,
                                  ex_cpn, ex_lstrike, ex_lprice, ex_zeta, one2F,
                                  nlam, lam_time, lam, alpha, gamma, rho,
                                  skip_last, min_fact, max_fact, use_jumps);

  if (err) {
    goto FREE_RETURN;
  }

  /*	3.)	Transform into sigma */

  *num_sig = nex;
  *sig_time = (double *)calloc(nex, sizeof(double));
  *sig = (double *)calloc(nex, sizeof(double));

  if (!sig_time || !sig) {
    err = "Allocation error (3) in cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  (*sig_time)[0] = ex_time[0];
  exp_fact = static_lgmcalcexpfact_tauts(0.0, ex_time[0], nlam, lam_time, lam);
  (*sig)[0] = sqrt(ex_zeta[0] / exp_fact);

  for (i = 1; i < nex; i++) {
    (*sig_time)[i] = ex_time[i];
    if (ex_zeta[i] > ex_zeta[i - 1]) {
      exp_fact = static_lgmcalcexpfact_tauts(ex_time[i - 1], ex_time[i], nlam,
                                             lam_time, lam);
      (*sig)[i] = sqrt((ex_zeta[i] - ex_zeta[i - 1]) / exp_fact);
    } else {
      smessage("Diagonal calibration failed at exercise year %.2f - "
               "Calibration stopped",
               ex_time[i]);
      for (j = i; j < nex; j++) {
        (*sig)[j] = (*sig)[i - 1];
      }
      i = nex;
    }
  }

  /*	4.)	Save instrument data if required */
  if (inst_data) {
    inst_data->num_inst = nex;
    inst_data->num_insts = nex;
    inst_data->start_dates = (long *)calloc(nex, sizeof(long));
    inst_data->end_dates = (long *)calloc(nex, sizeof(long));
    inst_data->short_strikes = (double *)calloc(nex, sizeof(double));
    inst_data->long_strikes = (double *)calloc(nex, sizeof(double));

    if (!inst_data->start_dates || !inst_data->end_dates ||
        !inst_data->short_strikes || !inst_data->long_strikes) {
      err = "Allocation error (4) in cpd_calib_diagonal";
      goto FREE_RETURN;
    }

    for (i = 0; i < nex; i++) {
      inst_data->start_dates[i] = cpn_date[ex_cpn[i]];
      inst_data->end_dates[i] = cpn_date[ncpn - 1];
      inst_data->short_strikes[i] = 0.0;
      inst_data->long_strikes[i] = ex_lstrike[i];
    }
  }

FREE_RETURN:

  if (err) {
    if (*sig_time)
      free(*sig_time);
    *sig_time = NULL;

    if (*sig)
      free(*sig);
    *sig = NULL;

    /*
    if (inst_data)
    {
            cpd_free_calib_inst_data (inst_data);
    }
    */
  }

  return err;
}

Err get_end_date(long ex_date, long struct_end_date, char *tenor, int theo_act,
                 long *end_date) {
  long start_date;
  Err err;

  start_date = add_unit(ex_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
  strupper(tenor);
  strip_white_space(tenor);
  if (*tenor == 'D') {
    *end_date = struct_end_date;
  } else {
    err = add_tenor(start_date, tenor, NO_BUSDAY_CONVENTION, end_date);
    if (err) {
      return err;
    }
  }

  if (theo_act) {
    *end_date = bus_date_method(*end_date, MODIFIED_SUCCEEDING);
  }

  return NULL;
}

#define ONE_MONTH 0.083333333
/*	Calibrate lgm: main function */
/*	New version: calibrates not necessarily to digonal
 *** no lambda calibration so far for this version *** */
Err cpd_calib_diagonal_2(
    /*	Market */
    char *yc_name,        /*	Name of the yield curve */
    char *vol_curve_name, /*	Name of the market vol curve */
    char *ref_rate_name,  /*	Name of the reference rate */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    char *instr_freq, /*	Frequency and basis of instruments */
    char *instr_basis,
    /*	If ex_date is NULL        ,
    exercise dates will be generated 2bd before start */
    /*	Structure */
    int num_ex_dates, /*	Exercise dates        ,
                                                all supposed to be on or
                   after today */
    long *ex_date_,   /*	Supposed to be sorted */
    int *cal_date,    /*	1: use ex_date as calibration date        , 0: don't */
    char **end_tenor_, /*	Tenors of the underlying instruments
                                                         or "DIAG" */
    long end_date,     /*	End date for diagonal */
    double *strike_,   /*	Strikes
                                         0: ATM */
    /*	Model */
    double lambda, /*	Lambda: may NOT be changed in the process */
    int one2F,     /*	Number of factors */
    double alpha,  /*	Alpha        , Gamma        , Rho (2F only) */
    double gamma, double rho,
    /*	Output */
    int *num_sig, /*	Answer */
    double **sig_time, double **sig,
    /*	Parameters */
    CPD_DIAG_CALIB_PARAM param,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA
        inst_data) /*	NULL = don't save calibration instrument data */
{
  int i, j, k, l, nex, ncpn;
  SrtCompounding ifreq;
  SrtBasisCode ibasis;
  double ex_time[MAX_CPN], ex_lfwd[MAX_CPN], ex_llvl[MAX_CPN],
      ex_lstrike[MAX_CPN], ex_lvol[MAX_CPN], ex_lprice[MAX_CPN],
      ex_zeta[MAX_CPN];
  int ex_cpn[MAX_CPN], ex_endcpn[MAX_CPN];
  long cpn_date[MAX_CPN];
  double cpn_time[MAX_CPN], cpn_cvg[MAX_CPN], cpn_df[MAX_CPN];
  long tmplng1[MAX_CPN], tmplng2[MAX_CPN];
  long *theo_end_dates, *act_end_dates;
  long theo_date, act_date, temp_date, temp_date2;
  long today;
  double lvl, dfi, dff;
  double power;
  double swp_rte, spr;
  double std;
  SrtCurvePtr yc_ptr;
  Err err = NULL;
  long ex_date__[MAX_INST], *ex_date;
  double strike__[MAX_INST], *strike;
  char end_tenor__[MAX_INST * 256];
  String end_tenor___[MAX_INST], *end_tenor;

  /*	Copy data so as not to change the original */
  ex_date = &(ex_date__[0]);
  strike = &(strike__[0]);
  end_tenor = &(end_tenor___[0]);
  for (i = 0; i < num_ex_dates; i++) {
    end_tenor[i] = &(end_tenor__[i * 256]);
  }
  memcpy(ex_date, ex_date_, num_ex_dates * sizeof(long));
  memcpy(strike, strike_, num_ex_dates * sizeof(double));
  for (i = 0; i < num_ex_dates; i++) {
    strcpy(end_tenor[i], end_tenor_[i]);
  }

  theo_end_dates = &(tmplng1[0]);
  act_end_dates = &(tmplng2[0]);

  *sig_time = NULL;
  *sig = NULL;

  yc_ptr = lookup_curve(yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }
  today = get_today_from_curve(yc_ptr);

  if (one2F == 2) {
    HermiteStandard(x, w, NUM_HERMITE);
  }

  /*	1.)	Setup the bond schedule and its coupons */

  /*	Coupons */

  err = interp_compounding(instr_freq, &ifreq);
  if (err) {
    goto FREE_RETURN;
  }

  err = interp_basis(instr_basis, &ibasis);
  if (err) {
    goto FREE_RETURN;
  }

  /*	Find the end date as the longest total maturity */
  theo_date = end_date;
  act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
  for (i = 0; i < num_ex_dates; i++) {
    err = get_end_date(ex_date[i], end_date, end_tenor[i], 0,
                       &(theo_end_dates[i]));
    if (err) {
      goto FREE_RETURN;
    }
    act_end_dates[i] = bus_date_method(theo_end_dates[i], MODIFIED_SUCCEEDING);
  }
  for (i = 0; i < num_ex_dates; i++) {
    if (theo_end_dates[i] > theo_date || act_end_dates[i] > act_date) {
      theo_date = theo_end_dates[i];
      act_date = act_end_dates[i];
    }
  }
  ncpn = 1;
  temp_date = theo_date;
  temp_date2 = act_date;

  while (act_date > today) {
    theo_date =
        add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
    act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    ncpn++;
  }
  ncpn--;

  if (ncpn < 2) {
    err = "Not enough coupons in cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  theo_date = temp_date;
  act_date = temp_date2;
  i = ncpn - 1;

  while (i >= 0) {
    cpn_time[i] = (act_date - today) * YEARS_IN_DAY;
    cpn_date[i] = act_date;

    theo_date =
        add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

    temp_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    cpn_cvg[i] = coverage(temp_date, act_date, ibasis);
    cpn_df[i] = swp_f_df(today, act_date, yc_name);
    act_date = temp_date;

    i--;
  }
  cpn_cvg[0] = 0.0;

  /*	Exercise */

  /*	Remove non-calibration dates */
  for (i = num_ex_dates - 1; i >= 0; i--) {
    if (cal_date[i] == 0) {
      for (k = i - 1; k >= 0; k--) {
        ex_date[k + 1] = ex_date[k];
        strcpy(end_tenor[k + 1], end_tenor[k]);
        theo_end_dates[k + 1] = theo_end_dates[k];
        act_end_dates[k + 1] = act_end_dates[k];
        strike[k + 1] = strike[k];
      }

      ex_date++;
      end_tenor++;
      theo_end_dates++;
      act_end_dates++;
      strike++;
      num_ex_dates--;
      if (num_ex_dates < 1) {
        err = "All exercise dates are past in cpd_calib_diagonal";
        goto FREE_RETURN;
      }
    }
  }

  /*	Remove redundant dates */
  j = ncpn - 1;
  l = ncpn + 1;
  for (i = num_ex_dates - 1; i >= 0; i--) {
    while (j > 0 && cpn_date[j] > ex_date[i]) {
      j--;
    }
    if (cpn_date[j] < ex_date[i]) {
      j++;
    }

    if (j >= ncpn - 1 || j == l) {
      for (k = i - 1; k >= 0; k--) {
        ex_date[k + 1] = ex_date[k];
        strcpy(end_tenor[k + 1], end_tenor[k]);
        theo_end_dates[k + 1] = theo_end_dates[k];
        act_end_dates[k + 1] = act_end_dates[k];
        strike[k + 1] = strike[k];
      }

      ex_date++;
      end_tenor++;
      theo_end_dates++;
      act_end_dates++;
      strike++;
      num_ex_dates--;
      if (num_ex_dates < 1) {
        err = "All exercise dates are past in cpd_calib_diagonal";
        goto FREE_RETURN;
      }
    } else {
      l = j;
    }
  }

  /*	Remove close dates */
  j = num_ex_dates - 1;
  for (i = num_ex_dates - 2; i >= 0; i--) {
    if ((ex_date[j] - ex_date[i]) * YEARS_IN_DAY <
        param->min_time - ONE_MONTH) {
      for (k = i - 1; k >= 0; k--) {
        ex_date[k + 1] = ex_date[k];
        strcpy(end_tenor[k + 1], end_tenor[k]);
        theo_end_dates[k + 1] = theo_end_dates[k];
        act_end_dates[k + 1] = act_end_dates[k];
        strike[k + 1] = strike[k];
      }

      ex_date++;
      end_tenor++;
      theo_end_dates++;
      act_end_dates++;
      strike++;
      num_ex_dates--;
      j--;
      if (num_ex_dates < 1) {
        err = "All exercise dates are past in cpd_calib_diagonal";
        goto FREE_RETURN;
      }
    } else {
      j = i;
    }
  }

  /*	Remove last? */
  if (param->skip_last && num_ex_dates > 1) {
    num_ex_dates--;
  }

  /*	Remove past dates */
  while (ex_date[0] <= today) {
    ex_date++;
    end_tenor++;
    theo_end_dates++;
    act_end_dates++;
    strike++;
    num_ex_dates--;
    if (num_ex_dates == 0) {
      "All exercise dates are past in cpd_calib_diagonal";
      goto FREE_RETURN;
    }
  }

  nex = num_ex_dates;
  j = 0;
  for (i = 0; i < nex; i++) {
    while (cpn_date[j] < ex_date[i]) {
      j++;
    }

    ex_cpn[i] = j;
    ex_time[i] = (ex_date[i] - today) * YEARS_IN_DAY;

    k = j;
    while (cpn_date[k] < act_end_dates[i]) {
      k++;
    }
    if (k > 0 &&
        cpn_date[k] - act_end_dates[i] > act_end_dates[i] - cpn_date[k - 1]) {
      k--;
    }

    if (k <= j) {
      k = j + 1;
    }
    ex_endcpn[i] = k;

    if (j >= ncpn || k >= ncpn) {
      err = "Coupon date bug in cpd_calib_diagonal";
      goto FREE_RETURN;
    }
  }

  /*	Underlyings */

  /*	Long */

  for (i = 0; i < nex; i++) {
    j = ex_cpn[i];
    l = ex_endcpn[i];

    lvl = 0.0;
    for (k = j + 1; k <= l; k++) {
      lvl += cpn_cvg[k] * cpn_df[k];
    }
    dfi = swp_f_df(today, cpn_date[j], yc_name);
    dff = swp_f_df(today, cpn_date[l], yc_name);

    ex_llvl[i] = lvl;
    ex_lfwd[i] = (dfi - dff) / lvl;

    /*	ATM std */
    err = get_cash_vol(
        vol_curve_name, add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
        theo_end_dates[i], ex_lfwd[i], 0, ref_rate_name, &std, &power);
    if (err) {
      goto FREE_RETURN;
    }
    std += (param->shift_type == 1 ? std * param->vol_shift : param->vol_shift);
    if (power > 0.5) {
      power = srt_f_optblksch(ex_lfwd[i], ex_lfwd[i], std, ex_time[i], 1.0,
                              SRT_CALL, PREMIUM);
      err = srt_f_optimpvol(power, ex_lfwd[i], ex_lfwd[i], ex_time[i], 1.0,
                            SRT_CALL, SRT_NORMAL, &std);
    }
    std *= sqrt(ex_time[i]);

    /*	Strike */
    if (param->strike_type == 0 || strike[i] < 1.0e-04) {
      ex_lstrike[i] = ex_lfwd[i];
    } else if (param->strike_type == 1) {
      ex_lstrike[i] = strike[i];
    } else if (param->strike_type == 2) {
      if (err = swp_f_ForwardRate(cpn_date[j], theo_end_dates[i], instr_freq,
                                  instr_basis, yc_name, ref_rate_name,
                                  &swp_rte)) {
        goto FREE_RETURN;
      }

      spr = swp_rte - ex_lfwd[i];

      ex_lstrike[i] = strike[i] - spr;
    } else if (param->strike_type == 3) {
      ex_lstrike[i] = ex_lfwd[i] + strike[i] * std;
    }

    /*	Apply max std */
    if (ex_lstrike[i] > ex_lfwd[i] + param->max_std * std) {
      ex_lstrike[i] = ex_lfwd[i] + param->max_std * std;
    } else if (ex_lstrike[i] < ex_lfwd[i] - param->max_std * std) {
      ex_lstrike[i] = ex_lfwd[i] - param->max_std * std;
    }

    /*	Make sure strikes are positive (actually more than 1bp)
                    otherwise use ATM	*/
    if (ex_lstrike[i] < 1.0e-04) {
      ex_lstrike[i] = ex_lfwd[i];
    }

    err = get_cash_vol(vol_curve_name,
                       add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                       theo_end_dates[i], ex_lstrike[i], 0, ref_rate_name,
                       &(ex_lvol[i]), &power);
    if (err) {
      goto FREE_RETURN;
    }
    ex_lvol[i] += (param->shift_type == 1 ? ex_lvol[i] * param->vol_shift
                                          : param->vol_shift);

    if (power > 0.5) {
      ex_lprice[i] = srt_f_optblksch(ex_lfwd[i], ex_lstrike[i], ex_lvol[i],
                                     ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
    } else {
      ex_lprice[i] = srt_f_optblknrm(ex_lfwd[i], ex_lstrike[i], ex_lvol[i],
                                     ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
    }
  }

  /*	2.)	Calibrate zeta */

  err = lgmprcapgivenlambda_2(ncpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time,
                              ex_cpn, ex_endcpn, ex_lstrike, ex_lprice, ex_zeta,
                              lambda, one2F, alpha, gamma, rho, 0);
  if (err) {
    goto FREE_RETURN;
  }

  /*	3.)	Transform into sigma */

  *num_sig = nex;
  *sig_time = (double *)calloc(nex, sizeof(double));
  *sig = (double *)calloc(nex, sizeof(double));

  if (!sig_time || !sig) {
    err = "Allocation error (3) in cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  (*sig_time)[0] = ex_time[0];
  (*sig)[0] =
      sqrt(ex_zeta[0] * 2 * (lambda) / (exp(2 * (lambda)*ex_time[0]) - 1.0));

  for (i = 1; i < nex; i++) {
    (*sig_time)[i] = ex_time[i];
    if (ex_zeta[i] > ex_zeta[i - 1]) {
      (*sig)[i] = sqrt(
          (ex_zeta[i] - ex_zeta[i - 1]) * 2 * (lambda) /
          (exp(2 * (lambda)*ex_time[i]) - exp(2 * (lambda)*ex_time[i - 1])));
    } else {
      smessage("Diagonal calibration failed at exercise year %.2f - "
               "Calibration stopped",
               ex_time[i]);
      for (j = i; j < nex; j++) {
        (*sig)[j] = (*sig)[i - 1];
      }
      i = nex;
    }
  }

  /*	4.)	Save instrument data if required */
  if (inst_data) {
    inst_data->num_inst = nex;
    inst_data->num_insts = nex;
    inst_data->start_dates = (long *)calloc(nex, sizeof(long));
    inst_data->end_dates = (long *)calloc(nex, sizeof(long));
    inst_data->start_datess = (long *)calloc(nex, sizeof(long));
    inst_data->end_datess = (long *)calloc(nex, sizeof(long));
    inst_data->short_strikes = (double *)calloc(nex, sizeof(double));
    inst_data->short_weights = (double *)calloc(nex, sizeof(double));
    inst_data->long_strikes = (double *)calloc(nex, sizeof(double));

    if (!inst_data->start_dates || !inst_data->end_dates ||
        !inst_data->start_datess || !inst_data->end_datess ||
        !inst_data->short_strikes || !inst_data->short_weights ||
        !inst_data->long_strikes) {
      err = "Allocation error (4) in cpd_calib_diagonal";
      goto FREE_RETURN;
    }

    for (i = 0; i < nex; i++) {
      inst_data->start_dates[i] = cpn_date[ex_cpn[i]];
      inst_data->end_dates[i] = cpn_date[ex_endcpn[i]];
      inst_data->start_datess[i] = cpn_date[ex_cpn[i]];
      inst_data->end_datess[i] = cpn_date[ex_endcpn[i]];
      inst_data->short_strikes[i] = 0.0;
      inst_data->short_weights[i] = 0.0;
      inst_data->long_strikes[i] = ex_lstrike[i];
    }
  }

FREE_RETURN:

  if (err) {
    if (*sig_time)
      free(*sig_time);
    *sig_time = NULL;

    if (*sig)
      free(*sig);
    *sig = NULL;

    /*
    if (inst_data)
    {
            cpd_free_calib_inst_data (inst_data);
    }
    */
  }

  return err;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//	All three routines below        , initialisation        , copying and
//freeing modified	// 	to allow the handling of long and short market
// prices. PMc 17Nov03		//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void cpd_init_calib_inst_data(CPD_CALIB_INST_DATA inst_data) {
  if (inst_data) {
    inst_data->num_inst = 0;
    inst_data->exer_dates_long = NULL;
    inst_data->start_dates = NULL;
    inst_data->end_dates = NULL;
    inst_data->long_strikes = NULL;
    inst_data->market_prices_long = NULL;

    inst_data->num_insts = 0;
    inst_data->exer_dates_short = NULL;
    inst_data->start_datess = NULL;
    inst_data->end_datess = NULL;
    inst_data->short_strikes = NULL;
    inst_data->short_weights = NULL;
    inst_data->market_prices_short = NULL;

    inst_data->num_inst_smile = 0;
    inst_data->start_dates_smile = NULL;
    inst_data->end_dates_smile = NULL;
    inst_data->alpha = NULL;
    inst_data->rho = NULL;
  }
}

/*	Copy calibration instrument data */
void cpd_copy_calib_inst_data(CPD_CALIB_INST_DATA dest,
                              CPD_CALIB_INST_DATA src) {
  if (dest && src) {
    if (src->num_inst > 0) {
      dest->num_inst = src->num_inst;

      if (src->exer_dates_long) {
        dest->exer_dates_long = (long *)calloc(dest->num_inst, sizeof(long));
        memcpy(dest->exer_dates_long, src->exer_dates_long,
               dest->num_inst * sizeof(long));
      } else {
        dest->exer_dates_long = NULL;
      }

      if (src->start_dates) {
        dest->start_dates = (long *)calloc(dest->num_inst, sizeof(long));
        memcpy(dest->start_dates, src->start_dates,
               dest->num_inst * sizeof(long));
      } else {
        dest->start_dates = NULL;
      }

      if (src->end_dates) {
        dest->end_dates = (long *)calloc(dest->num_inst, sizeof(long));
        memcpy(dest->end_dates, src->end_dates, dest->num_inst * sizeof(long));
      } else {
        dest->end_dates = NULL;
      }

      if (src->long_strikes) {
        dest->long_strikes = (double *)calloc(dest->num_inst, sizeof(double));
        memcpy(dest->long_strikes, src->long_strikes,
               dest->num_inst * sizeof(double));
      } else {
        dest->long_strikes = NULL;
      }

      if (src->market_prices_long) {
        dest->market_prices_long =
            (double *)calloc(dest->num_inst, sizeof(double));
        memcpy(dest->market_prices_long, src->market_prices_long,
               dest->num_inst * sizeof(double));
      } else {
        dest->market_prices_long = NULL;
      }
    }

    if (src->num_insts > 0) {
      dest->num_insts = src->num_insts;

      if (src->exer_dates_long) {
        dest->exer_dates_short = (long *)calloc(dest->num_insts, sizeof(long));
        memcpy(dest->exer_dates_short, src->exer_dates_short,
               dest->num_insts * sizeof(long));
      } else {
        dest->exer_dates_short = NULL;
      }

      if (src->start_datess) {
        dest->start_datess = (long *)calloc(dest->num_insts, sizeof(long));
        memcpy(dest->start_datess, src->start_datess,
               dest->num_insts * sizeof(long));
      } else {
        dest->start_datess = NULL;
      }

      if (src->end_datess) {
        dest->end_datess = (long *)calloc(dest->num_insts, sizeof(long));
        memcpy(dest->end_datess, src->end_datess,
               dest->num_insts * sizeof(long));
      } else {
        dest->end_datess = NULL;
      }

      if (src->short_strikes) {
        dest->short_strikes = (double *)calloc(dest->num_insts, sizeof(double));
        memcpy(dest->short_strikes, src->short_strikes,
               dest->num_insts * sizeof(double));
      } else {
        dest->short_strikes = NULL;
      }

      if (src->short_weights) {
        dest->short_weights = (double *)calloc(dest->num_insts, sizeof(double));
        memcpy(dest->short_weights, src->short_weights,
               dest->num_insts * sizeof(double));
      } else {
        dest->short_weights = NULL;
      }

      if (src->market_prices_short) {
        dest->market_prices_short =
            (double *)calloc(dest->num_insts, sizeof(double));
        memcpy(dest->market_prices_short, src->market_prices_short,
               dest->num_insts * sizeof(double));
      } else {
        dest->market_prices_short = NULL;
      }
    }

    if (src->num_inst_smile > 0) {
      dest->num_inst_smile = src->num_inst_smile;

      if (src->start_dates_smile) {
        dest->start_dates_smile =
            (long *)calloc(dest->num_inst_smile, sizeof(long));
        memcpy(dest->start_dates_smile, src->start_dates_smile,
               dest->num_inst_smile * sizeof(long));
      } else {
        dest->start_dates_smile = NULL;
      }

      if (src->end_dates_smile) {
        dest->end_dates_smile =
            (long *)calloc(dest->num_inst_smile, sizeof(long));
        memcpy(dest->end_dates_smile, src->end_dates_smile,
               dest->num_inst_smile * sizeof(long));
      } else {
        dest->end_dates_smile = NULL;
      }

      if (src->alpha) {
        dest->alpha = (double *)calloc(dest->num_inst_smile, sizeof(double));
        memcpy(dest->alpha, src->alpha, dest->num_inst_smile * sizeof(double));
      } else {
        dest->alpha = NULL;
      }

      if (src->rho) {
        dest->rho = (double *)calloc(dest->num_inst_smile, sizeof(double));
        memcpy(dest->rho, src->rho, dest->num_inst_smile * sizeof(double));
      } else {
        dest->rho = NULL;
      }
    }
  }
}

/*	Free calibration instrument data */
void cpd_free_calib_inst_data(CPD_CALIB_INST_DATA inst_data) {
  if (inst_data) {
    if (inst_data->num_inst > 0) {
      inst_data->num_inst = 0;

      if (inst_data->exer_dates_long) {
        free(inst_data->exer_dates_long);
        inst_data->exer_dates_long = NULL;
      }

      if (inst_data->start_dates) {
        free(inst_data->start_dates);
        inst_data->start_dates = NULL;
      }

      if (inst_data->end_dates) {
        free(inst_data->end_dates);
        inst_data->end_dates = NULL;
      }

      if (inst_data->long_strikes) {
        free(inst_data->long_strikes);
        inst_data->long_strikes = NULL;
      }

      if (inst_data->market_prices_long) {
        free(inst_data->market_prices_long);
        inst_data->market_prices_long = NULL;
      }
    }

    if (inst_data->num_insts > 0) {
      inst_data->num_insts = 0;

      if (inst_data->exer_dates_short) {
        free(inst_data->exer_dates_short);
        inst_data->exer_dates_short = NULL;
      }

      if (inst_data->start_datess) {
        free(inst_data->start_datess);
        inst_data->start_datess = NULL;
      }

      if (inst_data->end_datess) {
        free(inst_data->end_datess);
        inst_data->end_datess = NULL;
      }

      if (inst_data->short_strikes) {
        free(inst_data->short_strikes);
        inst_data->short_strikes = NULL;
      }

      if (inst_data->short_weights) {
        free(inst_data->short_weights);
        inst_data->short_weights = NULL;
      }

      if (inst_data->market_prices_short) {
        free(inst_data->market_prices_short);
        inst_data->market_prices_short = NULL;
      }
    }

    if (inst_data->num_inst_smile > 0) {
      inst_data->num_inst_smile = 0;

      if (inst_data->start_dates_smile) {
        free(inst_data->start_dates_smile);
        inst_data->start_dates_smile = NULL;
      }

      if (inst_data->end_dates_smile) {
        free(inst_data->end_dates_smile);
        inst_data->end_dates_smile = NULL;
      }

      if (inst_data->alpha) {
        free(inst_data->alpha);
        inst_data->alpha = NULL;
      }

      if (inst_data->rho) {
        free(inst_data->rho);
        inst_data->rho = NULL;
      }
    }
  }
}

/*	Calibrate a 3f model: main function */
Err cpd_calib_all(
    /*	Today */
    long today,
    /*	Get Cash Vol function */
    Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Domestic market */
    char *dom_yc_name,        /*	Name of the yield curve */
    char *dom_vol_curve_name, /*	Name of the market vol curve */
    char *dom_ref_rate_name,  /*	Name of the reference rate */
    char *dom_instr_freq,     /*	Frequency and basis of instruments */
    char *dom_instr_basis,
    /*	Domestic Model */
    double dom_lambda, /*	Lambda: may NOT be changed in the process */
    /*	Foreign market */
    char *for_yc_name,        /*	Name of the yield curve */
    char *for_vol_curve_name, /*	Name of the market vol curve */
    char *for_ref_rate_name,  /*	Name of the reference rate */
    char *for_instr_freq,     /*	Frequency and basis of instruments */
    char *for_instr_basis,
    /*	Domestic Model */
    double for_lambda, /*	Lambda: may NOT be changed in the process */
    /*	Fx market */
    long *fx_mkt_vol_date, /*	Option maturity dates */
    double *fx_mkt_vol,    /*	Option BS vol */
    int num_fx_mkt_vol,    /*	Number of options */
    /*	Fx model */
    double *corr_times, double *correl_dom_for, double *correl_dom_fx,
    double *correl_for_fx, long corr_n_times,
    /*	Structure */
    /*	If ex_date is NULL        ,
    exercise dates will be generated 2bd before start */
    int num_ex_dates, /*	Exercise dates        ,
                                                all supposed to be on or
                   after today */
    long *ex_date,    /*	Supposed to be sorted */
    int *cal_date,    /*	1: use ex_date as calibration date        , 0: don't */
    char **dom_end_tenor, /*	Tenors of the underlying instruments or "DIAG"
                           */
    char **for_end_tenor, long end_date, /*	End date for diagonal */
    double *dom_strike,                  /*	Domestic strikes 0: ATM */
    double *for_strike,                  /*	Foreign strikes 0: ATM */
    /*	Output */
    int *num_sig, /*	Answer */
    double **sig_time, double **dom_sig, double **for_sig, int *num_fx_vol,
    double **fx_vol_time, double **fx_vol,
    /*	Parameters */
    CPD_DIAG_CALIB_PARAM param,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA
        dom_inst_data, /*	NULL = don't save calibration instrument data */
    CPD_CALIB_INST_DATA
        for_inst_data) /*	NULL = don't save calibration instrument data */
{
  Err err = NULL;
  int tmp_num_sig1, tmp_num_sig2;
  double *tmp_sig_time1 = NULL, *tmp_sig_time2 = NULL;
  double *tmp_sig1 = NULL, *tmp_sig2 = NULL;
  int i, n;
  double *fx_mkt_vol_time = NULL;
  long tmp_ex_date;
  int free_dom_inst_data = 0, free_for_inst_data = 0;

  /*	Initialisation */
  *sig_time = NULL;
  *dom_sig = NULL;
  *for_sig = NULL;
  *fx_vol_time = NULL;
  *fx_vol = NULL;

  /*	1)	Calibrate domestic model */
  err = cpd_calib_diagonal_2(
      dom_yc_name, dom_vol_curve_name, dom_ref_rate_name, get_cash_vol,
      dom_instr_freq, dom_instr_basis, num_ex_dates, ex_date, cal_date,
      dom_end_tenor, end_date, dom_strike, dom_lambda, 1, 0.0, 0.0, 0.0,
      &tmp_num_sig1, &tmp_sig_time1, &tmp_sig1, param, dom_inst_data);
  if (err)
    goto FREE_RETURN;
  if (dom_inst_data)
    free_dom_inst_data = 1;

  /*	2)	Calibrate foreign model */
  err = cpd_calib_diagonal_2(
      for_yc_name, for_vol_curve_name, for_ref_rate_name, get_cash_vol,
      for_instr_freq, for_instr_basis, num_ex_dates, ex_date, cal_date,
      for_end_tenor, end_date, for_strike, for_lambda, 1, 0.0, 0.0, 0.0,
      &tmp_num_sig2, &tmp_sig_time2, &tmp_sig2, param, for_inst_data);
  if (err)
    goto FREE_RETURN;
  if (for_inst_data)
    free_for_inst_data = 1;

  /*	3)	Merge term structures */
  err = merge_rates_ts(tmp_sig_time1, tmp_sig1, tmp_num_sig1, tmp_sig_time2,
                       tmp_sig2, tmp_num_sig2, sig_time, dom_sig, for_sig,
                       num_sig);
  if (err)
    goto FREE_RETURN;
  if (tmp_sig_time1) {
    free(tmp_sig_time1);
    tmp_sig_time1 = NULL;
  }
  if (tmp_sig1) {
    free(tmp_sig1);
    tmp_sig1 = NULL;
  }
  if (tmp_sig_time2) {
    free(tmp_sig_time2);
    tmp_sig_time2 = NULL;
  }
  if (tmp_sig2) {
    free(tmp_sig2);
    tmp_sig2 = NULL;
  }

  /*	4)	Calibrate fx */

  /*	Cut vol curve to end of structure */
  n = num_fx_mkt_vol - 1;
  while (n > 0 && fx_mkt_vol_date[n] > end_date) {
    n--;
  }
  if (n < num_fx_mkt_vol - 1 && fx_mkt_vol_date[n] < end_date) {
    n++;
  }
  num_fx_mkt_vol = n + 1;

  /*	Calibrate */
  fx_mkt_vol_time = (double *)calloc(num_fx_mkt_vol, sizeof(double));
  *fx_vol_time = (double *)calloc(num_fx_mkt_vol, sizeof(double));
  if (!fx_mkt_vol_time || !*fx_vol_time) {
    err = "Allocation error in cpd_calib_all";
    goto FREE_RETURN;
  }

  for (i = 0; i < num_fx_mkt_vol; i++) {
    fx_mkt_vol_time[i] = (fx_mkt_vol_date[i] - today) * YEARS_IN_DAY;
    fx_mkt_vol[i] += param->fx_vol_shift;
    tmp_ex_date =
        add_unit(fx_mkt_vol_date[i], -2, SRT_BDAY, MODIFIED_SUCCEEDING);
    (*fx_vol_time)[i] = (tmp_ex_date - today) * YEARS_IN_DAY;
  }

  err = Fx3DtsCalibration_corr(
      *fx_vol_time, fx_mkt_vol_time, fx_mkt_vol, num_fx_mkt_vol, *sig_time,
      *num_sig, *dom_sig, dom_lambda, *for_sig, for_lambda, corr_times,
      correl_dom_for, correl_dom_fx, correl_for_fx, corr_n_times, fx_vol);
  *num_fx_vol = num_fx_mkt_vol;
  if (err) {
    goto FREE_RETURN;
  }
  if (fx_mkt_vol_time) {
    free(fx_mkt_vol_time);
    fx_mkt_vol_time = NULL;
  }

FREE_RETURN:

  if (tmp_sig_time1)
    free(tmp_sig_time1);
  if (tmp_sig1)
    free(tmp_sig1);
  if (tmp_sig_time2)
    free(tmp_sig_time2);
  if (tmp_sig2)
    free(tmp_sig2);

  if (fx_mkt_vol_time)
    free(fx_mkt_vol_time);

  if (err) {
    if (dom_inst_data && free_dom_inst_data) {
      cpd_free_calib_inst_data(dom_inst_data);
      dom_inst_data = NULL;
    }
    if (for_inst_data && free_dom_inst_data) {
      cpd_free_calib_inst_data(for_inst_data);
      for_inst_data = NULL;
    }

    if (*sig_time) {
      free(*sig_time);
      *sig_time = NULL;
    }
    if (*dom_sig) {
      free(*dom_sig);
      *dom_sig = NULL;
    }
    if (*for_sig) {
      free(*for_sig);
      *for_sig = NULL;
    }
    if (*fx_vol_time) {
      free(*fx_vol_time);
      *fx_vol_time = NULL;
    }
    if (*fx_vol) {
      free(*fx_vol);
      *fx_vol = NULL;
    }
  }

  return err;
}

#include "srtaccess.h"

/*	Make underlying out of the results of the previous function */
Err cpd_calib_all_makeund(long today, char *dom_ccy, char *dom_und_name,
                          char *dom_yc_name, double *dom_sig, double dom_lambda,
                          char *for_ccy, char *for_und_name, char *for_yc_name,
                          double *for_sig, double for_lambda,
                          double *corr_times, double *correl_dom_for,
                          double *correl_dom_fx, double *correl_for_fx,
                          long corr_n_times, int num_sig, double *sig_time,
                          int num_fx_vol, double spot_fx, double *fx_vol_time,
                          double *fx_vol,
                          /*	Output only */
                          char fx_und_name[]) {
  TermStruct *ts1 = NULL, *ts2 = NULL, *ts3 = NULL;
  double sig_dates[256];
  double *vol_crv[2], **tau_crv;
  double tau, *tmp;
  SrtUndListPtr und_list;
  double **corr = NULL, *corrdate = NULL;
  char ***corrundtab = NULL;
  int i;
  int dom_init = 0, for_init = 0, corr_init = 0, fx_init = 0;

  Err err = NULL;

  /*	Get the underlying list and check if not empty */
  und_list = get_underlying_list();
  if (und_list == NULL) {
    err = "No Underlying list defined: call SrtInit";
    goto FREE_RETURN;
  }

  /*	1)	Dom und */

  /*	Term Structure */
  for (i = 0; i < num_sig; i++) {
    sig_dates[i] = today + (long)(1.0e-08 + DAYS_IN_YEAR * sig_time[i]);
  }
  vol_crv[0] = &(sig_dates[0]);
  vol_crv[1] = dom_sig;
  tau = 1.0 / dom_lambda;
  tmp = &tau;
  tau_crv = &tmp;

  err = srt_f_init_IRM_OneFac_TermStruct(&ts1, today, vol_crv, 2, num_sig,
                                         tau_crv, 1, 1, LGM, 0.0, 0.0, 0.0, 0.0,
                                         0.0, 0, 0, NULL);
  if (err)
    goto FREE_RETURN;

  /*	Und */

  /*	Makes the Underlying Name UpperCase */
  strupper(dom_und_name);

  /*	Puts the Underlying in the Market List (updating its ticker) */
  err = srt_f_addundtolist(und_list, dom_und_name, "IR_UND", dom_ccy, "LGM",
                           dom_yc_name, NULL, NULL, ts1, 0.0);
  if (err)
    goto FREE_RETURN;
  dom_init = 1;

  /*	2)	For und */

  /*	Term Structure */
  vol_crv[1] = for_sig;
  tau = 1.0 / for_lambda;
  tmp = &tau;
  tau_crv = &tmp;

  err = srt_f_init_IRM_OneFac_TermStruct(&ts2, today, vol_crv, 2, num_sig,
                                         tau_crv, 1, 1, LGM, 0.0, 0.0, 0.0, 0.0,
                                         0.0, 0, 0, NULL);
  if (err)
    goto FREE_RETURN;

  /*	Und */

  /*	Makes the Underlying Name UpperCase */
  strupper(for_und_name);

  /*	Puts the Underlying in the Market List (updating its ticker) */
  err = srt_f_addundtolist(und_list, for_und_name, "IR_UND", for_ccy, "LGM",
                           for_yc_name, NULL, NULL, ts2, 0.0);
  if (err)
    goto FREE_RETURN;
  for_init = 1;

  /*	3)	Correlation matrix */

  /*	Define the fx underlying name */
  strcpy(fx_und_name, for_ccy);
  strcat(fx_und_name, "/");
  strcat(fx_und_name, dom_ccy);
  strupper(fx_und_name);

  corrundtab = smatrix_size(0, 2, 0, 1, 256);
  corr = dmatrix(0, 2, 0, corr_n_times - 1);
  corrdate = dvector(0, corr_n_times - 1);

  if (!corrundtab || !corr || !corrdate) {
    err = "Memory allocation failure in cpd_calib_all_makeund";
    goto FREE_RETURN;
  }

  strcpy(corrundtab[0][0], dom_und_name);
  strcpy(corrundtab[0][1], for_und_name);
  strcpy(corrundtab[1][0], dom_und_name);
  strcpy(corrundtab[1][1], fx_und_name);
  strcpy(corrundtab[2][0], for_und_name);
  strcpy(corrundtab[2][1], fx_und_name);

  for (i = 0; i < corr_n_times; i++) {
    corr[0][i] = correl_dom_for[i];
    corr[1][i] = correl_dom_fx[i];
    corr[2][i] = correl_for_fx[i];
    corrdate[i] = today + DAYS_IN_YEAR * corr_times[i] + 1.0E-08;
  }

  err = SrtInitCorrelationMatrix(corr_n_times, 3, corr, corrdate, corrundtab);

  corr_init = 1;

  if (err) {
    goto FREE_RETURN;
  }

  /*	4)	Fx und */

  /*	Init term struct */
  for (i = 0; i < num_fx_vol; i++) {
    sig_dates[i] = today + (long)(1.0e-08 + DAYS_IN_YEAR * fx_vol_time[i]);
  }
  vol_crv[1] = fx_vol;

  err = srt_f_init_FX_TermStruct(today, vol_crv, 2, num_fx_vol, FX_STOCH_RATES,
                                 0.0, 0.0, fx_und_name, dom_und_name,
                                 for_und_name, &ts3);
  if (err)
    goto FREE_RETURN;

  /*	Und */
  err = srt_f_addundtolist(und_list, fx_und_name, "FX_UND", dom_ccy,
                           "FX_STOCH_RATES", dom_und_name, for_und_name, NULL,
                           ts3, spot_fx);
  fx_init = 1;

FREE_RETURN:

  if (corrundtab) {
    free_smatrix_size(corrundtab, 0, 2, 0, 1, 256);
  }
  if (corr) {
    free_dmatrix(corr, 0, 2, 0, corr_n_times - 1);
  }
  if (corrdate) {
    free_dvector(corrdate, 0, corr_n_times - 1);
  }

  if (err) {
    if (dom_init) {
      srt_f_destroy_und(dom_und_name);
    }
    if (for_init) {
      srt_f_destroy_und(for_und_name);
    }
    if (fx_init) {
      srt_f_destroy_und(fx_und_name);
    }
    if (corr_init) {
      destroy_correlation_list();
    }
  }

  return err;
}

/*	Calibrate lgm: main function */
/*	New version: calibrates not necessarily to digonal
                with lambda calibration */
Err cpd_calib_diagonal_3(
    /*	Market */
    char *yc_name,        /*	Name of the yield curve */
    char *vol_curve_name, /*	Name of the market vol curve */
    char *ref_rate_name,  /*	Name of the reference rate */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    char *instr_freq, /*	Frequency and basis of instruments */
    char *instr_basis,

    /*	If ex_date is NULL        ,
    exercise dates will be generated 2bd before start */
    /*	Structure */
    int num_ex_dates, /*	Exercise dates        ,
                                                all supposed to be on or
                   after today */
    long *ex_date_,   /*	Supposed to be sorted */
    int *cal_date,    /*	1: use ex_date as calibration date        , 0: don't */
    char **end_tenorl_, /*	Tenors of the underlying instruments
                                                  or "DIAG" */
    char **end_tenors_, /*	Tenors of the underlying instruments
                                                  or "DIAG" */
    long end_date,      /*	End date for diagonal */
    double *strikel_,   /*	Strikes
                                          0: ATM */
    double *strikes_,   /*	Strikes
                                          0: ATM */
    /*	Model */
    int fix_lambda,
    int nlam, /*	Lambda TS: may NOT be changed in the process */
    double lam_time[], double lam[], int one2F, /*	Number of factors */
    double alpha, /*	Alpha        , Gamma        , Rho (2F only) */
    double gamma, double rho,
    /*	Output */
    int *num_sig, /*	Answer */
    double **sig_time, double **sig,
    /*	Parameters */
    CPD_DIAG_CALIB_PARAM param, DIAG_CALIB_LM_PARAMS lm_params,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA
        inst_data) /*	NULL = don't save calibration instrument data */
{
  int i, j, k, l, nex, ncpn;
  SrtCompounding ifreq;
  SrtBasisCode ibasis;
  double ex_time[MAX_CPN], ex_lfwd[MAX_CPN], ex_llvl[MAX_CPN],
      ex_lstrike[MAX_CPN], ex_lvol[MAX_CPN], ex_lprice[MAX_CPN],
      ex_sfwd[MAX_CPN], ex_slvl[MAX_CPN], ex_sstrike[MAX_CPN], ex_svol[MAX_CPN],
      ex_sprice[MAX_CPN], ex_svega[MAX_CPN], ex_zeta[MAX_CPN];
  int ex_cpn[MAX_CPN], ex_lendcpn[MAX_CPN], ex_sendcpn[MAX_CPN];
  long cpn_date[MAX_CPN];
  double cpn_time[MAX_CPN], cpn_cvg[MAX_CPN], cpn_df[MAX_CPN];
  long tmplngl1[MAX_CPN], tmplngl2[MAX_CPN];
  long tmplngs1[MAX_CPN], tmplngs2[MAX_CPN];
  long *theo_end_datesl, *act_end_datesl, *theo_end_datess, *act_end_datess;
  long theo_date, act_date, temp_date, temp_date2;
  long today;
  double lvl, dfi, dff;
  double power, exp_fact;
  double swp_rte, spr;
  double std;
  SrtCurvePtr yc_ptr;
  Err err = NULL;
  long ex_date__[MAX_INST], *ex_date;
  double strikel__[MAX_INST], *strikel;
  char end_tenorl__[MAX_INST * 256];
  String end_tenorl___[MAX_INST], *end_tenorl;
  double strikes__[MAX_INST], *strikes;
  char end_tenors__[MAX_INST * 256];
  String end_tenors___[MAX_INST], *end_tenors;

  /*	Copy data so as not to change the original */
  ex_date = &(ex_date__[0]);

  strikel = &(strikel__[0]);
  end_tenorl = &(end_tenorl___[0]);

  end_tenors = &(end_tenors___[0]);
  strikes = &(strikes__[0]);

  for (i = 0; i < num_ex_dates; i++) {
    end_tenorl[i] = &(end_tenorl__[i * 256]);
    end_tenors[i] = &(end_tenors__[i * 256]);
  }
  memcpy(ex_date, ex_date_, num_ex_dates * sizeof(long));
  memcpy(strikel, strikel_, num_ex_dates * sizeof(double));
  memcpy(strikes, strikes_, num_ex_dates * sizeof(double));

  for (i = 0; i < num_ex_dates; i++) {
    strcpy(end_tenorl[i], end_tenorl_[i]);
    strcpy(end_tenors[i], end_tenors_[i]);
  }

  theo_end_datesl = &(tmplngl1[0]);
  act_end_datesl = &(tmplngl2[0]);

  theo_end_datess = &(tmplngs1[0]);
  act_end_datess = &(tmplngs2[0]);

  *sig_time = NULL;
  *sig = NULL;

  yc_ptr = lookup_curve(yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }
  today = get_today_from_curve(yc_ptr);

  if (one2F == 2) {
    HermiteStandard(x, w, NUM_HERMITE);
  }

  /*	1.)	Setup the bond schedule and its coupons */

  /*	Coupons */

  err = interp_compounding(instr_freq, &ifreq);
  if (err) {
    goto FREE_RETURN;
  }

  err = interp_basis(instr_basis, &ibasis);
  if (err) {
    goto FREE_RETURN;
  }

  /*	Find the end date as the longest total maturity */
  theo_date = end_date;
  act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
  for (i = 0; i < num_ex_dates; i++) {
    err = get_end_date(ex_date[i], end_date, end_tenorl[i], 0,
                       &(theo_end_datesl[i]));
    if (err) {
      goto FREE_RETURN;
    }
    act_end_datesl[i] =
        bus_date_method(theo_end_datesl[i], MODIFIED_SUCCEEDING);

    err = get_end_date(ex_date[i], end_date, end_tenors[i], 0,
                       &(theo_end_datess[i]));
    if (err) {
      goto FREE_RETURN;
    }
    act_end_datess[i] =
        bus_date_method(theo_end_datess[i], MODIFIED_SUCCEEDING);
  }
  for (i = 0; i < num_ex_dates; i++) {
    if (theo_end_datesl[i] > theo_date || act_end_datesl[i] > act_date) {
      theo_date = theo_end_datesl[i];
      act_date = act_end_datesl[i];
    }

    if (theo_end_datess[i] > theo_date || act_end_datess[i] > act_date) {
      theo_date = theo_end_datess[i];
      act_date = act_end_datess[i];
    }
  }
  ncpn = 1;
  temp_date = theo_date;
  temp_date2 = act_date;

  while (act_date > today) {
    theo_date =
        add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
    act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    ncpn++;
  }
  ncpn--;

  if (ncpn < 2) {
    err = "Not enough coupons in cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  theo_date = temp_date;
  act_date = temp_date2;
  i = ncpn - 1;

  while (i >= 0) {
    cpn_time[i] = (act_date - today) * YEARS_IN_DAY;
    cpn_date[i] = act_date;

    theo_date =
        add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

    temp_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    cpn_cvg[i] = coverage(temp_date, act_date, ibasis);
    cpn_df[i] = swp_f_df(today, act_date, yc_name);
    act_date = temp_date;

    i--;
  }
  cpn_cvg[0] = 0.0;

  /*	Exercise */

  /*	Remove non-calibration dates */
  for (i = num_ex_dates - 1; i >= 0; i--) {
    if (cal_date[i] == 0) {
      for (k = i - 1; k >= 0; k--) {
        ex_date[k + 1] = ex_date[k];
        strcpy(end_tenorl[k + 1], end_tenorl[k]);
        theo_end_datesl[k + 1] = theo_end_datesl[k];
        act_end_datesl[k + 1] = act_end_datesl[k];
        strikel[k + 1] = strikel[k];

        strcpy(end_tenors[k + 1], end_tenors[k]);
        theo_end_datess[k + 1] = theo_end_datess[k];
        act_end_datess[k + 1] = act_end_datess[k];
        strikes[k + 1] = strikes[k];
      }

      ex_date++;
      end_tenorl++;
      theo_end_datesl++;
      act_end_datesl++;
      strikel++;
      end_tenors++;
      theo_end_datess++;
      act_end_datess++;
      strikes++;
      num_ex_dates--;
      if (num_ex_dates < 1) {
        err = "All exercise dates are past in cpd_calib_diagonal";
        goto FREE_RETURN;
      }
    }
  }

  /*	Remove redundant dates */
  j = ncpn - 1;
  l = ncpn + 1;
  for (i = num_ex_dates - 1; i >= 0; i--) {
    while (j > 0 && cpn_date[j] > ex_date[i]) {
      j--;
    }
    if (cpn_date[j] < ex_date[i]) {
      j++;
    }

    if (j >= ncpn - 1 || j == l) {
      for (k = i - 1; k >= 0; k--) {
        ex_date[k + 1] = ex_date[k];
        strcpy(end_tenorl[k + 1], end_tenorl[k]);
        theo_end_datesl[k + 1] = theo_end_datesl[k];
        act_end_datesl[k + 1] = act_end_datesl[k];
        strikel[k + 1] = strikel[k];
        strcpy(end_tenors[k + 1], end_tenors[k]);
        theo_end_datess[k + 1] = theo_end_datess[k];
        act_end_datess[k + 1] = act_end_datess[k];
        strikes[k + 1] = strikes[k];
      }

      ex_date++;
      end_tenorl++;
      theo_end_datesl++;
      act_end_datesl++;
      strikel++;
      end_tenors++;
      theo_end_datess++;
      act_end_datess++;
      strikes++;
      num_ex_dates--;
      if (num_ex_dates < 1) {
        err = "All exercise dates are past in cpd_calib_diagonal";
        goto FREE_RETURN;
      }
    } else {
      l = j;
    }
  }

  /*	Remove close dates */
  j = num_ex_dates - 1;
  for (i = num_ex_dates - 2; i >= 0; i--) {
    if ((ex_date[j] - ex_date[i]) * YEARS_IN_DAY <
        param->min_time - ONE_MONTH) {
      for (k = i - 1; k >= 0; k--) {
        ex_date[k + 1] = ex_date[k];
        strcpy(end_tenorl[k + 1], end_tenorl[k]);
        theo_end_datesl[k + 1] = theo_end_datesl[k];
        act_end_datesl[k + 1] = act_end_datesl[k];
        strikel[k + 1] = strikel[k];
        strcpy(end_tenors[k + 1], end_tenors[k]);
        theo_end_datess[k + 1] = theo_end_datess[k];
        act_end_datess[k + 1] = act_end_datess[k];
        strikes[k + 1] = strikes[k];
      }

      ex_date++;
      end_tenorl++;
      theo_end_datesl++;
      act_end_datesl++;
      strikel++;
      end_tenors++;
      theo_end_datess++;
      act_end_datess++;
      strikes++;
      num_ex_dates--;
      j--;
      if (num_ex_dates < 1) {
        err = "All exercise dates are past in cpd_calib_diagonal";
        goto FREE_RETURN;
      }
    } else {
      j = i;
    }
  }

  /*	Remove last? */
  if (param->skip_last && num_ex_dates > 1) {
    num_ex_dates--;
  }

  /*	Remove past dates */
  while (ex_date[0] <= today) {
    ex_date++;
    end_tenorl++;
    theo_end_datesl++;
    act_end_datesl++;
    strikel++;
    end_tenors++;
    theo_end_datess++;
    act_end_datess++;
    strikes++;
    num_ex_dates--;
    if (num_ex_dates == 0) {
      "All exercise dates are past in cpd_calib_diagonal";
      goto FREE_RETURN;
    }
  }

  nex = num_ex_dates;
  j = 0;
  for (i = 0; i < nex; i++) {
    while (cpn_date[j] < ex_date[i]) {
      j++;
    }

    ex_cpn[i] = j;
    ex_time[i] = (ex_date[i] - today) * YEARS_IN_DAY;

    /* Long Underlying */
    k = j;
    while (cpn_date[k] < act_end_datesl[i]) {
      k++;
    }
    if (k > 0 &&
        cpn_date[k] - act_end_datesl[i] > act_end_datesl[i] - cpn_date[k - 1]) {
      k--;
    }

    if (k <= j) {
      k = j + 1;
    }
    ex_lendcpn[i] = k;

    if (j >= ncpn || k >= ncpn) {
      err = "Coupon date bug in long underlyings in cpd_calib_diagonal";
      goto FREE_RETURN;
    }

    /* Short Underlying */
    k = j;
    while (cpn_date[k] < act_end_datess[i]) {
      k++;
    }
    if (k > 0 &&
        cpn_date[k] - act_end_datess[i] > act_end_datess[i] - cpn_date[k - 1]) {
      k--;
    }

    if (k <= j) {
      k = j + 1;
    }
    ex_sendcpn[i] = k;

    if (j >= ncpn || k >= ncpn) {
      err = "Coupon date bug in short underlying in cpd_calib_diagonal";
      goto FREE_RETURN;
    }
  }

  /*	Underlyings */

  /*	Long */

  for (i = 0; i < nex; i++) {
    j = ex_cpn[i];
    l = ex_lendcpn[i];

    lvl = 0.0;
    for (k = j + 1; k <= l; k++) {
      lvl += cpn_cvg[k] * cpn_df[k];
    }
    dfi = swp_f_df(today, cpn_date[j], yc_name);
    dff = swp_f_df(today, cpn_date[l], yc_name);

    ex_llvl[i] = lvl;
    ex_lfwd[i] = (dfi - dff) / lvl;

    /*	ATM std */
    err = get_cash_vol(
        vol_curve_name, add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
        theo_end_datesl[i], ex_lfwd[i], 0, ref_rate_name, &std, &power);
    if (err) {
      goto FREE_RETURN;
    }
    std += (param->shift_type == 1 ? std * param->vol_shift : param->vol_shift);
    if (power > 0.5) {
      power = srt_f_optblksch(ex_lfwd[i], ex_lfwd[i], std, ex_time[i], 1.0,
                              SRT_CALL, PREMIUM);
      err = srt_f_optimpvol(power, ex_lfwd[i], ex_lfwd[i], ex_time[i], 1.0,
                            SRT_CALL, SRT_NORMAL, &std);
    }
    std *= sqrt(ex_time[i]);

    /*	Strike */
    if (param->strike_type == 0 || strikel[i] < 1.0e-04) {
      ex_lstrike[i] = ex_lfwd[i];
    } else if (param->strike_type == 1) {
      ex_lstrike[i] = strikel[i];
    } else if (param->strike_type == 2) {
      if (err = swp_f_ForwardRate(cpn_date[j], theo_end_datesl[i], instr_freq,
                                  instr_basis, yc_name, ref_rate_name,
                                  &swp_rte)) {
        goto FREE_RETURN;
      }

      spr = swp_rte - ex_lfwd[i];

      ex_lstrike[i] = strikel[i] - spr;
    } else if (param->strike_type == 3) {
      ex_lstrike[i] = ex_lfwd[i] + strikel[i] * std;
    }

    /*	Apply max std */
    if (ex_lstrike[i] > ex_lfwd[i] + param->max_std * std) {
      ex_lstrike[i] = ex_lfwd[i] + param->max_std * std;
    } else if (ex_lstrike[i] < ex_lfwd[i] - param->max_std * std) {
      ex_lstrike[i] = ex_lfwd[i] - param->max_std * std;
    }

    /*	Make sure strikes are positive (actually more than 1bp)
                    otherwise use ATM	*/
    if (ex_lstrike[i] < 1.0e-04) {
      ex_lstrike[i] = ex_lfwd[i];
    }

    err = get_cash_vol(vol_curve_name,
                       add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                       theo_end_datesl[i], ex_lstrike[i], 0, ref_rate_name,
                       &(ex_lvol[i]), &power);
    if (err) {
      goto FREE_RETURN;
    }
    ex_lvol[i] += (param->shift_type == 1 ? ex_lvol[i] * param->vol_shift
                                          : param->vol_shift);

    if (power > 0.5) {
      ex_lprice[i] = srt_f_optblksch(ex_lfwd[i], ex_lstrike[i], ex_lvol[i],
                                     ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
    } else {
      ex_lprice[i] = srt_f_optblknrm(ex_lfwd[i], ex_lstrike[i], ex_lvol[i],
                                     ex_time[i], ex_llvl[i], SRT_PUT, PREMIUM);
    }
  }

  /*	Short */

  for (i = 0; i < nex; i++) {
    j = ex_cpn[i];
    l = ex_sendcpn[i];

    lvl = 0.0;
    for (k = j + 1; k <= l; k++) {
      lvl += cpn_cvg[k] * cpn_df[k];
    }
    dfi = swp_f_df(today, cpn_date[j], yc_name);
    dff = swp_f_df(today, cpn_date[l], yc_name);

    ex_slvl[i] = lvl;
    ex_sfwd[i] = (dfi - dff) / lvl;

    /*	ATM std */
    err = get_cash_vol(
        vol_curve_name, add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
        theo_end_datess[i], ex_sfwd[i], 0, ref_rate_name, &std, &power);
    if (err) {
      goto FREE_RETURN;
    }
    std += (param->shift_type == 1 ? std * param->vol_shift : param->vol_shift);
    if (power > 0.5) {
      power = srt_f_optblksch(ex_sfwd[i], ex_sfwd[i], std, ex_time[i], 1.0,
                              SRT_CALL, PREMIUM);
      err = srt_f_optimpvol(power, ex_sfwd[i], ex_sfwd[i], ex_time[i], 1.0,
                            SRT_CALL, SRT_NORMAL, &std);
    }
    std *= sqrt(ex_time[i]);

    /*	Strike */
    if (param->strike_type == 0 || strikes[i] < 1.0e-04) {
      ex_sstrike[i] = ex_sfwd[i];
    } else if (param->strike_type == 1) {
      ex_sstrike[i] = strikes[i];
    } else if (param->strike_type == 2) {
      if (err = swp_f_ForwardRate(cpn_date[j], theo_end_datess[i], instr_freq,
                                  instr_basis, yc_name, ref_rate_name,
                                  &swp_rte)) {
        goto FREE_RETURN;
      }

      spr = swp_rte - ex_sfwd[i];

      ex_sstrike[i] = strikes[i] - spr;
    } else if (param->strike_type == 3) {
      ex_sstrike[i] = ex_sfwd[i] + strikes[i] * std;
    }

    /*	Apply max std */
    if (ex_sstrike[i] > ex_sfwd[i] + param->max_std * std) {
      ex_sstrike[i] = ex_sfwd[i] + param->max_std * std;
    } else if (ex_sstrike[i] < ex_sfwd[i] - param->max_std * std) {
      ex_sstrike[i] = ex_sfwd[i] - param->max_std * std;
    }

    /*	Make sure strikes are positive (actually more than 1bp)
                    otherwise use ATM	*/
    if (ex_sstrike[i] < 1.0e-04) {
      ex_sstrike[i] = ex_sfwd[i];
    }

    err = get_cash_vol(vol_curve_name,
                       add_unit(ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                       theo_end_datess[i], ex_sstrike[i], 0, ref_rate_name,
                       &(ex_svol[i]), &power);
    if (err) {
      goto FREE_RETURN;
    }
    ex_svol[i] += (param->shift_type == 1 ? ex_svol[i] * param->vol_shift
                                          : param->vol_shift);

    if (power > 0.5) {
      ex_sprice[i] = srt_f_optblksch(ex_sfwd[i], ex_sstrike[i], ex_svol[i],
                                     ex_time[i], ex_slvl[i], SRT_PUT, PREMIUM);

      ex_svega[i] = srt_f_optblksch(ex_sfwd[i], ex_sstrike[i], ex_svol[i],
                                    ex_time[i], ex_slvl[i], SRT_PUT, VEGA);
    } else {
      ex_sprice[i] = srt_f_optblknrm(ex_sfwd[i], ex_sstrike[i], ex_svol[i],
                                     ex_time[i], ex_slvl[i], SRT_PUT, PREMIUM);

      ex_svega[i] = srt_f_optblknrm(ex_sfwd[i], ex_sstrike[i], ex_svol[i],
                                    ex_time[i], ex_slvl[i], SRT_PUT, VEGA);
    }
  }

  /*	2.)	Calibrate zeta */

  err = lgmcalibzetalambda_tauts(ncpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time,
                                 ex_cpn, ex_lendcpn, ex_lstrike, ex_lprice,
                                 ex_sendcpn, ex_sstrike, ex_sprice, ex_svega,
                                 ex_zeta, fix_lambda, nlam, lam_time, lam,
                                 one2F, alpha, gamma, rho, 0, lm_params);

  if (err) {
    goto FREE_RETURN;
  }

  /*	3.)	Transform into sigma */

  *num_sig = nex;
  *sig_time = (double *)calloc(nex, sizeof(double));
  *sig = (double *)calloc(nex, sizeof(double));

  if (!sig_time || !sig) {
    err = "Allocation error (3) in cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  (*sig_time)[0] = ex_time[0];
  exp_fact = static_lgmcalcexpfact_tauts(0.0, ex_time[0], nlam, lam_time, lam);
  (*sig)[0] = sqrt(ex_zeta[0] / exp_fact);

  for (i = 1; i < nex; i++) {
    (*sig_time)[i] = ex_time[i];
    if (ex_zeta[i] > ex_zeta[i - 1]) {
      exp_fact = static_lgmcalcexpfact_tauts(ex_time[i - 1], ex_time[i], nlam,
                                             lam_time, lam);
      (*sig)[i] = sqrt((ex_zeta[i] - ex_zeta[i - 1]) / exp_fact);
    } else {
      smessage("Diagonal calibration failed at exercise year %.2f - "
               "Calibration stopped",
               ex_time[i]);
      for (j = i; j < nex; j++) {
        (*sig)[j] = (*sig)[i - 1];
      }
      i = nex;
    }
  }

  /*	4.)	Save instrument data if required */
  if (inst_data) {
    inst_data->num_inst = nex;
    inst_data->num_insts = nex;
    inst_data->start_dates = (long *)calloc(nex, sizeof(long));
    inst_data->end_dates = (long *)calloc(nex, sizeof(long));
    inst_data->start_datess = (long *)calloc(nex, sizeof(long));
    inst_data->end_datess = (long *)calloc(nex, sizeof(long));
    inst_data->short_strikes = (double *)calloc(nex, sizeof(double));
    inst_data->short_weights = (double *)calloc(nex, sizeof(double));
    inst_data->long_strikes = (double *)calloc(nex, sizeof(double));

    if (!inst_data->start_dates || !inst_data->end_dates ||
        !inst_data->start_datess || !inst_data->end_datess ||
        !inst_data->short_strikes || !inst_data->short_weights ||
        !inst_data->long_strikes) {
      err = "Allocation error (4) in cpd_calib_diagonal";
      goto FREE_RETURN;
    }

    for (i = 0; i < nex; i++) {
      inst_data->start_dates[i] = cpn_date[ex_cpn[i]];
      inst_data->end_dates[i] = cpn_date[ex_lendcpn[i]];
      inst_data->start_datess[i] = cpn_date[ex_cpn[i]];
      inst_data->end_datess[i] = cpn_date[ex_lendcpn[i]];
      inst_data->short_strikes[i] = 0.0;
      inst_data->short_weights[i] = 1.0;
      inst_data->long_strikes[i] = ex_lstrike[i];
    }
  }

FREE_RETURN:

  if (err) {
    if (*sig_time)
      free(*sig_time);
    *sig_time = NULL;

    if (*sig)
      free(*sig);
    *sig = NULL;

    if (inst_data) {
      cpd_free_calib_inst_data(inst_data);
    }
  }

  return err;
}