/* ===================================================================================
   FILENAME:      srt_f_betaetaclsdfrm.c

   PURPOSE:       Compute swaption  , bond option  ,caps/floors prices in
   betaeta  , via a closed form analytical solution.
   ===================================================================================
 */

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_betaetaclsdfrm.h"

#define MAXITER 20

static double betaeta_discount_bond_option(
    double fixing_time, double period_start_time, double period_end_time,
    double period_start_df, double period_end_df, double bond_strike,
    TermStruct *ts, SrtReceiverType pay_rec, SrtMdlDim mdl_dim);

/* Evaluates cash flows for a swaption as sum(CF(t  ,x)/N(t  ,x))*density(x) */
static double swaption_cfs(double y_bar);

/* Evaluates cash flows for a cap/floor as sum(CF(t  ,x)/N(t  ,x))*density(x) */
static double capfl_cfs(double y_bar);

/* Density corresponding to power = 0<=<.5 */
static double density05(double y_bar);

/* Density corresponding to power =.5<=<1 */
static double density51(double y_bar);

/* Density corresponding to power = 1 */
static double density1(double y_bar);

/* Compute the range of integration for the int^b_a(cf's*density) */
static Err int_range(double *a, double *b, double pay_time, TermStruct *ts,
                     double (*g_x)(double));
/* -------------------------------------------------------------------------
   FUNCTION: srt_f_betaetaclsdfrm

   PURPOSE:  betaeta model --- closed form prices for swaptions  , bond options
   , caps and floors. When dealing with a swaption  , a cap or a floor  , the
   bond strike has no impact. Please note that if it is an option on a bond  ,
   the bond strike is assumed to be a clean strike. The coupons are calculated
                         with the coverage when generating the fixed leg: it is
   a clean bond generated
   ------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
        Static variables for use with Simpson's rule
   ----------------------------------------------------------------------------
 */
int STATIC_n;
double STATIC_bond_strike;
TermStruct *STATIC_ts;
double STATIC_fixing_time;
double *STATIC_coupon;
double *STATIC_pay_time;
double *STATIC_df;
SrtReceiverType STATIC_pay_rec;
SrtMdlDim STATIC_mdl_dim;
int STATIC_rp;
double STATIC_beta;
double STATIC_eta;
double STATIC_deltau;
double STATIC_x_k;
double STATIC_lambda;
double STATIC_period_start_time;
double STATIC_period_end_time;
double STATIC_period_start_df;
double STATIC_period_end_df;
int STATIC_density_off;
/* ----------------------------------------------------------------------------
        End Static variables for use with Simpson's rule
   ----------------------------------------------------------------------------
 */

/* ==================================================================

    srt_f_betaeta_coupon_bond_option()

        prices a bond option or a swaption (option on bond at par) using the
        betaeta model

        the inputs:
        n		     dates/coupons of underlying go from [0] to [n]
        bond_strike	 strike of swaption/bond option
        *ts		     term structure of volatility
        *fixing      time    time in years from today to fixing date
        *coupon      array of n coupons of underlying swap/bond  ,
                                 starting with coupon[1]
        *pay_time	 array of n times of corresponding coupons (in years  ,
                                 from today);
                                 Index of first coupon is 1.  pay_time[0]
                                 is the maturity of the option and not used.
                                 (Today is 0.0  , hence there is no
                                 variable "today")

        *df	array of n+1 dfs for n+1 times in pay_time.
                        (indexed as above)

      pay_rec        rec means call on bond  , which is equivalent
                        to put on rates; receive fixed (receiver);
                        pay means put on bond; call on rates; pay fixed

        returns price.
   ======================================================================== */

double srt_f_etabeta_coupon_bond_option(int n, double bond_strike,
                                        TermStruct *ts, double fixing_time,
                                        double *coupon, double *pay_time,
                                        double *df, SrtReceiverType pay_rec,
                                        SrtMdlDim mdl_dim) {
  double a, b;      /* integration bounds */
  double precision; /* precision for simpson */
  double price;
  int type; /* used for int_range 1 for swaption 0 for cap */
  Err err = NULL;

  /* Check if option has expired */
  if (pay_time[0] < 0.0)
    return 0.0;

  STATIC_rp = (pay_rec == SRT_RECEIVER ? 1 : -1);

  /* Assign values to static global variables for simpson integration */

  STATIC_n = n;
  STATIC_bond_strike = bond_strike;
  STATIC_ts = ts;
  STATIC_fixing_time = fixing_time;
  STATIC_coupon = coupon;
  STATIC_pay_time = pay_time;
  STATIC_df = df;
  STATIC_pay_rec = pay_rec;
  STATIC_mdl_dim = mdl_dim;

  /* Assign integration range */

  a = 0.;
  b = 0.;
  type = 1; /* swaption */

  /* Finds the critical x* that corresponds to the exercise ( intrinsic(x*) =
   * strike ) */
  err = int_range(&a, &b, pay_time[0], ts, swaption_cfs);
  /* Need to check for errors
  if (err)
  {
  }
  */

  /* Set precision for simpson  */
  precision = .001;

  /* Compute expected value by integrating cash flows against the density */
  if (a < b)
    price = sm_qsimp(swaption_cfs, a, b, precision);
  else
    price = 0;

  return price;
} /* END betaeta_coupon_bond_option(...) */

/* ---------------------------------------------------------------- */

static Err int_range(double *a, double *b, double pay_time, TermStruct *ts,
                     double (*g_x)(double)) {
  /*
  Computes integration range for simpson
          1) Find range which covers density out to 5 S.D.'s
          2) Change one of the bounds to put it on the strike
  */
  Err err = NULL;
  double beta, eta;
  double zeta;
  double y;
  double bump;
  double f_x, d_f_x;
  double i;
  int NotFound;
  double mean;

  beta = find_beta(pay_time, ts);
  eta = find_eta(pay_time, ts);
  zeta = Zeta_func(pay_time, ts);

  /* Compute mean of the density i.e. y in densityXX */
  if (eta != 1.) { /* mean is the same except for power = 1 */
    mean = 1. / (beta * (1. - eta));
  } else
    mean = .5 * beta * zeta;

  /* Compute range which it out to 5 S.D.'s */
  *a = mean - 5 * sqrt(zeta);
  *b = mean + 5 * sqrt(zeta);

  if (*a > *b)
    return ("Error");

  if ((eta != 1.) && (*a < 0.))
    *a = .00001;

  /* Turn off multiplying by the density in swaption_cfs */
  STATIC_density_off = 1;

  /* Move the appropriate bound onto the strike using Newton */
  y = mean;

  bump = .0000001;
  i = 0;
  NotFound = 1;
  while ((i < 10) && (NotFound)) { /* Newton loop */
    f_x = g_x(y);
    d_f_x = (g_x(y + bump) - f_x) / bump;
    y += -f_x / d_f_x;

    /* Check if in range with new y */
    f_x = g_x(y);
    if (fabs(f_x) < .000001)
      NotFound = 0;
    i++;
  }

  if ((!NotFound) && (y < *b) && (y > *a)) { /* Newton converged */
    if (STATIC_pay_rec == 0)
      *a = y;
    else
      *b = y;
  }
  /* else just use a  ,b as 5 S.D.'s in either direction */

  /* Turn on multiplying by the density in swaption_cfs */
  STATIC_density_off = 0;

  return err;
}

/* ---------------------------------------------------------------- */

/* Returns the value of the swaption*density in state y_bar */
static double swaption_cfs(double y_bar) {
  Err err = NULL;
  int i;
  /* model vars */
  double **M;
  double lambda_t;
  double zeta_t, zeta_0;
  double A_t;
  double M_txT;
  double s, theta;
  double xb_yb;
  double nu;
  /* swaption values */
  double intrinsic;
  double answer;
  double beta, eta;

  /* Compute intrisic value */
  intrinsic = 0.0;
  for (i = 1; i <= STATIC_n;
       i++) { /* loop over cash flows after the initial reverse flow */

    /* compute model parameters */
    beta = find_beta(STATIC_pay_time[i], STATIC_ts);
    eta = find_eta(STATIC_pay_time[i], STATIC_ts);
    M = find_M_eta_beta(STATIC_pay_time[i], STATIC_ts);
    if (eta != 1.) {
      nu = 1. / (2. * (1. - eta));
      xb_yb = (pow(fabs(beta * (1 - eta) * y_bar), 2 * nu) - 1.) / beta;
    } else
      xb_yb = (exp(beta * y_bar) - 1.) / beta;
    lambda_t = Psi_func(STATIC_pay_time[i], STATIC_ts);
    zeta_t = Zeta_func(STATIC_pay_time[i], STATIC_ts);

    s = beta * beta * zeta_t / pow(1. + beta * y_bar, 2 * (1. - eta));
    theta = lambda_t * (1. + beta * y_bar) / beta;
    A_t = M_eta_beta_func(s, theta, eta, M);

    zeta_0 = Zeta_func(STATIC_fixing_time, STATIC_ts); /* option expiry */
    s = beta * beta * (zeta_t - zeta_0) /
        pow(1. + beta * y_bar, 2 * (1. - eta));
    theta = lambda_t * (1. + beta * y_bar) / beta;
    M_txT = M_eta_beta_func(s, theta, eta, M);

    /* value cash flow over numeraire */
    intrinsic += STATIC_rp * STATIC_coupon[i] * STATIC_df[i] *
                 exp(-lambda_t * xb_yb - A_t + M_txT);

  } /* END of loop on cash flow dates */

  /* compute model parameters for initial cash flow (strike) */
  beta = find_beta(STATIC_pay_time[0], STATIC_ts);
  eta = find_eta(STATIC_pay_time[0], STATIC_ts);
  M = find_M_eta_beta(STATIC_pay_time[0], STATIC_ts);
  if (eta != 1.) {
    nu = 1. / (2. * (1. - eta));
    xb_yb = (pow(fabs(beta * (1 - eta) * y_bar), 2 * nu) - 1.) / beta;
  } else
    xb_yb = (exp(beta * y_bar) - 1.) / beta;
  lambda_t = Psi_func(STATIC_pay_time[0], STATIC_ts);
  zeta_t = Zeta_func(STATIC_pay_time[0], STATIC_ts);

  s = beta * beta * zeta_t / pow(1. + beta * y_bar, 2 * (1. - eta));
  theta = lambda_t * (1. + beta * y_bar) / beta;
  A_t = M_eta_beta_func(s, theta, eta, M);

  zeta_0 = Zeta_func(STATIC_fixing_time, STATIC_ts);
  s = beta * beta * (zeta_t - zeta_0) / pow(1. + beta * y_bar, 2 * (1. - eta));
  theta = lambda_t * (1. + beta * y_bar) / beta;
  M_txT = M_eta_beta_func(s, theta, eta, M);

  /* value cash flow over numeraire */
  intrinsic -= STATIC_rp * STATIC_bond_strike * STATIC_df[0] *
               exp(-lambda_t * xb_yb - A_t + M_txT);

  /* Set static variables for use in density computation */
  STATIC_x_k = 0.;
  STATIC_deltau = zeta_t;
  STATIC_lambda = lambda_t;
  STATIC_beta = beta;
  STATIC_eta = eta;

  /* intrinsic times density at state y_bar */

  if (STATIC_pay_time[0] == 0.0)
    return intrinsic;

  if (eta == 0.)
    if (STATIC_density_off)
      answer = intrinsic;
    else {
      if (intrinsic < 0.0)
        intrinsic = 0.0;
      answer = intrinsic * density05(y_bar);
    }

  if (eta == .5)
    if (STATIC_density_off)
      answer = intrinsic;
    else {
      if (intrinsic < 0.0)
        intrinsic = 0.0;
      answer = intrinsic * density51(y_bar);
    }

  if (eta == 1.)
    if (STATIC_density_off)
      answer = intrinsic;
    else {
      if (intrinsic < 0.0)
        intrinsic = 0.0;
      answer = intrinsic * density1(y_bar);
    }

  return (answer);
}

/* ---------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
        For the pricing of a cap/floor in betaeta (with spreads already in the
   coupon)
   ----------------------------------------------------------------------------
 */

double srt_f_etabeta_capfloor(
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
      caplet = betaeta_discount_bond_option(
          fixing_time[i], period_time[i], period_time[i + 1], df[i], df[i + 1],
          1.0 / (1.0 + payment[i + 1]), ts, rec_pay, mdl_dim);
      caplet *= 1.0 + payment[i + 1];
      price += caplet;
    }
  }

  return price;

} /* END srt_f_etabeta_capfloor (...) */

/* ----------------------------------------------------------------------------
 */

/* ---------------------------------------------------------------- */

/* ------------------------------------------------------------------------------
 */
/* Price an option on a discount bond (zero coupon)  , using one factor beta-eta
 * model */

static double betaeta_discount_bond_option(
    double fixing_time, double period_start_time, double period_end_time,
    double period_start_df, double period_end_df, double bond_strike,
    TermStruct *ts, SrtReceiverType pay_rec, SrtMdlDim mdl_dim) {

  double a, b;      /* integration bounds */
  double precision; /* precision for simpson */
  double beta, zeta;
  double price;
  Err err = NULL;

  price = 0.;

  /* Check if option has expired */
  if (period_start_time < 0.0)
    return 0.0;

  STATIC_rp = (pay_rec == SRT_RECEIVER ? 1 : -1);

  /* Assign values to static global variables for simpson integration */

  STATIC_fixing_time = fixing_time;
  STATIC_period_start_time = period_start_time;
  STATIC_period_end_time = period_end_time;
  STATIC_period_start_df = period_start_df;
  STATIC_period_end_df = period_end_df;
  STATIC_bond_strike = bond_strike;
  STATIC_ts = ts;
  STATIC_pay_rec = pay_rec;
  STATIC_mdl_dim = mdl_dim;

  /* Assign integration range */

  beta = find_beta(STATIC_period_end_time, ts);
  zeta = Zeta_func(STATIC_period_end_time, ts);

  /* perform newton to find lower integration bound */
  err = int_range(&a, &b, period_start_time, ts, capfl_cfs);
  /* Need to check for errors
  if (err)
  {
  }
  */

  /* Set precision for simpson  */
  precision = .001;

  /* Compute expected value by integrating against the density */
  if (a < b)
    price = sm_qsimp(capfl_cfs, a, b, precision);
  else
    price = 0;

  return price;
}

/* ------------------------------------------------------------------------------
 */

/* Returns the value of the swaption*density in state y_bar */
static double capfl_cfs(double y_bar) {
  Err err = NULL;

  /* model vars */
  double **M;
  double lambda_t;
  double lambda_T;
  double zeta_t, zeta_0;
  double A_t;
  double A_T;
  double M_txT;
  double s, theta;
  double xb_yb;
  double nu;
  /* swaption values */
  double intrinsic;
  double answer;
  double beta, eta;

  /* Compute intrisic value */
  intrinsic = 0.0;

  zeta_0 = Zeta_func(STATIC_fixing_time, STATIC_ts);
  beta = find_beta(STATIC_period_end_time, STATIC_ts);
  eta = find_eta(STATIC_period_end_time, STATIC_ts);
  if (eta != 1.) {
    nu = 1. / (2. * (1. - eta));
    xb_yb = (pow(fabs(beta * (1 - eta) * y_bar), 2 * nu) - 1.) / beta;
  } else
    xb_yb = (exp(beta * y_bar) - 1.) / beta;

  M = find_M_eta_beta(STATIC_period_end_time, STATIC_ts);

  lambda_T = Psi_func(STATIC_period_end_time, STATIC_ts);
  zeta_t = Zeta_func(STATIC_period_end_time, STATIC_ts);

  s = beta * beta * zeta_t / pow(1. + beta * y_bar, 2 * (1. - eta));
  theta = lambda_T * (1. + beta * y_bar) / beta;
  A_T = M_eta_beta_func(s, theta, eta, M);

  s = beta * beta * (zeta_t - zeta_0) / pow(1. + beta * y_bar, 2 * (1. - eta));
  theta = lambda_T * (1. + beta * y_bar) / beta;
  M_txT = M_eta_beta_func(s, theta, eta, M);

  /* value the forward bond over numeraire */
  intrinsic =
      STATIC_rp * STATIC_period_end_df * exp(-lambda_T * xb_yb - A_T + M_txT);

  /* compute model parameters for initial cash flow (strike) */
  beta = find_beta(STATIC_period_start_time, STATIC_ts);
  eta = find_eta(STATIC_period_start_time, STATIC_ts);
  if (eta != 1) {
    nu = 1. / (2. * (1. - eta));
    xb_yb = (pow(fabs(beta * (1 - eta) * y_bar), 2 * nu) - 1.) / beta;
  } else
    xb_yb = (exp(beta * y_bar) - 1.) / beta;

  M = find_M_eta_beta(STATIC_period_start_time, STATIC_ts);

  lambda_t = Psi_func(STATIC_period_start_time, STATIC_ts);
  zeta_t = Zeta_func(STATIC_period_start_time, STATIC_ts);

  s = beta * beta * zeta_t / pow(1. + beta * y_bar, 2 * (1. - eta));
  theta = lambda_t * (1. + beta * y_bar) / beta;
  A_t = M_eta_beta_func(s, theta, eta, M);

  zeta_0 = Zeta_func(STATIC_fixing_time, STATIC_ts); /* option expiry */
  s = beta * beta * (zeta_t - zeta_0) / pow(1. + beta * y_bar, 2 * (1. - eta));
  theta = lambda_t * (1. + beta * y_bar) / beta;
  M_txT = M_eta_beta_func(s, theta, eta, M);

  /* value cash flow over numeraire for the initial cash flow (strike) */
  intrinsic -= STATIC_rp * STATIC_bond_strike * STATIC_period_start_df *
               exp(-lambda_t * xb_yb - A_t + M_txT);

  if (STATIC_period_start_time == 0.0)
    return intrinsic;

  /* Set static variables for use in density computation */
  STATIC_x_k = 0.;
  STATIC_deltau = Zeta_func(STATIC_period_start_time, STATIC_ts);
  STATIC_lambda = Psi_func(STATIC_period_start_time, STATIC_ts);
  STATIC_beta = beta;
  STATIC_eta = eta;

  /* intrinsic times density at state y_bar */

  if (eta == 0.)
    if (STATIC_density_off)
      answer = intrinsic;
    else {
      if (intrinsic < 0.0)
        intrinsic = 0.0;
      answer = intrinsic * density05(y_bar);
    }

  if (eta == .5)
    if (STATIC_density_off)
      answer = intrinsic;
    else {
      if (intrinsic < 0.0)
        intrinsic = 0.0;
      answer = intrinsic * density51(y_bar);
    }

  if (eta == 1.)
    if (STATIC_density_off)
      answer = intrinsic;
    else {
      if (intrinsic < 0.0)
        intrinsic = 0.0;
      answer = intrinsic * density1(y_bar);
    }

  return (answer);
}

/* ------------------------------------------------------------------------------
 */

static double density05(double y_bar)
/* Computes the integrand for Case 1: 0<=eta<.5  .5 <= nu < 1
   No barrier b.c.'s
 */
/* NB: Case 1a now contains both 1a and 1b so the 1b integrand
   is no longer needed */
{
  double ans;
  double nu;
  double y;
  double xb_yb;
  double dxb_yb;
  double xt_yb;
  double i_nu, i_minus_nu, k_nu;
  Err err = NULL;

  nu = 1. / (2. * (1. - STATIC_eta));
  y = (pow(fabs(1. + STATIC_beta * STATIC_x_k), 1 - STATIC_eta)) /
      (STATIC_beta * (1. - STATIC_eta));
  xb_yb = (pow(fabs(STATIC_beta * (1 - STATIC_eta) * y_bar), 2 * nu) - 1.) /
          STATIC_beta;
  dxb_yb = 2 * nu * pow((1. - STATIC_eta) * STATIC_beta * y_bar, 2 * nu - 1);
  xt_yb = (-pow(fabs(STATIC_beta * (1 - STATIC_eta) * y_bar), 2 * nu) - 1.) /
          STATIC_beta;

  err = I_nu(nu, y * y_bar / STATIC_deltau, &i_nu);
  err = I_nu(nu, y * y_bar / STATIC_deltau, &i_minus_nu);
  err = K_nu(nu, y * y_bar / STATIC_deltau, &k_nu);

  ans = (/* 1a */
         pow(y / y_bar, nu) * (y_bar / STATIC_deltau) *
         (.5 * i_nu + .5 * i_minus_nu) *
         /* exp(-y*y_bar/deltau)*  taken out due to e^{-x} from N.R. I_nu */
         exp(-((y_bar - y) * (y_bar - y)) / (2 * STATIC_deltau))) +
        (/* 1b */
         pow(y / y_bar, nu) * (y_bar / STATIC_deltau) *
         (sin(SRT_PI * nu) / SRT_PI) * k_nu *
         exp(-(y_bar + y) * (y_bar + y) / (2 * STATIC_deltau)));

  return ans;
}
/* ---------------------------------------------------------------- */

static double density51(double y_bar) {

  /* Computes the integrand for Case 3: .5 <=eta< 1   nu >= 1 */

  double y;
  double xb_yb;
  double dxb_yb;
  double nu;
  double ans; /* value of integrand */
  double F;
  double F_xb;
  double i_nu;
  Err err = NULL;

  nu = 1 / (2 * (1 - STATIC_eta));
  y = (pow(fabs(1 + STATIC_beta * STATIC_x_k), 1 - STATIC_eta)) /
      (STATIC_beta * (1 - STATIC_eta));
  xb_yb = (pow(fabs(STATIC_beta * (1 - STATIC_eta) * y_bar), 2 * nu) - 1) /
          STATIC_beta;
  dxb_yb = pow(STATIC_beta * (1 - STATIC_eta) * y_bar, 2 * nu - 1);
  F = exp(-STATIC_lambda * (xb_yb - STATIC_x_k));
  F_xb = -STATIC_lambda * F;

  err = I_nu(nu - 1., y * y_bar / STATIC_deltau, &i_nu);

  ans = ((y_bar / STATIC_deltau) * F - F_xb * dxb_yb) * pow(y_bar / y, 1 - nu) *
        i_nu *
        /*	  exp(-y*y_bar/STATIC_DELTAU)*    taken out due to NR scaling */
        exp(-(y_bar - y) * (y_bar - y) / (2 * STATIC_deltau));

  return ans;
}

/* ---------------------------------------------------------------- */

/* Computes the integrand for Case 4: eta = 1 (bdy at - infty) */
static double density1(double y_bar) {
  double y;
  double xb_yb;
  double ans; /* value of integrand */
  double F;

  y = log(1. + STATIC_beta * STATIC_x_k);
  xb_yb = (exp(STATIC_beta * y_bar) - 1.) / STATIC_beta;
  F = exp(-STATIC_lambda * (xb_yb - STATIC_x_k));

  ans = F *
        exp(-(y_bar - y + .5 * STATIC_beta * STATIC_deltau) *
            (y_bar - y + .5 * STATIC_beta * STATIC_deltau) /
            (2 * STATIC_deltau)) /
        sqrt(2 * SRT_PI * STATIC_deltau);

  return ans;
}

/* =============================================================================
 */