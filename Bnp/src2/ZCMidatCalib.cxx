
#include "AmortMidatCalib.h"
#include "CPDCalib.h"
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
#define CALPRES 1.0e-08
#define NITER 10
#define BUMP_LAM_LM 5.0E-04

#define MAX_FACT 0.25

extern double x[NUM_HERMITE + 1], w[NUM_HERMITE + 1];

double ZCMidat_lgmsopval1F(int ncpn,     /*	Number of cash-flow dates        ,
                                            including     start and end date */
                           double cpn[], /*	Notional */
                           double df[],  /*	Df to cash flow dates */
                           double cvg[], /*	cvg from i-1 to i */
                           double cpn_G1[], /*	G1 at cash-flow dates */
                           double ex_zeta1, /*	Z1 at exercise date */
                           double ex_G1)    /*	G1 at exercise date */
{
  int i;
  static double coupon[MAX_CPN];

  for (i = 0; i < ncpn; i++) {
    coupon[i] = cpn[i];
  }

  return lgmopval1F(ncpn, coupon, cpn_G1, ex_zeta1, ex_G1);
}

double ZCMidat_lgmsopval2F(int ncpn,     /*	Number of cash-flow dates        ,
                                            including     start and end date */
                           double cpn[], /*	Notional */
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
  static double coupon[MAX_CPN];

  for (i = 0; i < ncpn; i++) {
    coupon[i] = cpn[i];
  }

  return lgmopval2F(ncpn, coupon, cpn_G1, cpn_G2, ex_zeta1, ex_zeta2, ex_zeta12,
                    ex_G1, ex_G2);
}

double ZCMidat_lgmZCsopval2F(int ncpn,         /*	Number of cash-flow dates         ,
                                            including         start and end date */
                             double not [],    /*	Notional */
                             double df[],      /*	Df to cash flow dates */
                             double cvg[],     /*	cvg from i-1 to i */
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

  cpn[0] = -not [0] * df[0];

  for (i = 1; i < ncpn; i++) {
    cpn[i] = not [i - 1] * df[i] * cvg[i] * strike;
  }

  cpn[ncpn - 1] += not [ncpn - 2] * df[ncpn - 1];

  return lgmopval2F(ncpn, cpn, cpn_G1, cpn_G2, ex_zeta1, ex_zeta2, ex_zeta12,
                    ex_G1, ex_G2);
}

double ZCMidat_lgmZCopval2F(int ncpn,         /*	Number of cash-flow dates         ,
                                           including         start and end date */
                            double df[],      /*	Df to cash flow dates */
                            double cvg[],     /*	cvg from i-1 to i */
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
  double sum;
  double Beta;

  sum = 0.0;
  Beta = 1.0;

  cpn[0] = -df[0];

  for (i = 1; i < ncpn; i++) {
    Beta = Beta * (1 + cvg[i] * strike);
    sum = cvg[i] * strike * Beta;
  }

  cpn[ncpn - 1] = sum * df[ncpn - 1];

  return lgmopval2F(ncpn, cpn, cpn_G1, cpn_G2, ex_zeta1, ex_zeta2, ex_zeta12,
                    ex_G1, ex_G2);
}

static int NCPN_LM, NEX_LM, *EX_CPN_LM, *EX_LENDCPN_LM, *EX_SENDCPN_LM,
    ONE2F_LM, SKIP_LAST_LM, NB_MOM_LM, FREQ_SHORT_LM, SHIFT_FREQ_LM;

static double *CPN_TIME_LM, *CPN_DF_LM, *CPN_CVG_LM, *EX_TIME_LM,
    *EX_LSTRIKE_LM, *EX_LPRICE_LM, *EX_SSTRIKE_LM, *EX_SPRICE_LM, *EX_ZETA_LM,
    *LAM_TIME_LM, LAM_LM[MAX_CPN], ALPHA_LM, GAMMA_LM, RHO_LM,
    CPN_G_LM[MAX_CPN], CPN_G2_LM[MAX_CPN], EX_G_LM[MAX_CPN], EX_G2_LM[MAX_CPN],
    ZETA2_LM[MAX_CPN], ZETA12_LM[MAX_CPN], PRICE_LM[MAX_CPN],
    SENSI_LM[MAX_CPN][MAX_CPN];

/*	Calibrate zeta to diagonal given G: 1F case */
Err ZCMidat_lgmcalibzeta1F(
    int ncpn,               /*	Total number of cash-flow dates */
    double cpn[],           /*	Discounted Cash-Flow */
    double cpn_time[],      /*	Cash-Flow times */
    double cpn_df[],        /*	Df to cash-flow dates */
    double cpn_cvg[],       /*	cvg from i-1 to i */
    double cpn_G[],         /*	G at cash-flow dates */
    int nex,                /*	Total number of exercise dates */
    double ex_time[],       /*	Exercise times */
    int ex_cpn[],           /*	Index of the first cash-flow to be exercised */
    double ex_G[],          /*	G at exercise date */
    double strike[],        /*	Strikes */
    double mkt_price[],     /*	Market prices */
    double ex_zeta[],       /*	Output: zetas */
    double exerFixCoupon[], /*	Exercise Fix Coupon	*/
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
  static double coupon[MAX_CPN];

  for (i = 0; i < ncpn; ++i) {
    coupon[i] = cpn[i];
  }

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

    coupon[j] = -cpn_df[j];
    coupon[ncpn - j - 1] =
        coupon[ncpn - j - 1] + exerFixCoupon[j] * cpn_df[ncpn - j - 1];

    pr1 = amortMidat_lgmamortsopval1F(ncpn - j, coupon + j, cpn_df + j,
                                      cpn_cvg + j, cpn_G + j, z1, ex_G[i]);

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

    pr2 = amortMidat_lgmamortsopval1F(ncpn - j, coupon + j, cpn_df + j,
                                      cpn_cvg + j, cpn_G + j, z2, ex_G[i]);

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

      pr2 = amortMidat_lgmamortsopval1F(ncpn - j, coupon + j, cpn_df + j,
                                        cpn_cvg + j, cpn_G + j, z2, ex_G[i]);

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

/*	Calibrate zeta to ZC diagonal swaption given G: 2F case */
Err ZCMidat_lgmcalibzeta2F(
    int ncpn,               /*	Total number of cash-flow dates */
    double cpn[],           /*	Discounted Cash-Flow */
    double cpn_time[],      /*	Cash-Flow times */
    double cpn_df[],        /*	Df to cash-flow dates */
    double cpn_cvg[],       /*	cvg from i-1 to i */
    double cpn_G1[],        /*	G1 at cash-flow dates */
    double cpn_G2[],        /*	G2 at cash-flow dates */
    int nex,                /*	Total number of exercise dates */
    double ex_time[],       /*	Exercise times */
    int ex_cpn[],           /*	Index of the first cash-flow to be exercised */
    double ex_G1[],         /*	G1 at exercise date */
    double ex_G2[],         /*	G2 at exercise date */
    double strike[],        /*	Strikes */
    double mkt_price[],     /*	Market prices */
    double ex_zeta[],       /*	Output: zetas (1) */
    double exerFixCoupon[], /*	Float Notional	*/
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
  static double coupon[MAX_CPN];

  for (i = 0; i < ncpn; ++i) {
    coupon[i] = cpn[i];
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

    coupon[j] = -cpn_df[j];
    coupon[ncpn - j - 1] =
        coupon[ncpn - j - 1] + exerFixCoupon[j] * cpn_df[ncpn - j - 1];

    pr1 = amortMidat_lgmamortsopval2F(ncpn - j, coupon + j, cpn_df + j,
                                      cpn_cvg + j, cpn_G1 + j, cpn_G2 + j, z1,
                                      z2, z12, ex_G1[i], ex_G2[i]);

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

    pr2 = amortMidat_lgmamortsopval2F(ncpn - j, coupon + j, cpn_df + j,
                                      cpn_cvg + j, cpn_G1 + j, cpn_G2 + j, q1,
                                      q2, q12, ex_G1[i], ex_G2[i]);

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

      pr2 = amortMidat_lgmamortsopval2F(ncpn - j, coupon + j, cpn_df + j,
                                        cpn_cvg + j, cpn_G1 + j, cpn_G2 + j, q1,
                                        q2, q12, ex_G1[i], ex_G2[i]);

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
Err ZCMidat_lgmprcapgivenlambda(
    int ncpn,            /*	Total number of cash-flow dates */
    double cpn[],        /*	Discounted Cash-Flows */
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
    double ex_sweight[], double ex_zeta[], /*	Output: zetas */
    double exerFixCoupon[],                /*	Float Notional	*/
    double lambda,                         /*	Lambda */
    int one2F,                             /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    double alpha, double gamma, double rho,
    int skip_last,     /*	If 1        , the last option is disregarded
                                         and the forward volatility is flat
                    from option     n-1 */
    int price_cap,     /*	0: just calibrate */
    double *ex_sprice) /*	Cap price as output */
{
  static double cpn_G[MAX_CPN], cpn_G2[MAX_CPN], ex_G[MAX_CPN], ex_G2[MAX_CPN],
      zeta2[MAX_CPN], zeta12[MAX_CPN];
  Err err;

  /*	Setup G function */
  export_lgmsetupG(lambda, ncpn, cpn_time, cpn_G, nex, ex_time, ex_G);

  if (one2F == 2) {
    export_lgmsetupG2(lambda, gamma, ncpn, cpn_time, cpn_G2, nex, ex_time,
                      ex_G2);
  }

  if (one2F == 1) {
    err = ZCMidat_lgmcalibzeta1F(
        ncpn, cpn, cpn_time, cpn_df, cpn_cvg, cpn_G, nex, ex_time, ex_cpn, ex_G,
        ex_lstrike, ex_lprice, ex_zeta, exerFixCoupon, lambda, skip_last);

    if (err) {
      return err;
    }

    if (price_cap) {
      *ex_sprice = lgmcapval1F(ncpn, cpn_df, cpn_cvg, cpn_G, nex, ex_cpn,
                               ex_sncpn, ex_sweight, ex_zeta, ex_G, ex_sstrike);
    }
  } else {
    err = ZCMidat_lgmcalibzeta2F(ncpn, cpn, cpn_time, cpn_df, cpn_cvg, cpn_G,
                                 cpn_G2, nex, ex_time, ex_cpn, ex_G, ex_G2,
                                 ex_lstrike, ex_lprice, ex_zeta, exerFixCoupon,
                                 lambda, alpha, gamma, rho, skip_last);

    if (err) {
      return err;
    }

    export_lgmcalczeta2zeta12(nex, ex_time, ex_zeta, lambda, alpha, gamma, rho,
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
Err ZCMidat_lgmcalibzetalambda(
    int ncpn,            /*	Total number of cash-flow dates */
    double cpn[],        /*	Discounted Cash-Flows */
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
    double ex_sweight[], double ex_sprice, /*	Market price for cap */
    double ex_zeta[],                      /*	Output: zetas */
    double exerFixCoupon[],                /*	Float Notional	*/
    int fix_lambda, /*	0: calib lambda to cap        , 1: fix lambda calib
                                      to diagonal */
    double *lambda, /*	Lambda: may be changed in the process */
    int one2F,      /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    double alpha, double gamma, double rho,
    int skip_last) /*	If 1        , the last option is disregarded
                                           and the forward volatility is flat
                      from option n-1 */
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

  err = ZCMidat_lgmprcapgivenlambda(
      ncpn, cpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time, ex_cpn, ex_sncpn,
      ex_lstrike, ex_lprice, ex_sstrike, ex_sweight, ex_zeta, exerFixCoupon,
      *lambda, one2F, alpha, gamma, rho, skip_last, (fix_lambda ? 0 : 1), &pr1);

  if (err || fix_lambda) {
    return err;
  }

  if (fabs(ex_sprice - pr1) < CALPRES) {
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

    if (err = ZCMidat_lgmprcapgivenlambda(
            ncpn, cpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time, ex_cpn,
            ex_sncpn, ex_lstrike, ex_lprice, ex_sstrike, ex_sweight, ex_zeta,
            exerFixCoupon, lam1, one2F, alpha, gamma, rho, skip_last, 1,
            &pr1)) {
      return err;
    }

    if (fabs(ex_sprice - pr1) < CALPRES) {
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
    lam2 = 0.5 * (lam1 + lam2);

    if (err = ZCMidat_lgmprcapgivenlambda(
            ncpn, cpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time, ex_cpn,
            ex_sncpn, ex_lstrike, ex_lprice, ex_sstrike, ex_sweight, ex_zeta,
            exerFixCoupon, lam2, one2F, alpha, gamma, rho, skip_last, 1,
            &pr2)) {
      return err;
    }

    if (fabs(ex_sprice - pr2) < CALPRES) {
      *lambda = lam2;
      return NULL;
    }

    if (pr2 < ex_sprice) {
      pr1 = pr2;
      lam1 = lam2;
      lam2 = temp;
      pr2 = temppr;
    }

    if (fabs(lam2 - lam1) < CALPRES || fabs(pr2 - pr1) < CALPRES)
      break;

    while (pr1 > pr2 && it < 10) {
      it++;

      temp = lam1;
      lam1 = lam2;
      pr1 = pr2;
      lam2 += lam2 - temp;

      if (err = ZCMidat_lgmprcapgivenlambda(
              ncpn, cpn, cpn_time, cpn_df, cpn_cvg, nex, ex_time, ex_cpn,
              ex_sncpn, ex_lstrike, ex_lprice, ex_sstrike, ex_sweight, ex_zeta,
              exerFixCoupon, lam2, one2F, alpha, gamma, rho, skip_last, 1,
              &pr2)) {
        return err;
      }

      if (fabs(ex_sprice - pr2) < CALPRES) {
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

/*	Calibrate lgm: main function */
Err ZCMidat_cpd_calib_diagonal(
    char *yc_name,               /*	Name of the yield curve */
    char *vol_curve_name,        /*	Name of the market vol curve */
    char *default_ref_rate_name, /*	Name of the reference rate */
    Err (*get_cash_vol)(         /*	Function to get cash vol from the market */
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
    long start_date,                  /*	Start date of the amortised swap */
    long end_date,                    /*	End date for the amortised swap */
    double *long_strike,              /*	Diagonal swaption strikes
                                                                        NULL = ATM */
    double *short_strike,             /*	Short swaption strikes
                                                                NULL = ATM */
    int strike_type,                  /*	0: ATM
                                                        1: CASH
                                                        2: SWAP
                                                        3: STD */
    double *diag_prices,              /* Diagonal Prices */
    double max_std_long, double max_std_short, char *ref_rate_name,
    char *swaption_freq, /*	Frequency        , basis and ref. rate of
                      underlying swaptions */
    char *swaption_basis, int fix_lambda, /*	0: calib lambda to cap        ,
                                 1: fix lambda calib to diagonal */
    int one_f_equi,                       /*	1F equivalent flag:
                                                            if set to 1        , then 2F lambda will
                                       calibrate                       to the cap priced within
                                       calibrated 1F                       with the given lambda */
    int skip_last,  /*	If 1        , the last option is disregarded
                                      and the forward volatility is flat
                 from option  n-1 */
    double *lambda, /*	Lambda: may be changed in the process */
    int one2F,      /*	Number of factors */
    /*	Alpha        , Gamma        , Rho (2F only) */
    double alpha, double gamma, double rho, int *num_sig, /*	Answer */
    double **sig_time, double **sig,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA
        inst_data) /*	NULL = don't save calibration instrument data */
{
  int i, j, k, ncpn, n, float_index, fix_index;
  int nbcoupon;
  SrtCompounding ifreq;
  SrtBasisCode ibasis;
  double ex_lfwd[MAX_CPN], ex_llvl[MAX_CPN], ex_lstrike[MAX_CPN],
      ex_lvol[MAX_CPN], ex_lprice[MAX_CPN], ex_sfwd[MAX_CPN], ex_slvl[MAX_CPN],
      ex_sstrike[MAX_CPN], ex_svol[MAX_CPN], ex_sweight[MAX_CPN],
      ex_sprice[MAX_CPN], ex_zeta[MAX_CPN], ex_G[MAX_CPN];
  int ex_sncpn[MAX_CPN];
  double cpn_G[MAX_CPN];
  long theo_date, act_date;
  long today;
  double lvl, dfi, dff;
  double power;
  double swp_rte, spr;
  double cap_price;
  double std;
  SrtCurvePtr yc_ptr;
  Err err = NULL;
  double coupon[MAX_CPN];
  double coupon_time[MAX_CPN];
  double coupon_df[MAX_CPN];
  double coupon_cvg[MAX_CPN];
  int ex_coupon[MAX_CPN];

  double temp;

  //---------------------------------------
  //------------SWAPDP---------------------
  //---------------------------------------
  SwapDP swapdp;
  long float_nb_dates, float_nb_pay_dates;
  long *float_fixing_dates = NULL, *float_start_dates = NULL,
       *float_end_dates = NULL, *float_pay_dates = NULL;
  double *float_cvgs = NULL, *float_spreads = NULL, *float_pay_times = NULL;
  long fix_nb_dates, fix_nb_pay_dates;
  long *fix_start_dates = NULL, *fix_end_dates = NULL, *fix_pay_dates = NULL;
  double *fix_cvgs = NULL, *fix_pay_times = NULL;

  double *coupon_dates = NULL;

  int nbexercise;
  double exer_time[MAX_CPN];
  long exer_date[MAX_CPN];

  double Cumcoupon;

  double *exerFixCoupon = NULL;

  //---------------------------------------
  //---------------------------------------

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
    err = "Not enough coupons in amortMidat_cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  theo_date = end_date;
  act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
  i = ncpn - 1;

  //--------------------------------------------
  //---------SWAPDP-----------------------------
  //---------Build Swap Schedule----------------
  //--------------------------------------------

  err = swp_f_initSwapDP(start_date, theo_date, swaption_freq, swaption_basis,
                         &swapdp);
  if (err) {
    goto FREE_RETURN;
  }

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &swapdp, today, ref_rate_name, &float_pay_dates, &float_nb_pay_dates,
      &float_fixing_dates, &float_start_dates, &float_end_dates, &float_cvgs,
      &float_spreads, &float_nb_dates);
  if (err) {
    goto FREE_RETURN;
  }
  float_pay_times = dvector(0, float_nb_pay_dates - 1);
  for (n = 0; n < float_nb_pay_dates; ++n) {
    float_pay_times[n] = (float_pay_dates[n] - today) * YEARS_IN_DAY;
  }

  err = swp_f_make_FixedLegDatesAndCoverages(
      &swapdp, today, &fix_pay_dates, &fix_nb_pay_dates, &fix_start_dates,
      &fix_end_dates, &fix_cvgs, &fix_nb_dates);
  if (err) {
    goto FREE_RETURN;
  }
  fix_pay_times = dvector(0, fix_nb_pay_dates - 1);
  for (n = 0; n < fix_nb_pay_dates; ++n) {
    fix_pay_times[n] = (fix_pay_dates[n] - today) * YEARS_IN_DAY;
  }

  nbcoupon = float_nb_pay_dates;

  coupon_dates = (double *)calloc(nbcoupon, sizeof(double));
  memcpy(coupon_dates, float_pay_times, nbcoupon * sizeof(double));
  num_f_concat_vector(&nbcoupon, &coupon_dates, fix_nb_pay_dates,
                      fix_pay_times);
  num_f_sort_vector(nbcoupon, coupon_dates);
  num_f_unique_vector(&nbcoupon, coupon_dates);

  for (k = 0; k < nbcoupon; ++k) {
    coupon[k] = 0;
  }

  Cumcoupon = 1.0;
  nbexercise = fix_nb_pay_dates - 1;
  float_index = 1;
  fix_index = 1;
  coupon[0] = -swp_f_df(today, float_pay_dates[0], yc_name);
  coupon_time[0] = coupon_dates[0];
  coupon_df[0] = swp_f_df(today, float_pay_dates[0], yc_name);
  coupon_cvg[0] = 0;
  for (k = 1; k < nbcoupon; ++k) {
    coupon_time[k] = coupon_dates[k];
    coupon_df[k] = swp_f_df(today, float_pay_dates[k], yc_name);
    coupon_cvg[k] = float_cvgs[k - 1];
    if ((coupon_dates[k] == fix_pay_times[fix_index]) &&
        (coupon_dates[k] != float_pay_times[float_index])) {
      Cumcoupon = Cumcoupon *
                  (1 + long_strike[fix_index - 1] * fix_cvgs[fix_index - 1]);
      Cumcoupon += long_strike[fix_index - 1];
      //			coupon[k] += long_strike[fix_index-1] *
      // fix_cvgs[fix_index-1] * swp_f_df (today        ,
      // fix_pay_dates[fix_index]        , yc_name);
      ex_coupon[fix_index - 1] =
          k - (float_nb_pay_dates - 1) / (fix_nb_pay_dates - 1);
      exer_time[fix_index - 1] = fix_pay_times[fix_index - 1];
      exer_date[fix_index - 1] = fix_pay_dates[fix_index - 1];
      fix_index = fix_index + 1;
    } else if ((coupon_dates[k] != fix_pay_times[fix_index]) &&
               (coupon_dates[k] == float_pay_times[float_index])) {
      temp = swp_f_df(today, float_pay_dates[float_index], yc_name);
      coupon[k] +=
          -(float_spreads[float_index - 1] * float_cvgs[float_index - 1]) *
          swp_f_df(today, float_pay_dates[float_index], yc_name);
      float_index = float_index + 1;
    } else {
      if (float_index < float_nb_pay_dates - 1) {
        Cumcoupon = Cumcoupon *
                    (1 + long_strike[fix_index - 1] * fix_cvgs[fix_index - 1]);
        Cumcoupon += long_strike[fix_index - 1];
        //				coupon[k] += long_strike[fix_index-1] *
        // fix_cvgs[fix_index-1] * swp_f_df (today        ,
        // fix_pay_dates[fix_index]        , yc_name);
        ex_coupon[fix_index - 1] =
            k - (float_nb_pay_dates - 1) / (fix_nb_pay_dates - 1);
        exer_time[fix_index - 1] = fix_pay_times[fix_index - 1];
        exer_date[fix_index - 1] = fix_pay_dates[fix_index - 1];
        fix_index = fix_index + 1;
        coupon[k] +=
            -(float_spreads[float_index - 1] * float_cvgs[float_index - 1]) *
            swp_f_df(today, float_pay_dates[float_index], yc_name);
        float_index = float_index + 1;
      } else {
        Cumcoupon = Cumcoupon *
                    (1 + long_strike[fix_index - 1] * fix_cvgs[fix_index - 1]);
        Cumcoupon += long_strike[fix_index - 1];
        //				coupon[k] += long_strike[fix_index-1] *
        // fix_cvgs[fix_index-1] * swp_f_df (today        ,
        // fix_pay_dates[fix_index]        , yc_name);
        ex_coupon[fix_index - 1] =
            k - (float_nb_pay_dates - 1) / (fix_nb_pay_dates - 1);
        exer_time[fix_index - 1] = fix_pay_times[fix_index - 1];
        exer_date[fix_index - 1] = fix_pay_dates[fix_index - 1];
        fix_index = fix_index + 1;
        coupon[k] += -(-1 + float_spreads[float_index - 1] *
                                float_cvgs[float_index - 1]) *
                     swp_f_df(today, float_pay_dates[float_index], yc_name);
        float_index = float_index + 1;
      }
    }
  }

  if (ex_date) {
    nbexercise = num_ex_dates;
    for (i = 0; i < num_ex_dates; ++i) {
      exer_date[i] = ex_date[i];
      exer_time[i] = (ex_date[i] - today) * YEARS_IN_DAY;
    }
  }

  exerFixCoupon = (double *)calloc(fix_nb_pay_dates - 1, sizeof(double));
  temp = 1;
  for (k = fix_nb_pay_dates - 2; k >= 0; --k) {
    temp = temp * (1 + long_strike[k] * fix_cvgs[k]);
    exerFixCoupon[k] = temp;
  }

  //-----------------------------------------------------------------
  //----------------End Of Build Swap Schedule-----------------------
  //-----------------------------------------------------------------

  /*	Underlyings */

  /*	Long */

  dff = swp_f_df(today, float_pay_dates[nbcoupon - 1], yc_name);
  for (i = 0; i < nbexercise; i++) {
    j = ex_coupon[i];

    lvl = 0.0;
    for (k = j + 1; k < nbcoupon; k++) {
      lvl += coupon_cvg[k] * coupon_df[k];
    }
    dfi = swp_f_df(today, float_pay_dates[j], yc_name);

    ex_llvl[i] = lvl;
    ex_lfwd[i] = (dfi - dff) / lvl;

    /*	ATM std */
    err = get_cash_vol(vol_curve_name,
                       add_unit(exer_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                       float_pay_dates[nbcoupon - 1], ex_lfwd[i], 0,
                       ref_rate_name, &std, &power);
    if (err) {
      goto FREE_RETURN;
    }
    std += (shift_type == 1 ? std * vol_shift : vol_shift);
    if (power > 0.5) {
      power = srt_f_optblksch(ex_lfwd[i], ex_lfwd[i], std, exer_time[i], 1.0,
                              SRT_CALL, PREMIUM);
      err = srt_f_optimpvol(power, ex_lfwd[i], ex_lfwd[i], exer_time[i], 1.0,
                            SRT_CALL, SRT_NORMAL, &std);
    }
    std *= sqrt(exer_time[i]);

    /*	Strike */
    if ((!long_strike) || (!strike_type)) {
      ex_lstrike[i] = ex_lfwd[i];
    } else if (strike_type == 1) {
      ex_lstrike[i] = long_strike[i];
    } else if (strike_type == 2) {
      if (err = swp_f_ForwardRate(
              float_pay_dates[j], float_pay_dates[nbcoupon - 1], swaption_freq,
              swaption_basis, yc_name, ref_rate_name, &swp_rte)) {
        goto FREE_RETURN;
      }

      spr = swp_rte - ex_lfwd[i];

      ex_lstrike[i] = long_strike[i];
      //			ex_lstrike[i] = long_strike[i] - spr;
    } else if (strike_type == 3) {
      ex_lstrike[i] = ex_lfwd[i] + long_strike[i] * std;
    }

    /*	Apply max std */
    if (ex_lstrike[i] > ex_lfwd[i] + max_std_long * std) {
      //			ex_lstrike[i] = ex_lfwd[i] + max_std_long * std;
    } else if (ex_lstrike[i] < ex_lfwd[i] - max_std_long * std) {
      //			ex_lstrike[i] = ex_lfwd[i] - max_std_long * std;
    }

    /*	Make sure strikes are positive (actually more than 1bp)
                    otherwise use ATM	*/
    if (ex_lstrike[i] < 1.0e-04) {
      ex_lstrike[i] = ex_lfwd[i];
    }

    err = get_cash_vol(vol_curve_name,
                       add_unit(exer_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                       float_pay_dates[nbcoupon - 1], ex_lstrike[i], 0,
                       ref_rate_name, &(ex_lvol[i]), &power);
    if (err) {
      goto FREE_RETURN;
    }
    ex_lvol[i] += (shift_type == 1 ? ex_lvol[i] * vol_shift : vol_shift);

    if (power > 0.5) {
      ex_lprice[i] =
          srt_f_optblksch(ex_lfwd[i], ex_lstrike[i], ex_lvol[i], exer_time[i],
                          ex_llvl[i], SRT_PUT, PREMIUM);
    } else {
      ex_lprice[i] =
          srt_f_optblknrm(ex_lfwd[i], ex_lstrike[i], ex_lvol[i], exer_time[i],
                          ex_llvl[i], SRT_PUT, PREMIUM);
    }
  }

  /*	Short */

  for (i = 0; i < nbexercise; i++) {
    ex_sweight[i] = 1.0;
  }

  if (!fix_lambda) {
    cap_price = 0.0;
    for (i = 0; i < nbexercise; i++) {
      if (i < nbexercise - 1) {
        ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
      } else {
        ex_sncpn[i] = nbexercise - ex_coupon[i];
      }

      if (ex_sncpn[i] < 2) {
        err = "One exercise date controls less than 2 coupons in "
              "cpd_calib_diagonal";
        goto FREE_RETURN;
      }

      lvl = 0.0;
      for (k = ex_coupon[i] + 1; k < ex_coupon[i] + ex_sncpn[i]; k++) {
        lvl += coupon_cvg[k] * coupon_df[k];
      }
      dfi = swp_f_df(today, float_pay_dates[ex_coupon[i]], yc_name);
      dff = swp_f_df(today, float_pay_dates[ex_coupon[i] + ex_sncpn[i] - 1],
                     yc_name);

      ex_slvl[i] = lvl;
      ex_sfwd[i] = (dfi - dff) / lvl;

      /*	ATM std */
      err =
          get_cash_vol(vol_curve_name,
                       add_unit(exer_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                       float_pay_dates[ex_coupon[i] + ex_sncpn[i] - 1],
                       ex_sfwd[i], 0, ref_rate_name, &std, &power);
      if (err) {
        goto FREE_RETURN;
      }
      std += (shift_type == 1 ? std * vol_shift : vol_shift);
      if (power > 0.5) {
        power = srt_f_optblksch(ex_sfwd[i], ex_sfwd[i], std, exer_time[i], 1.0,
                                SRT_CALL, PREMIUM);
        err = srt_f_optimpvol(power, ex_sfwd[i], ex_sfwd[i], exer_time[i], 1.0,
                              SRT_CALL, SRT_NORMAL, &std);
      }
      std *= sqrt(exer_time[i]);

      /*	Strike */
      if ((!short_strike) || (!strike_type)) {
        ex_sstrike[i] = ex_sfwd[i];
      } else if (strike_type == 1) {
        ex_sstrike[i] = short_strike[i];
      } else if (strike_type == 2) {
        if (err = swp_f_ForwardRate(
                float_pay_dates[ex_coupon[i]],
                float_pay_dates[ex_coupon[i] + ex_sncpn[i] - 1], swaption_freq,
                swaption_basis, yc_name, ref_rate_name, &swp_rte)) {
          goto FREE_RETURN;
        }

        spr = swp_rte - ex_sfwd[i];

        ex_sstrike[i] = short_strike[i] - spr;
      } else if (strike_type == 3) {
        ex_sstrike[i] = ex_sfwd[i] + short_strike[i] * std;
      }

      /*	Apply max std */
      if (ex_sstrike[i] > ex_sfwd[i] + max_std_short * std) {
        //				ex_sstrike[i] = ex_sfwd[i] + max_std_short
        //* std;
      } else if (ex_sstrike[i] < ex_sfwd[i] - max_std_short * std) {
        //				ex_sstrike[i] = ex_sfwd[i] - max_std_short
        //* std;
      }

      /*	Make sure strikes are positive (actually more than 1bp)
                      otherwise use ATM	*/
      if (ex_sstrike[i] < 1.0e-04) {
        ex_sstrike[i] = ex_sfwd[i];
      }

      err =
          get_cash_vol(vol_curve_name,
                       add_unit(exer_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
                       float_pay_dates[ex_coupon[i] + ex_sncpn[i] - 1],
                       ex_sstrike[i], 0, ref_rate_name, &(ex_svol[i]), &power);
      if (err) {
        goto FREE_RETURN;
      }
      ex_svol[i] += (shift_type == 1 ? ex_svol[i] * vol_shift : vol_shift);

      if (power > 0.5) {
        ex_sprice[i] =
            srt_f_optblksch(ex_sfwd[i], ex_sstrike[i], ex_svol[i], exer_time[i],
                            ex_slvl[i], SRT_PUT, PREMIUM);
      } else {
        ex_sprice[i] =
            srt_f_optblknrm(ex_sfwd[i], ex_sstrike[i], ex_svol[i], exer_time[i],
                            ex_slvl[i], SRT_PUT, PREMIUM);
      }

      cap_price += ex_sprice[i];
    }
  }

  /*	The 1F equivalent case */
  if (one2F == 2 && fix_lambda && one_f_equi) {
    cap_price = 0.0;
    for (i = 0; i < nbexercise; i++) {
      if (i < nbexercise - 1) {
        ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
      } else {
        ex_sncpn[i] = nbcoupon - ex_coupon[i];
      }

      if (ex_sncpn[i] < 2) {
        err = "One exercise date controls less than 2 coupons in "
              "cpd_calib_diagonal";
        goto FREE_RETURN;
      }

      lvl = 0.0;
      for (k = ex_coupon[i] + 1; k < ex_coupon[i] + ex_sncpn[i]; k++) {
        lvl += coupon_cvg[k] * coupon_df[k];
      }
      dfi = swp_f_df(today, float_pay_dates[ex_coupon[i]], yc_name);
      dff = swp_f_df(today, float_pay_dates[ex_coupon[i] + ex_sncpn[i] - 1],
                     yc_name);

      ex_slvl[i] = lvl;
      ex_sfwd[i] = (dfi - dff) / lvl;
      ex_sstrike[i] = ex_sfwd[i];
      //			ex_sstrike[i] = long_strike[i];
    }

    err = ZCMidat_lgmcalibzetalambda(
        nbcoupon, coupon, coupon_time, coupon_df, coupon_cvg, nbexercise,
        exer_time, ex_coupon, ex_sncpn, ex_lstrike, ex_lprice, ex_sstrike,
        ex_sweight, 0.0, ex_zeta, exerFixCoupon, 1, lambda, 1, 0.0, 0.0, 0.0,
        skip_last);

    if (err) {
      goto FREE_RETURN;
    }

    export_lgmsetupG(*lambda, nbcoupon, coupon_time, cpn_G, nbexercise,
                     exer_time, ex_G);

    cap_price =
        lgmcapval1F(nbcoupon, coupon_df, coupon_cvg, cpn_G, nbexercise,
                    ex_coupon, ex_sncpn, ex_sweight, ex_zeta, ex_G, ex_sstrike);

    fix_lambda = 0;
  }

  /*	2.)	Calibrate lambda and zeta */

  err = ZCMidat_lgmcalibzetalambda(
      nbcoupon, coupon, coupon_time, coupon_df, coupon_cvg, nbexercise,
      exer_time, ex_coupon, ex_sncpn, ex_lstrike, ex_lprice, ex_sstrike,
      ex_sweight, cap_price, ex_zeta, exerFixCoupon, fix_lambda, lambda, one2F,
      alpha, gamma, rho, skip_last);

  if (err) {
    goto FREE_RETURN;
  }

  /*	3.)	Transform into sigma */

  *num_sig = nbexercise;
  *sig_time = (double *)calloc(nbexercise, sizeof(double));
  *sig = (double *)calloc(nbexercise, sizeof(double));

  if (!sig_time || !sig) {
    err = "Allocation error (3) in cpd_calib_diagonal";
    goto FREE_RETURN;
  }

  (*sig_time)[0] = exer_time[0];
  (*sig)[0] = sqrt(ex_zeta[0] * 2 * (*lambda) /
                   (exp(2 * (*lambda) * exer_time[0]) - 1.0));

  for (i = 1; i < nbexercise; i++) {
    (*sig_time)[i] = exer_time[i];
    if (ex_zeta[i] > ex_zeta[i - 1]) {
      (*sig)[i] = sqrt((ex_zeta[i] - ex_zeta[i - 1]) * 2 * (*lambda) /
                       (exp(2 * (*lambda) * exer_time[i]) -
                        exp(2 * (*lambda) * exer_time[i - 1])));
    } else {
      smessage("Diagonal calibration failed at exercise year %.2f - "
               "Calibration stopped",
               exer_time[i]);
      for (j = i; j < nbexercise; j++) {
        (*sig)[j] = (*sig)[i - 1];
      }
      i = nbexercise;
    }
  }

  /*	4.)	Save instrument data if required */
  if (inst_data) {
    inst_data->num_inst = nbexercise;
    inst_data->start_dates = (long *)calloc(nbexercise, sizeof(long));
    inst_data->end_dates = (long *)calloc(nbexercise, sizeof(long));
    inst_data->short_strikes = (double *)calloc(nbexercise, sizeof(double));
    inst_data->long_strikes = (double *)calloc(nbexercise, sizeof(double));

    if (!inst_data->start_dates || !inst_data->end_dates ||
        !inst_data->short_strikes || !inst_data->long_strikes) {
      err = "Allocation error (4) in cpd_calib_diagonal";
      goto FREE_RETURN;
    }

    for (i = 0; i < nbexercise; i++) {
      inst_data->start_dates[i] = float_pay_dates[ex_coupon[i]];
      inst_data->end_dates[i] = float_pay_dates[nbcoupon - 1];
      if (!fix_lambda) {
        inst_data->short_strikes[i] = ex_sstrike[i];
      } else {
        inst_data->short_strikes[i] = 0.0;
      }
      inst_data->long_strikes[i] = ex_lstrike[i];
    }
  }

FREE_RETURN:

  if (float_fixing_dates)
    free(float_fixing_dates);
  if (float_start_dates)
    free(float_start_dates);
  if (float_end_dates)
    free(float_end_dates);
  if (float_pay_dates)
    free(float_pay_dates);
  if (float_cvgs)
    free(float_cvgs);
  if (float_spreads)
    free(float_spreads);
  if (float_pay_times)
    free_dvector(float_pay_times, 0, float_nb_pay_dates - 1);

  if (fix_start_dates)
    free(fix_start_dates);
  if (fix_end_dates)
    free(fix_end_dates);
  if (fix_pay_dates)
    free(fix_pay_dates);
  if (fix_cvgs)
    free(fix_cvgs);
  if (fix_pay_times)
    free_dvector(fix_pay_times, 0, fix_nb_pay_dates - 1);
  if (exerFixCoupon)
    free(exerFixCoupon);

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
