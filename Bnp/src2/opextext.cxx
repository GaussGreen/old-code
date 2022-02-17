/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

static double RTSAFE_STRIKE1;
static double RTSAFE_STRIKE2;
static double RTSAFE_BARRIER;
static double RTSAFE_VOL;
static double RTSAFE_MAT;
static double RTSAFE_DISC;
static int RTSAFE_CP1;
static SrtCallPutType RTSAFE_CALLPUT2;
static SrtBarrierType RTSAFE_SCTYPE2;

/******************************************************************************
 *
 * FUNCTION     	: srt_f_extoptcombar
 *
 * PURPOSE      	: EXTINGUISHABLE COMPOUND BARRIER Option
 *
 * DESCRIPTION  	: an option on a barrier option as underlying
 *
 *******************************************************************************/

/******************************************************************************/
static Err funk(double spot, double *f, double *df) {
  /* function to be used with the Newton-Raphson routine - rtsafe - */
  /* return the value of price - f(v) */

  double fwd1;

  fwd1 = spot / RTSAFE_DISC;

  *f = srt_f_optexting(fwd1, spot, RTSAFE_STRIKE2, RTSAFE_BARRIER, RTSAFE_VOL,
                       RTSAFE_MAT, RTSAFE_DISC, RTSAFE_CALLPUT2, RTSAFE_SCTYPE2,
                       PREMIUM) -
       RTSAFE_STRIKE1;

  *df = srt_f_optexting(fwd1, spot, RTSAFE_STRIKE2, RTSAFE_BARRIER, RTSAFE_VOL,
                        RTSAFE_MAT, RTSAFE_DISC, RTSAFE_CALLPUT2,
                        RTSAFE_SCTYPE2, DELTA_FWD);

  *f *= RTSAFE_CP1;
  *df *= RTSAFE_CP1;

  return NULL;
}

/*****************************************************************************/

static double com_bar_fn3(double xd, double xu, double l_b, double mat,
                          double mu, double sig, int cp) {
  double result = 0.0;
  double n1, n2;
  double sigrt;

  sigrt = sig * sqrt(mat);

  n1 = norm((-cp * xd + mu * mat) / sigrt) -
       exp(2.0 * l_b * mu / (sig * sig)) *
           norm((-cp * xd + 2.0 * l_b + mu * mat) / sigrt);

  n2 = norm((-cp * xu + mu * mat) / sigrt) -
       exp(2.0 * l_b * mu / (sig * sig)) *
           norm((-cp * xu + 2.0 * l_b + mu * mat) / sigrt);

  result = n1 - n2;

  return result;
}

/*********************************************************************************/

static double com_bar_fn2(double x1, double x2, double b1, double b2, double t1,
                          double t2, double mu1, double mu2, double sig1,
                          double sig2) {
  double result = 0.0;
  double n1, n2, n3, n4;
  double rtsig2, rtsig1;

  rtsig2 = sig2 * sqrt(t2 - t1);
  rtsig1 = sig1 * sqrt(t1);

  if (x1 >= 100.0 || x2 >= 100.0) {
    return 0.0;
  } else if (x1 <= -100.0 && x2 <= -100) {
    return 1.0;
  } else {
    n1 = srt_f_intnrm_i1(x1 - x2 + mu2 * (t2 - t1), rtsig2, mu1 * t1 - x1,
                         rtsig1);

    n2 = srt_f_intnrm_i1(x1 - x2 + mu2 * (t2 - t1), rtsig2,
                         mu1 * t1 - x1 + 2.0 * b1, rtsig1);

    n3 = srt_f_intnrm_j2(-2.0 * mu2 / (sig2 * sig2),
                         2.0 * b2 - (x1 + x2) + mu2 * (t2 - t1), rtsig2,
                         mu1 * t1 - x1, rtsig1);

    n4 = srt_f_intnrm_j2(-2.0 * mu2 / (sig2 * sig2),
                         2.0 * b2 - (x1 + x2) + mu2 * (t2 - t1), rtsig2,
                         mu1 * t1 - x1 + 2.0 * b1, rtsig1);

    result = n1 - n2 * exp(2.0 * mu1 * b1 / (sig1 * sig1)) -
             exp(2.0 * mu2 * (b2 - x1) / (sig2 * sig2)) *
                 (n3 - exp(2.0 * mu1 * b1 / (sig1 * sig1)) * n4);

    return result;
  }
}

/******************************************************************************/
static double com_bar_fn1(double xd1, double xu1, double xd2, double xu2,
                          double b1, double b2, double t1, double t2,
                          double mu1, double mu2, double sig1, double sig2) {
  double result = 0.0;
  double n1, n2, n3, n4;

  n1 = com_bar_fn2(xd1, xd2, b1, b2, t1, t2, mu1, mu2, sig1, sig2);

  n2 = com_bar_fn2(xu1, xd2, b1, b2, t1, t2, mu1, mu2, sig1, sig2);

  n3 = com_bar_fn2(xd1, xu2, b1, b2, t1, t2, mu1, mu2, sig1, sig2);

  n4 = com_bar_fn2(xu1, xu2, b1, b2, t1, t2, mu1, mu2, sig1, sig2);

  result = n1 - n2 - (n3 - n4);

  return result;
}

/*********************************************************************************/

/*********************************************************************************/

static double z_brak(double x1, double x2, double xacc, int niter, double fwd,
                     double strike1, double strike2, double barrier, double sig,
                     double mat, double disc, int cp1, SrtCallPutType ot2,
                     SrtBarrierType sc_type2)
/* routine from Numerical Recipes        , pg. 354 */
{
  int i, j, jmax;
  double dx, dxold, f, fwd1;
  double rts;

  dxold = fabs(x2 - x1);

  for (i = 0; i <= niter; i++) {
    jmax = (int)pow(2.0, (double)i);
    dx = dxold / ((double)jmax);

    for (j = 0; j < jmax; j++) {
      rts = x2 - (j + .5) * dx;

      fwd1 = rts / disc;

      f = srt_f_optexting(fwd1, rts, strike2, barrier, sig, mat, disc, ot2,
                          sc_type2, PREMIUM) -
          strike1;

      if (f > xacc) {
        return rts;
      }
    }
  }

  return -1.0;
}

/*********************************************************************************/
static double
srt_f_opt_extcombar_pri(double fwd1, double fwd2, double spot, double strike1,
                        double strike2, double barrier1, double barrier2,
                        double sig1, double sig2, double mat1, double mat2,
                        double disc1, double disc2, SrtCallPutType call_put1,
                        SrtCallPutType call_put2, SrtBarrierType scud_type1,
                        SrtBarrierType scud_type2) {
  double l_b1, l_b2;
  double l_k2;
  double l_xstar_u, l_xstar_d;
  double premium;
  double spot_mid;
  double spot_star_up, spot_star_down;
  double xacc, xacc1;
  double disc12;
  double mu1, mu2;
  double rate1, rate2;
  double n1, n2;
  double x1, x2;
  double xd1, xu1, xd2, xu2, bar1, bar2;

  int cp1, cp2;
  int du1, du2;

  Err err;

  /* Accuracy for Newton-Rhapson rtsafe algorithm */
  xacc = 1.0e-7;

  rate1 = log(fwd1 / spot) / mat1;

  if (mat2 > mat1) {
    rate2 = log(fwd2 / fwd1) / (mat2 - mat1);
  } else {
    rate2 = 0.0;
  }

  mu1 = rate1 - (sig1 * sig1 / 2);
  mu2 = rate2 - (sig2 * sig2 / 2);

  disc12 = disc2 / disc1;

  l_b1 = log(barrier1 / spot);
  l_b2 = log(barrier2 / spot);

  l_k2 = log(strike2 / spot);

  cp1 = (call_put1 == SRT_CALL) ? 1 : -1;
  cp2 = (call_put2 == SRT_CALL) ? 1 : -1;

  du1 = (scud_type1 == SRT_DOWN) ? 1 : -1;
  du2 = (scud_type2 == SRT_DOWN) ? 1 : -1;

  if (l_b1 * du1 > 0)
    return 0;

  if (l_b2 * du2 > 0)
    return 0;

  if (scud_type2 == SRT_DOWN) /*** underlying barrier  down ***/
  {
    x1 = barrier2;

    if (call_put2 == SRT_CALL) /*** compound call  ***/
    {
      x2 = DMAX(1.1 * strike2 + 5 * sig2 * sqrt(mat2 - mat1) + barrier2, spot);
    } else if (call_put2 == SRT_PUT) /*** compound put  ***/
    {
      x2 = DMAX(strike2 + 10 * sig2 * sqrt(mat2 - mat1), spot);
    }
  } else if (scud_type2 == SRT_UP) /***  underlying barrier up  ***/
  {
    x1 = DMAX(0.01 * spot,
              DMIN(0.5 * spot, spot - 5 * sig2 * sqrt(mat2 - mat1)));

    if (call_put2 == SRT_CALL) /*** compound call  ***/
    {
      x2 = barrier2;
    } else if (call_put2 == SRT_PUT) /*** compound put  ***/
    {
      x2 = DMAX(strike2 + 10 * sig2 * sqrt(mat2 - mat1), spot);
    }
  }

  RTSAFE_STRIKE1 = strike1;
  RTSAFE_STRIKE2 = strike2;
  RTSAFE_BARRIER = barrier2;
  RTSAFE_VOL = sig2;
  RTSAFE_MAT = mat2 - mat1;
  RTSAFE_DISC = disc12;
  RTSAFE_CP1 = cp1;
  RTSAFE_CALLPUT2 = call_put2;
  RTSAFE_SCTYPE2 = scud_type2;

  if ((scud_type2 == SRT_UP) &&
      (call_put2 == SRT_CALL)) /***  underlying barrier up  ***/
  {
    xacc1 = 10.0 * xacc;

    spot_mid = z_brak(x1, x2, xacc1, 10, fwd2, strike1, strike2, barrier2, sig2,
                      (mat2 - mat1), disc12, cp1, call_put2, scud_type2);

    if (spot_mid == -1.0) {
      return 0.0;
    }

    /*  no value of the spot which gives this precise strike */
    err = rtsafe(funk, x1, spot_mid, xacc, 100, &spot_star_down);
    if (err)
      return INFINITY;

    l_xstar_d = log(spot_star_down / spot);

    x1 = spot_mid;
  }

  err = rtsafe(funk, x1, x2, xacc, 100, &spot_star_up);
  if (err)
    return INFINITY;

  l_xstar_u = log(spot_star_up / spot);

  if (scud_type1 == SRT_DOWN) /*** underlying barrier  down ***/
  {
    if (scud_type2 == SRT_DOWN) /*** underlying barrier  down ***/
    {
      bar1 = l_b1;
      bar2 = l_b2;

      if (call_put1 == SRT_CALL) /*** compound call  ***/
      {
        if (call_put2 == SRT_CALL) /*** underlying call  ***/
        {
          xd1 = DMAX(l_xstar_u, l_b1);
          xu1 = 100;
          xd2 = l_k2;
          xu2 = 100;
        }
      }
    }
  }

  if (scud_type1 == SRT_UP) /*** underlying barrier  up ***/
  {
    if (scud_type2 == SRT_UP) /*** underlying barrier  up ***/
    {
      bar1 = -l_b1;
      bar2 = -l_b2;

      if (call_put1 == SRT_CALL) /*** compound call  ***/
      {
        if (call_put2 == SRT_CALL) /*** underlying call  ***/
        {
          xd1 = -DMIN(l_xstar_u, l_b1);
          xu1 = -l_xstar_d;

          xd2 = -l_b2;
          xu2 = -l_k2;

          mu1 = -mu1;
          mu2 = -mu2;
        }
      }
    }
  }

  n1 =
      com_bar_fn1(xd1, xu1, xd2, xu2, bar1, bar2, mat1, mat2,
                  mu1 + du1 * sig1 * sig1, mu2 + du2 * sig2 * sig2, sig1, sig2);

  premium = cp2 * fwd2 * disc2 * n1;

  n2 = com_bar_fn1(xd1, xu1, xd2, xu2, bar1, bar2, mat1, mat2, mu1, mu2, sig1,
                   sig2);

  premium -= cp2 * strike2 * disc2 * n2;

  n2 = com_bar_fn3(xd1, xu1, bar1, mat1, mu1, sig1, cp1);

  premium -= strike1 * disc1 * n2;
  premium *= cp1;

  return premium;

} /* END srt_f_opt_extcombar_pri ()  */

/********************************************************************************/
double srt_f_opt_extcombar(double fwd1, double fwd2, double spot,
                           double strike1, double strike2, double barrier1,
                           double barrier2, double sig1, double sig2,
                           double mat1, double mat2, double disc1, double disc2,
                           SrtCallPutType call_put1, SrtCallPutType call_put2,
                           SrtBarrierType scud_type1, SrtBarrierType scud_type2,
                           SrtGreekType greek) {
  double fval;
  double x1, x2;
  double f1, f2;
  double delta, gamma, vega, theta;

  if ((greek == PREMIUM) || (greek == GAMMA) || (greek == GAMMA_FWD1) ||
      (greek == GAMMA_FWD2)) {
    fval = srt_f_opt_extcombar_pri(
        fwd1, fwd2, spot, strike1, strike2, barrier1, barrier2, sig1, sig2,
        mat1, mat2, disc1, disc2, call_put1, call_put2, scud_type1, scud_type2);
  }

  if (greek == PREMIUM) { /*  return premium  */
    return fval;
  } else if (greek == DELTA_FWD1) { /*  delta forward1  */
    x1 = fwd1 * 1.01;
    x2 = fwd1 * 0.99;

    f1 = srt_f_opt_extcombar_pri(x1, fwd2, spot, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2);

    f2 = srt_f_opt_extcombar_pri(x2, fwd2, spot, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2);

    delta = (f1 - f2) / (x1 - x2);

    return delta;
  } else if (greek == DELTA_FWD2) { /*  delta fwd2  */
    x1 = fwd2 * 1.01;
    x2 = fwd2 * 0.99;

    delta =
        (srt_f_opt_extcombar_pri(fwd1, x1, spot, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2)

         - srt_f_opt_extcombar_pri(fwd1, x2, spot, strike1, strike2, barrier1,
                                   barrier2, sig1, sig2, mat1, mat2, disc1,
                                   disc2, call_put1, call_put2, scud_type1,
                                   scud_type2)) /
        (x1 - x2);

    return delta;
  } else if (greek == DELTA) /*  delta spot  */
  {
    x1 = spot * 1.01;
    x2 = spot * 0.99;

    delta =
        (srt_f_opt_extcombar_pri(fwd1, fwd2, x1, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2)

         -
         srt_f_opt_extcombar_pri(fwd1, fwd2, x2, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2))

        / (x1 - x2);

    return delta;
  } else if (greek == DELTAX) { /*  delta strike1  */
    x1 = strike1 * 1.01;
    x2 = strike1 * 0.99;

    delta =
        (srt_f_opt_extcombar_pri(fwd1, fwd2, spot, x1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2)

         -
         srt_f_opt_extcombar_pri(fwd1, fwd2, spot, x2, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2))

        / (x1 - x2);

    return delta;
  } else if (greek == DELTAY) { /*  delta strike2  */
    x1 = strike2 * 1.01;
    x2 = strike2 * 0.99;

    delta =
        (srt_f_opt_extcombar_pri(fwd1, fwd2, spot, strike1, x1, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2)

         - srt_f_opt_extcombar_pri(fwd1, fwd2, spot, strike1, x2, barrier1,
                                   barrier2, sig1, sig2, mat1, mat2, disc1,
                                   disc2, call_put1, call_put2, scud_type1,
                                   scud_type2)) /
        (x1 - x2);

    return delta;
  } else if (greek == GAMMA_FWD1) { /*  gamma fwd1  */
    x1 = 1.01 * fwd1;
    x2 = 0.99 * fwd1;

    gamma =
        (srt_f_opt_extcombar_pri(x1, fwd2, spot, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2) -
         2.0 * fval +
         srt_f_opt_extcombar_pri(x2, fwd2, spot, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2))

        / ((x1 - fwd1) * (fwd1 - x2));

    return gamma;
  } else if (greek == GAMMA_FWD2) /*  gamma fwd2   */

  {

    x1 = fwd2 * 1.01;
    x2 = fwd2 * 0.99;

    gamma =
        (srt_f_opt_extcombar_pri(fwd1, x1, spot, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2)

         - 2.0 * fval +
         srt_f_opt_extcombar_pri(fwd1, x2, spot, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2))

        / ((x1 - fwd2) * (fwd2 - x2));

    return gamma;
  } else if (greek == GAMMA) { /*  gamma spot   */
    x1 = spot * 1.01;
    x2 = spot * 0.99;

    gamma =
        (srt_f_opt_extcombar_pri(fwd1, fwd2, x1, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2) -
         2.0 * fval +
         srt_f_opt_extcombar_pri(fwd1, fwd2, x2, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2))

        / ((x1 - spot) * (spot - x2));

    return gamma;
  } else if (greek == VEGA1) /*  vega sig1  */

  {

    x1 = sig1 * 1.01;
    x2 = sig1 * 0.99;

    vega =
        (srt_f_opt_extcombar_pri(fwd1, fwd2, spot, strike1, strike2, barrier1,
                                 barrier2, x1, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2)

         -
         srt_f_opt_extcombar_pri(fwd1, fwd2, spot, strike1, strike2, barrier1,
                                 barrier2, x2, sig2, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2))

        / (x1 - x2);

    return vega;
  } else if (greek == VEGA2) /*  vega sig2  */
  {
    x1 = sig2 * 1.01;
    x2 = sig2 * 0.99;

    vega =
        (srt_f_opt_extcombar_pri(fwd1, fwd2, spot, strike1, strike2, barrier1,
                                 barrier2, sig1, x1, mat1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2)

         - srt_f_opt_extcombar_pri(fwd1, fwd2, spot, strike1, strike2, barrier1,
                                   barrier2, sig1, x2, mat1, mat2, disc1, disc2,
                                   call_put1, call_put2, scud_type1,
                                   scud_type2)) /
        (x1 - x2);

    return vega;
  } else if (greek == THETA) { /*  theta1 = theta mat1  */
    x1 = mat1 * 1.01;
    x2 = mat1 * 0.99;

    theta =
        (srt_f_opt_extcombar_pri(fwd1, fwd2, spot, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, x1, mat2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2)

         - srt_f_opt_extcombar_pri(fwd1, fwd2, spot, strike1, strike2, barrier1,
                                   barrier2, sig1, sig2, x2, mat2, disc1, disc2,
                                   call_put1, call_put2, scud_type1,
                                   scud_type2)) /
        (x1 - x2);

    return theta;
  } else if (greek == THETA_1D) { /*  theta2 = theta mat2  */
    x1 = mat2 * 1.01;
    x2 = mat2 * 0.99;

    theta =
        (srt_f_opt_extcombar_pri(fwd1, fwd2, spot, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, x1, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2)

         -
         srt_f_opt_extcombar_pri(fwd1, fwd2, spot, strike1, strike2, barrier1,
                                 barrier2, sig1, sig2, mat1, x2, disc1, disc2,
                                 call_put1, call_put2, scud_type1, scud_type2))

        / (x1 - x2);

    return theta;
  } else
    return UNKNOWN_GREEK;

} /* END srt_f_opt_extcombar() */
