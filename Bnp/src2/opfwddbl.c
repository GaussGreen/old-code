/*******************************************************************************
**
**	srt_f_optfwddbl.c
**
*******************************************************************************/

/* ==========================================================================
   include files
   ========================================================================== */

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/* ==========================================================================
   part_bar_fn(...)
        - Gaussian Method
   ========================================================================== */

static double probafunc(double T1, double T2, double sigma01, double sigma12,
                        double mu01, double mu12, double l, double u, double k,
                        int nbterms);

/* ==========================================================================
   srt_f_optfwddbl(...)
        - Gaussian Method for a partial double barrier
   ========================================================================== */

double srt_f_optfwddbl(double fwd1, double fwd2, double spot, double strike,
                       double b_do, double b_up, double vol1, double vol2,
                       double mat1, double mat2, double disc,
                       SrtCallPutType call_put, int nb_term,
                       SrtGreekType greek) {
  double up, dwn, k, mu1, mu2, bound_dwn, bound_up, premium, shift, answer;

  double drift1, drift2, price_k, price_s, sig1, sig2, sig1_2, sig2_2;

  int cp;

  if (mat1 <= 0) {
    if ((spot <= b_do) || (spot >= b_up)) {
      return 0.0;
    }

    premium = srt_f_optdblbar(spot, fwd2, strike, b_do, b_up, vol2, mat2, disc,
                              call_put, nb_term, greek);
  } else {
    cp = (call_put == SRT_CALL) ? 1 : -1;

    if ((cp == 1) && (strike >= b_up))
      return 0.0;
    if ((cp == -1) && (strike <= b_do))
      return 0.0;

    up = log(b_up / spot);
    dwn = log(b_do / spot);
    k = log(strike / spot);

    if (strike >= b_up)
      k = up;
    if (strike <= b_do)
      k = dwn;

    drift1 = log(fwd1 / spot) / mat1;
    sig1 = vol1;

    if (mat1 < mat2) {
      drift2 = log(fwd2 / fwd1) / (mat2 - mat1);
      sig2 = sqrt((vol2 * vol2 * mat2 - vol1 * vol1 * mat1) / (mat2 - mat1));
    } else {
      return 0.0;
    }

    sig1_2 = sig1 * sig1;
    sig2_2 = sig2 * sig2;

    mu1 = drift1 - sig1_2 / 2;
    mu2 = drift2 - sig2_2 / 2;
    premium = 0.0;

    if (call_put == SRT_CALL) {
      bound_up = up; /*  call  */
      bound_dwn = dwn;
    } else if (call_put == SRT_PUT) {
      bound_up = (-dwn); /*  put  */
      bound_dwn = (-up);
      k = (-k);
      mu1 = (-mu1);
      mu2 = (-mu2);
      sig1_2 = (-sig1_2);
      sig2_2 = (-sig2_2);
    }

    price_s = probafunc(mat1, mat2, sig1, sig2, mu1 + sig1_2, mu2 + sig2_2,
                        bound_dwn, bound_up, k, nb_term);
    /*			= part_dbar_fn( k   ,
                            bound_up  ,
                            bound_dwn  ,
                            mu1 + sig1_2  ,
                            mu2 + sig2_2   ,
                            sig1  ,
                            sig2  ,
                            mat1  ,
                            mat2  ,
                            nb_term ) ;     */

    /*	price_k = part_dbar_fn( k   ,
                            bound_up  ,
                            bound_dwn  ,
                            mu1  ,
                            mu2  ,
                            sig1  ,
                            sig2  ,
                            mat1  ,
                            mat2  ,
                            nb_term ) ; */

    price_k = probafunc(mat1, mat2, sig1, sig2, mu1, mu2, bound_dwn, bound_up,
                        k, nb_term);

    premium = fwd2 * price_s - strike * price_k;
    premium *= disc;
    premium *= cp;
  }

  switch (greek) {
  case PREMIUM: /*** PREMIUM ***/
    answer = premium;
    break;

  case DELTA_FWD1: /*** DELTA FWD at T1 ***/
    shift = fwd1 / 10000;
    answer =
        (srt_f_optfwddbl(fwd1 + shift, fwd2, spot, strike, b_do, b_up, vol1,
                         vol2, mat1, mat2, disc, call_put, nb_term, PREMIUM) -
         premium) /
        shift;
    break;

  case DELTA_FWD2: /*** DELTA FWD at T2 ***/
    shift = fwd2 / 10000;
    answer =
        (srt_f_optfwddbl(fwd1, fwd2 + shift, spot, strike, b_do, b_up, vol1,
                         vol2, mat1, mat2, disc, call_put, nb_term, PREMIUM) -
         premium) /
        shift;
    break;

  case DELTA: /*** DELTA SPOT + FWD1 + FWD2 ***/
    shift = spot / 10000;
    answer =
        (srt_f_optfwddbl(fwd1 * (1 + shift / spot), fwd2 * (1 + shift / spot),
                         spot + shift, strike, b_do, b_up, vol1, vol2, mat1,
                         mat2, disc, call_put, nb_term, PREMIUM) -
         premium) /
        shift;
    break;

  case GAMMA: /*** GAMMA ***/
    shift = spot / 10000;
    answer =
        srt_f_optfwddbl(fwd1 * (1 + shift / spot), fwd2 * (1 + shift / spot),
                        spot + shift, strike, b_do, b_up, vol1, vol2, mat1,
                        mat2, disc, call_put, nb_term, PREMIUM);
    answer =
        srt_f_optfwddbl(fwd1 * (1 - shift / spot), fwd2 * (1 - shift / spot),
                        spot - shift, strike, b_do, b_up, vol1, vol2, mat1,
                        mat2, disc, call_put, nb_term, PREMIUM);
    answer -= 2 * premium;
    answer /= shift * shift;
    break;

  case VEGA1: /*** VEGA at T1***/
    shift = GVOPT.vol_add;
    answer =
        (srt_f_optfwddbl(fwd1, fwd2, spot, strike, b_do, b_up, vol1 + shift,
                         vol2, mat1, mat2, disc, call_put, nb_term, PREMIUM) -
         premium) /
        shift;
    break;

  case VEGA2: /*** VEGA at T2***/
    shift = GVOPT.vol_add;
    answer = (srt_f_optfwddbl(fwd1, fwd2, spot, strike, b_do, b_up, vol1,
                              vol2 + shift, mat1, mat2, disc, call_put, nb_term,
                              PREMIUM) -
              premium) /
             shift;
    break;

  case THETA: /*** THETA  ***/
    shift = YEARS_IN_DAY;
    answer = srt_f_optfwddbl(fwd1, fwd2, spot, strike, b_do, b_up, vol1, vol2,
                             mat1 - shift, mat2 - shift,
                             disc * exp(-shift * log(disc) / mat2), call_put,
                             nb_term, PREMIUM) -
             premium;
    break;

  default:
    answer = UNKNOWN_GREEK;
    break;
  }

  return (answer);
}

double probafunc(double T1, double T2, double sigma01, double sigma12,
                 double mu01, double mu12, double l, double u, double k,
                 int nbterms) {

  double return1;
  double a1, b1, c1, d1, e1, f1, g1, h1;
  double xn;
  double t12 = T2 - T1;
  double one, two;
  int p1;

  return1 = 0.0;

  for (p1 = -nbterms; p1 <= nbterms; p1++) {
    xn = 2 * (u - l) * p1;
    a1 = exp(mu12 * xn / (sigma12 * sigma12)) *
         srt_f_intnrm_i2(xn - k + u + mu12 * t12, sigma12 * sqrt(t12),
                         u - mu01 * T1, sigma01 * sqrt(T1));
    b1 = exp(mu12 * xn / (sigma12 * sigma12)) *
         srt_f_intnrm_i2(xn + mu12 * t12, sigma12 * sqrt(t12), u - mu01 * T1,
                         sigma01 * sqrt(T1));
    c1 = exp(-mu12 * xn / (sigma12 * sigma12)) *
         srt_f_intnrm_j1(2 * mu12 / (sigma12 * sigma12),
                         -xn - k + u + mu12 * t12, sigma12 * sqrt(t12),
                         u - mu01 * T1, sigma01 * sqrt(T1));
    d1 =
        exp(-mu12 * xn / (sigma12 * sigma12)) *
        srt_f_intnrm_j1(2 * mu12 / (sigma12 * sigma12), -xn + mu12 * t12,
                        sigma12 * sqrt(t12), u - mu01 * T1, sigma01 * sqrt(T1));

    e1 = exp(mu12 * xn / (sigma12 * sigma12)) *
         srt_f_intnrm_i2(xn - k + l + mu12 * t12, sigma12 * sqrt(t12),
                         l - mu01 * T1, sigma01 * sqrt(T1));
    f1 = exp(mu12 * xn / (sigma12 * sigma12)) *
         srt_f_intnrm_i2(xn + l - u + mu12 * t12, sigma12 * sqrt(t12),
                         l - mu01 * T1, sigma01 * sqrt(T1));
    g1 = exp((2 * u - 2 * l - xn) * mu12 / (sigma12 * sigma12)) *
         srt_f_intnrm_j1(2 * mu12 / (sigma12 * sigma12),
                         2 * u - xn - k - l + mu12 * t12, sigma12 * sqrt(t12),
                         l - mu01 * T1, sigma01 * sqrt(T1));
    h1 = exp((2 * u - 2 * l - xn) * mu12 / (sigma12 * sigma12)) *
         srt_f_intnrm_j1(2 * mu12 / (sigma12 * sigma12),
                         -xn + u - l + mu12 * t12, sigma12 * sqrt(t12),
                         l - mu01 * T1, sigma01 * sqrt(T1));
    one = a1 - b1 - c1 + d1;
    two = -e1 + f1 + g1 - h1;

    return1 += one + two;
  }

  return return1;
}

/* Old mathematica version */
/*
double probafunc(double T1  , double T2  , double sigma01  ,double sigma12
,double mu01  ,double mu12   , double d   , double u  , double k   , int
nbterms)

        {
        double return1;
        double a;
        int p1;
        return1=0.0;
        for (p1=- nbterms ; p1<=nbterms;p1++)

        {
        double powsig12_2=Power(sigma12  ,2);
        double powsig01_2=Power(sigma01  ,2);
        double powsig12_2t1=powsig12_2*T1;
        double powsig01_2t2=powsig01_2*T2;
        double powsig12_2t2=powsig12_2*T2;
        double powsig01_2t1=powsig01_2*T1;
        double mu01t1=mu01*T1;
        double mu12t2=mu12*T2;
        double powsig12_4=powsig12_2*powsig12_2;
        double powsig01_4=powsig01_2*powsig01_2;


                a=
        -(Power(E  ,d*p1 - 2*d*mu12*p1/powsig12_2 - p1*u +
2*mu12*p1*u/powsig12_2)* CBND((-2*d + 2*mu01t1 -
powsig01_2t1)/(2*sigma01*Sqrt(T1))  , (2*k + 4*d*p1 - 2*mu01t1 + 2*mu12*T1 +
powsig01_2t1 - powsig12_2t1 - 2*mu12t2 + powsig12_2t2 - 4*p1*u)/
    (2*Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2)))) +
        Power(E  ,d*p1 - 2*d*mu12*p1/powsig12_2 - p1*u +
2*mu12*p1*u/powsig12_2)* CBND((-2*d + 2*mu01t1 -
powsig01_2t1)/(2*sigma01*Sqrt(T1))  , (4*d*p1 - 2*mu01t1 + 2*mu12*T1 +
powsig01_2t1 - powsig12_2t1 - 2*mu12t2 + powsig12_2t2 + 2*u - 4*p1*u)/
    (2*Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2))) -
        Power(E  ,-(d*p1) + 2*d*mu12*p1/powsig12_2 + mu01t1 + 2*Power(mu12
,2)*powsig01_2t1/powsig12_4 - 2*mu01*mu12*T1/powsig12_2 -
mu12*powsig01_2t1/powsig12_2 - u + p1*u + 2*mu12*u/powsig12_2 -
    2*mu12*p1*u/powsig12_2)*CBND((2*d*powsig12_2 + 4*mu12*powsig01_2t1 -
2*mu01*powsig12_2t1 - powsig01_2*powsig12_2t1)/(2*sigma01*powsig12_2*Sqrt(T1)) ,
    (2*k*powsig12_2 - 4*d*p1*powsig12_2 - 4*mu12*powsig01_2t1 +
2*mu01*powsig12_2t1 + 2*mu12*powsig12_2t1 + powsig01_2*powsig12_2t1 -
powsig12_4*T1 - 2*mu12*powsig12_2t2 + powsig12_4*T2 - 4*powsig12_2*u +
4*p1*powsig12_2*u)/ (2*powsig12_2*Sqrt(powsig01_2t1 - powsig12_2t1 +
powsig12_2t2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2))) +
        Power(E  ,-(d*p1) + 2*d*mu12*p1/powsig12_2 + mu01t1 + 2*Power(mu12
,2)*powsig01_2t1/powsig12_4 - 2*mu01*mu12*T1/powsig12_2 -
mu12*powsig01_2t1/powsig12_2 - u + p1*u + 2*mu12*u/powsig12_2 -
    2*mu12*p1*u/powsig12_2)*CBND((2*d*powsig12_2 + 4*mu12*powsig01_2t1 -
2*mu01*powsig12_2t1 - powsig01_2*powsig12_2t1)/(2*sigma01*powsig12_2*Sqrt(T1)) ,
    (-4*d*p1*powsig12_2 - 4*mu12*powsig01_2t1 + 2*mu01*powsig12_2t1 +
2*mu12*powsig12_2t1 + powsig01_2*powsig12_2t1 - powsig12_4*T1 -
2*mu12*powsig12_2t2 + powsig12_4*T2 - 2*powsig12_2*u + 4*p1*powsig12_2*u)/
    (2*powsig12_2*Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2))) +
        Power(E  ,d*p1 - 2*d*mu12*p1/powsig12_2 - p1*u +
2*mu12*p1*u/powsig12_2)* CBND((2*mu01t1 - powsig01_2t1 -
2*u)/(2*sigma01*Sqrt(T1))  , (2*k + 4*d*p1 - 2*mu01t1 + 2*mu12*T1 + powsig01_2t1
- powsig12_2t1 - 2*mu12t2 + powsig12_2t2 - 4*p1*u)/ (2*Sqrt(powsig01_2t1 -
powsig12_2t1 + powsig12_2t2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2))) -
        Power(E  ,d*p1 - 2*d*mu12*p1/powsig12_2 - p1*u +
2*mu12*p1*u/powsig12_2)* CBND((2*mu01t1 - powsig01_2t1 -
2*u)/(2*sigma01*Sqrt(T1))  , (4*d*p1 - 2*mu01t1 + 2*mu12*T1 + powsig01_2t1 -
powsig12_2t1 - 2*mu12t2 + powsig12_2t2 + 2*u - 4*p1*u)/ (2*Sqrt(powsig01_2t1 -
powsig12_2t1 + powsig12_2t2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2))) +
        Power(E  ,-(d*p1) + 2*d*mu12*p1/powsig12_2 + mu01t1 + 2*Power(mu12
,2)*powsig01_2t1/powsig12_4 - 2*mu01*mu12*T1/powsig12_2 -
mu12*powsig01_2t1/powsig12_2 - u + p1*u + 2*mu12*u/powsig12_2 -
    2*mu12*p1*u/powsig12_2)*CBND((4*mu12*powsig01_2t1 - 2*mu01*powsig12_2t1 -
    powsig01_2*powsig12_2t1 + 2*powsig12_2*u)/(2*sigma01*powsig12_2*Sqrt(T1))  ,
    (2*k*powsig12_2 - 4*d*p1*powsig12_2 - 4*mu12*powsig01_2t1 +
2*mu01*powsig12_2t1 + 2*mu12*powsig12_2t1 + powsig01_2*powsig12_2t1 -
powsig12_4*T1 - 2*mu12*powsig12_2t2 + powsig12_4*T2 - 4*powsig12_2*u +
4*p1*powsig12_2*u)/ (2*powsig12_2*Sqrt(powsig01_2t1 - powsig12_2t1 +
powsig12_2t2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2))) -
        Power(E  ,-(d*p1) + 2*d*mu12*p1/powsig12_2 + mu01t1 + 2*Power(mu12
,2)*powsig01_2t1/powsig12_4 - 2*mu01*mu12*T1/powsig12_2 -
mu12*powsig01_2t1/powsig12_2 - u + p1*u + 2*mu12*u/powsig12_2 -
    2*mu12*p1*u/powsig12_2)*CBND((4*mu12*powsig01_2t1 - 2*mu01*powsig12_2t1 -
    powsig01_2*powsig12_2t1 + 2*powsig12_2*u)/(2*sigma01*powsig12_2*Sqrt(T1))  ,
    (-4*d*p1*powsig12_2 - 4*mu12*powsig01_2t1 + 2*mu01*powsig12_2t1 +
2*mu12*powsig12_2t1 + powsig01_2*powsig12_2t1 - powsig12_4*T1 -
2*mu12*powsig12_2t2 + powsig12_4*T2 - 2*powsig12_2*u + 4*p1*powsig12_2*u)/
    (2*powsig12_2*Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(powsig01_2t1 - powsig12_2t1 + powsig12_2t2)));

        return1+=a;



/* this works   , in case of bugs use it again

        a=
        -(Power(E  ,d*p1 - 2*d*mu12*p1/Power(sigma12  ,2) - p1*u +
2*mu12*p1*u/Power(sigma12  ,2))* CBND((-2*d + 2*mu01*T1 - Power(sigma01
,2)*T1)/(2*sigma01*Sqrt(T1))  , (2*k + 4*d*p1 - 2*mu01*T1 + 2*mu12*T1 +
Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 - 2*mu12*T2 + Power(sigma12 ,2)*T2
- 4*p1*u)/ (2*Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 + Power(sigma12
,2)*T2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 +
Power(sigma12  ,2)*T2)))) + Power(E  ,d*p1 - 2*d*mu12*p1/Power(sigma12  ,2) -
p1*u + 2*mu12*p1*u/Power(sigma12  ,2))* CBND((-2*d + 2*mu01*T1 - Power(sigma01
,2)*T1)/(2*sigma01*Sqrt(T1))  , (4*d*p1 - 2*mu01*T1 + 2*mu12*T1 + Power(sigma01
,2)*T1 - Power(sigma12  ,2)*T1 - 2*mu12*T2 + Power(sigma12  ,2)*T2 + 2*u -
4*p1*u)/ (2*Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 + Power(sigma12
,2)*T2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 +
Power(sigma12  ,2)*T2))) - Power(E  ,-(d*p1) + 2*d*mu12*p1/Power(sigma12  ,2) +
mu01*T1 + 2*Power(mu12  ,2)*Power(sigma01  ,2)*T1/Power(sigma12  ,4) -
    2*mu01*mu12*T1/Power(sigma12  ,2) - mu12*Power(sigma01  ,2)*T1/Power(sigma12
,2) - u + p1*u + 2*mu12*u/Power(sigma12  ,2) - 2*mu12*p1*u/Power(sigma12
,2))*CBND((2*d*Power(sigma12  ,2) + 4*mu12*Power(sigma01  ,2)*T1 -
2*mu01*Power(sigma12  ,2)*T1 - Power(sigma01  ,2)*Power(sigma12
,2)*T1)/(2*sigma01*Power(sigma12  ,2)*Sqrt(T1))  , (2*k*Power(sigma12  ,2) -
4*d*p1*Power(sigma12  ,2) - 4*mu12*Power(sigma01  ,2)*T1 + 2*mu01*Power(sigma12
,2)*T1 + 2*mu12*Power(sigma12  ,2)*T1 + Power(sigma01  ,2)*Power(sigma12  ,2)*T1
- Power(sigma12  ,4)*T1 - 2*mu12*Power(sigma12  ,2)*T2 + Power(sigma12  ,4)*T2 -
4*Power(sigma12  ,2)*u + 4*p1*Power(sigma12  ,2)*u)/ (2*Power(sigma12
,2)*Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 + Power(sigma12  ,2)*T2))
,
    -(sigma01*Sqrt(T1)/Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 +
Power(sigma12  ,2)*T2))) + Power(E  ,-(d*p1) + 2*d*mu12*p1/Power(sigma12  ,2) +
mu01*T1 + 2*Power(mu12  ,2)*Power(sigma01  ,2)*T1/Power(sigma12  ,4) -
    2*mu01*mu12*T1/Power(sigma12  ,2) - mu12*Power(sigma01  ,2)*T1/Power(sigma12
,2) - u + p1*u + 2*mu12*u/Power(sigma12  ,2) - 2*mu12*p1*u/Power(sigma12
,2))*CBND((2*d*Power(sigma12  ,2) + 4*mu12*Power(sigma01  ,2)*T1 -
2*mu01*Power(sigma12  ,2)*T1 - Power(sigma01  ,2)*Power(sigma12
,2)*T1)/(2*sigma01*Power(sigma12  ,2)*Sqrt(T1))  ,
    (-4*d*p1*Power(sigma12  ,2) - 4*mu12*Power(sigma01  ,2)*T1 +
2*mu01*Power(sigma12  ,2)*T1 + 2*mu12*Power(sigma12  ,2)*T1 + Power(sigma01
,2)*Power(sigma12  ,2)*T1 - Power(sigma12  ,4)*T1 - 2*mu12*Power(sigma12  ,2)*T2
+ Power(sigma12  ,4)*T2 - 2*Power(sigma12  ,2)*u + 4*p1*Power(sigma12  ,2)*u)/
    (2*Power(sigma12  ,2)*Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 +
Power(sigma12  ,2)*T2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 +
Power(sigma12  ,2)*T2))) + Power(E  ,d*p1 - 2*d*mu12*p1/Power(sigma12  ,2) -
p1*u + 2*mu12*p1*u/Power(sigma12  ,2))* CBND((2*mu01*T1 - Power(sigma01  ,2)*T1
- 2*u)/(2*sigma01*Sqrt(T1))  , (2*k + 4*d*p1 - 2*mu01*T1 + 2*mu12*T1 +
Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 - 2*mu12*T2 + Power(sigma12 ,2)*T2
- 4*p1*u)/ (2*Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 + Power(sigma12
,2)*T2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 +
Power(sigma12  ,2)*T2))) - Power(E  ,d*p1 - 2*d*mu12*p1/Power(sigma12  ,2) -
p1*u + 2*mu12*p1*u/Power(sigma12  ,2))* CBND((2*mu01*T1 - Power(sigma01  ,2)*T1
- 2*u)/(2*sigma01*Sqrt(T1))  , (4*d*p1 - 2*mu01*T1 + 2*mu12*T1 + Power(sigma01
,2)*T1 - Power(sigma12  ,2)*T1 - 2*mu12*T2 + Power(sigma12  ,2)*T2 + 2*u -
4*p1*u)/ (2*Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 + Power(sigma12
,2)*T2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 +
Power(sigma12  ,2)*T2))) + Power(E  ,-(d*p1) + 2*d*mu12*p1/Power(sigma12  ,2) +
mu01*T1 + 2*Power(mu12  ,2)*Power(sigma01  ,2)*T1/Power(sigma12  ,4) -
    2*mu01*mu12*T1/Power(sigma12  ,2) - mu12*Power(sigma01  ,2)*T1/Power(sigma12
,2) - u + p1*u + 2*mu12*u/Power(sigma12  ,2) - 2*mu12*p1*u/Power(sigma12
,2))*CBND((4*mu12*Power(sigma01  ,2)*T1 - 2*mu01*Power(sigma12  ,2)*T1 -
    Power(sigma01  ,2)*Power(sigma12  ,2)*T1 + 2*Power(sigma12
,2)*u)/(2*sigma01*Power(sigma12  ,2)*Sqrt(T1))  , (2*k*Power(sigma12  ,2) -
4*d*p1*Power(sigma12  ,2) - 4*mu12*Power(sigma01  ,2)*T1 + 2*mu01*Power(sigma12
,2)*T1 + 2*mu12*Power(sigma12  ,2)*T1 + Power(sigma01  ,2)*Power(sigma12  ,2)*T1
- Power(sigma12  ,4)*T1 - 2*mu12*Power(sigma12  ,2)*T2 + Power(sigma12  ,4)*T2 -
4*Power(sigma12  ,2)*u + 4*p1*Power(sigma12  ,2)*u)/ (2*Power(sigma12
,2)*Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 + Power(sigma12  ,2)*T2))
,
    -(sigma01*Sqrt(T1)/Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 +
Power(sigma12  ,2)*T2))) - Power(E  ,-(d*p1) + 2*d*mu12*p1/Power(sigma12  ,2) +
mu01*T1 + 2*Power(mu12  ,2)*Power(sigma01  ,2)*T1/Power(sigma12  ,4) -
    2*mu01*mu12*T1/Power(sigma12  ,2) - mu12*Power(sigma01  ,2)*T1/Power(sigma12
,2) - u + p1*u + 2*mu12*u/Power(sigma12  ,2) - 2*mu12*p1*u/Power(sigma12
,2))*CBND((4*mu12*Power(sigma01  ,2)*T1 - 2*mu01*Power(sigma12  ,2)*T1 -
    Power(sigma01  ,2)*Power(sigma12  ,2)*T1 + 2*Power(sigma12
,2)*u)/(2*sigma01*Power(sigma12  ,2)*Sqrt(T1))  ,
    (-4*d*p1*Power(sigma12  ,2) - 4*mu12*Power(sigma01  ,2)*T1 +
2*mu01*Power(sigma12  ,2)*T1 + 2*mu12*Power(sigma12  ,2)*T1 + Power(sigma01
,2)*Power(sigma12  ,2)*T1 - Power(sigma12  ,4)*T1 - 2*mu12*Power(sigma12  ,2)*T2
+ Power(sigma12  ,4)*T2 - 2*Power(sigma12  ,2)*u + 4*p1*Power(sigma12  ,2)*u)/
    (2*Power(sigma12  ,2)*Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 +
Power(sigma12  ,2)*T2))  ,
    -(sigma01*Sqrt(T1)/Sqrt(Power(sigma01  ,2)*T1 - Power(sigma12  ,2)*T1 +
Power(sigma12  ,2)*T2)));

        return1+=a;

  */
/*
        }


return return1;

        }

        */
