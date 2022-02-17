/*******************************************************************************
**
**              Copyright (c) 1993 PARIBAS Capital Markets Group
**
********************************************************************************
**
**      SYSTEM:         SRT     SORT        , Fixed Income 2020 Addins
**      SUB_SYSTEM:     OPT     Option Tools
**
**      MODULE NAME:    SRT_F_OPTLIFOPT
**
**      PURPOSE:        LIFE OPTION AND PARTIAL LIFE OPTION
**
**      AUTHORS:        Julia Matsumoto
**
**      DATE:           1st June        , 1994
**
**      VERSION:        2.0
**
**      DESCRIPTION:    ??
**
********************************************************************************
**                      Amendment History
********************************************************************************
**
**      AMENDED BY:     Julia Matsumoto
**      DATE:           7 September        , 1994
**      VERSION:        2.1
**      REASON:         Documenting        , fixing bugs.
**      DESCRIPTION:    Remy: missing discount factor        , failing in
*extremal situations.
**
********************************************************************************
**
**      AMENDED BY:     Julia Matsumoto
**      DATE:           7 April        , 1995
**      VERSION:        2.2
**      REASON:         Separating mathematical functions into another file.
**      DESCRIPTION:
**
********************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/* ========================================================================
   FUNC: fct_life  (static)
   DESC: Function under numerical integration for the Life Option.
   MODIFIES:
   DECLARATION:
   ======================================================================== */

static double fct_life(double x, va_list argptr) {
  double m, lambda, vol, d;
  static double expo, result, u;

  m = va_arg(argptr, double);
  vol = va_arg(argptr, double);
  lambda = va_arg(argptr, double);

  d = vol * sqrt(x);
  u = (m - lambda * x) / d;
  result = exp(-u * u / 2) / d;

  return (result);
}

/* ========================================================================
   FUNC: fct_life2  (static)
   DESC: Function under numerical integration for the Life Option.
   STATUS: Not used
   MODIFIES:
   DECLARATION:
   ======================================================================== */

static double fct_life2(double u, va_list argptr) {
  double m, lambda, vol, d, f;
  static double result;

  m = va_arg(argptr, double);
  vol = va_arg(argptr, double);
  lambda = va_arg(argptr, double);

  f = m - lambda * u / vol / vol;
  d = -f * f / 4 / u;

  return (exp(d) / sqrt(u));
}

/* =========================================================================
   FUNC: srt_f_optlifopt
   DESC: Life Option Pricing function.
   MODIFIES:
   DECLARATION:
   ========================================================================= */

double srt_f_optlifopt(double fwd, double spot, double barrier, double vol,
                       double mat, double disc, SrtBarrierType down_up,
                       SrtGreekType greek) {
  double myu, alpha, lambda, m, d, d_u;
  double premium = 0.0, answer = 0.0;
  double shift;

  d_u = (down_up == SRT_UP) ? 1 : -1;

  if (mat <= 0) {
    premium = 0.0;
  } else if ((spot - barrier) * d_u > 0) {
    premium = disc * mat;
  } else if (vol == 0.0) {
    if ((fwd - barrier) * d_u > 0)
      premium = disc * mat * (1 - log(barrier / spot) / (log(fwd / spot)));
    else
      premium = 0.0;
  } else {
    myu = 1 / mat * log(fwd / spot);
    /*alpha = pow(barrier/spot        ,2*myu/(vol*vol)-1);*/
    lambda = myu - vol * vol / 2;
    m = log(barrier / spot);
    alpha = exp(2 * lambda * m / (vol * vol));
    d = vol * sqrt(mat);

    premium = sm_qsimp_list(fct_life, 0.000001, mat, m, vol, lambda);

    premium *= fabs(m) / sqrt(2 * SRT_PI);
    premium = -premium + mat * (norm(d_u * (-m + lambda * mat) / d) +
                                alpha * norm(d_u * (-m - lambda * mat) / d));

    premium *= disc;
  }

  switch (greek) {
  case PREMIUM:
    answer = premium;
    break;

  case DELTA_FWD:
    shift = fwd / 10000;
    answer = (srt_f_optlifopt(fwd + shift, spot, barrier, vol, mat, disc,
                              down_up, PREMIUM) -
              premium) /
             shift;
    break;

  case DELTA:
    shift = spot / 10000;
    answer = (srt_f_optlifopt(fwd * (1 + shift / spot), spot + shift, barrier,
                              vol, mat, disc, down_up, PREMIUM) -
              premium) /
             shift;
    break;

  case GAMMA:
    shift = spot / 10000;
    answer = (srt_f_optlifopt(fwd * (1 + shift / spot), spot + shift, barrier,
                              vol, mat, disc, down_up, PREMIUM) +
              srt_f_optlifopt(fwd * (1 - shift / spot), spot - shift, barrier,
                              vol, mat, disc, down_up, PREMIUM) -
              2 * premium) /
             (shift * shift);
    break;

  case VEGA:
    shift = GVOPT.vol_add;
    answer = (srt_f_optlifopt(fwd, spot, barrier, vol + shift, mat, disc,
                              down_up, PREMIUM) -
              premium) /
             shift;
    break;

  case THETA:
    shift = YEARS_IN_DAY;
    answer = (srt_f_optlifopt(fwd, spot, barrier, vol, mat - shift,
                              disc * exp(-shift * log(disc) / mat), down_up,
                              PREMIUM) -
              premium);
    break;

  default:
    answer = UNKNOWN_GREEK;
  }

  return (answer);

} /* END srt_f_optlifopt*/

/* ======================================================================== */
