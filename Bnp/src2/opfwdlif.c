/*******************************************************************************
**
**              Copyright (c) 1993 PARIBAS Capital Markets Group
**
********************************************************************************
**
**      SYSTEM:         SRT     SORT  , Fixed Income 2020 Addins
**      SUB_SYSTEM:     OPT     Option Tools
**
**      MODULE NAME:    SRT_F_OPTFWDLIF
**
**      PURPOSE:        FORWARD LIFE OPTION
**
**      AUTHORS:        Julia Matsumoto
**
**      DATE:           1st June  , 1994
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
**      DATE:           7 September  , 1994
**      VERSION:        2.1
**      REASON:         Documenting  , fixing bugs.
**      DESCRIPTION:    Remy: missing discount factor  , failing in extremal
*situations.
**
********************************************************************************
**
**      AMENDED BY:     Julia Matsumoto
**      DATE:           7 April  , 1995
**      VERSION:        2.2
**      REASON:         Separating mathematical functions into another file.
**      DESCRIPTION:    trapzd_  ,intgrl_smp_ will be in the srt_f_utlnew.c
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
   FUNC: fct_global  (static)
   DESC: Function under numerical integration for the Forward Life Option.
   MODIFIES:
   DECLARATION:
   ======================================================================== */

static double fct_global(double x, va_list argptr) {
  double lambda1, lambda2, spot, barrier, mat1, mat2, vol;
  double ret, m, d, u;

  spot = va_arg(argptr, double);
  barrier = va_arg(argptr, double);
  mat1 = va_arg(argptr, double);
  mat2 = va_arg(argptr, double);
  vol = va_arg(argptr, double);
  lambda1 = va_arg(argptr, double);
  lambda2 = va_arg(argptr, double);

  m = log(barrier) - x;
  d = vol * sqrt(mat1);
  ret = sm_qsimp_list(fct_life, 0.000001,
                      mat2 - mat1 /*life option itself duration*/, m, vol,
                      lambda1);
  ret *= fabs(m) / sqrt(2 * SRT_PI);
  u = (x - log(spot) - lambda2 * mat1) / d;
  ret *= exp(-u * u / 2) / sqrt(2 * SRT_PI) / d;
  return (ret);
}

/* ========================================================================
   FUNC: srt_f_optfwdlife
   DESC:  Forward Start Life Option function.
   MODIFIES:
   DECLARATION:
   ======================================================================== */

double srt_f_optfwdlife(double fwd1, double fwd2, double spot, double barrier,
                        double vol1, double vol2, double mat1, double mat2,
                        double disc1, double disc2, SrtBarrierType down_up,
                        SrtGreekType greek) {
  double myu1, lambda1, lambda2;
  double myu2, myud, d1, d2;
  double answer = 0.0, premium = 0.0;
  double shift;
  double f1, f2, f11, f12, a, intgrl;
  int du;

  du = (down_up == SRT_UP ? 1 : -1);

  if (mat2 == 0.0) {
    premium = 0.0;
  } else if (mat1 == 0.0) {
    premium = srt_f_optlifopt(fwd2, spot, barrier, vol2, mat2, disc2, down_up,
                              PREMIUM);
  } else if (vol1 == 0.0) /* CHECK THIS LOAD OF SHIT */
  {
    if (((fwd2 <= barrier) && (fwd1 < barrier) && (du == 1)) ||
        ((fwd2 >= barrier) && (fwd1 > barrier) && (du == -1)))
      return (0.0);
    else if (((fwd2 <= barrier) && (fwd1 < barrier) && (du == -1)) ||
             ((fwd2 >= barrier) && (fwd1 > barrier) && (du == 1)))
      return (disc1 * (mat2 - mat1));
    else {
      myu1 = log(fwd2 / fwd1) / (mat2 - mat1);
      answer = disc1 * (mat2 - mat1 - log(barrier / fwd1) / myu1);
      return (answer);
    }
  } else {
    myu1 = log(fwd2 / fwd1) / (mat2 - mat1);
    myu2 = log(fwd1 / spot) / mat1;
    myud = myu1 - myu2;

    lambda1 = myu1 - vol1 * vol1 / 2;
    lambda2 = myu2 - vol2 * vol2 / 2;
    d1 = vol1 * sqrt(mat1);
    d2 = vol2 * sqrt(mat2);

    intgrl = sm_qsimp_list(fct_global, (du ? 0 : log(barrier)),
                           (du ? log(barrier) : log(barrier) * 3), spot,
                           barrier, mat1, mat2, vol1, vol2, lambda1, lambda2);

    f11 = bivar(du * (log(barrier / fwd1) + vol1 * vol1 / 2 * mat1) / d1,
                du * (-(log(barrier / fwd2) + vol2 * vol2 / 2 * mat2) / d2),
                -vol1 / vol2 * sqrt(mat1 / mat2));

    f12 = bivar(du * (log(barrier / spot) + (myud + lambda1) * mat1) / d1,
                du * (-log(barrier / spot) - lambda1 * mat1 - myud * mat1) / d2,
                -vol1 / vol2 * sqrt(mat1 / mat2));
    a = pow(barrier / spot, 2 * lambda1 / vol1 / vol2) *
        exp(2 * lambda1 * myud * mat1 / vol1 / vol2);
    f1 = (mat2 - mat1) * (f11 + a * f12);
    f2 =
        (mat2 - mat1) * norm(du * (-log(barrier / spot) + lambda2 * mat1) / d1);

    premium = f1 - intgrl + f2;
    premium *= disc1;
  }

  switch (greek) {
  case PREMIUM:
    answer = premium;
    break;

  case DELTA_FWD1:
    shift = fwd1 / 10000;
    answer = (srt_f_optfwdlife(fwd1 + shift, fwd2, spot, barrier, vol1, vol2,
                               mat1, mat2, disc1, disc2, down_up, PREMIUM) -
              premium) /
             shift;
    break;

  case DELTA_FWD2:
    shift = fwd2 / 10000;
    answer = (srt_f_optfwdlife(fwd1, fwd2 + shift, spot, barrier, vol1, vol2,
                               mat1, mat2, disc1, disc2, down_up, PREMIUM) -
              premium) /
             shift;
    break;

  case DELTA:
    shift = spot / 10000;
    answer =
        (srt_f_optfwdlife(fwd1 * (1 + shift / spot), fwd2 * (1 + shift / spot),
                          spot + shift, barrier, vol1, vol2, mat1, mat2, disc1,
                          disc2, down_up, PREMIUM) -
         premium) /
        shift;
    break;

  case GAMMA:
    shift = spot / 10000;
    answer = srt_f_optfwdlife(
        fwd1 * (1 + shift / spot), fwd2 * (1 + shift / spot), spot + shift,
        barrier, vol1, vol2, mat1, mat2, disc1, disc2, down_up, PREMIUM);
    answer += srt_f_optfwdlife(
        fwd1 * (1 - shift / spot), fwd2 * (1 - shift / spot), spot - shift,
        barrier, vol1, vol2, mat1, mat2, disc1, disc2, down_up, PREMIUM);
    answer -= 2 * premium;
    answer /= shift * shift;
    break;

  case VEGA1:
    shift = GVOPT.vol_add;
    answer = (srt_f_optfwdlife(fwd1, fwd2, spot, barrier, vol1 + shift, vol2,
                               mat1, mat2, disc1, disc2, down_up, PREMIUM) -
              premium) /
             shift;
    break;

  case VEGA2:
    shift = GVOPT.vol_add;
    answer = (srt_f_optfwdlife(fwd1, fwd2, spot, barrier, vol1, vol2 + shift,
                               mat1, mat2, disc1, disc2, down_up, PREMIUM) -
              premium) /
             shift;
    break;

  case THETA:
    shift = YEARS_IN_DAY;
    answer = srt_f_optfwdlife(
                 fwd1, fwd2, spot, barrier, vol1, vol2, mat1 - shift,
                 mat2 - shift, disc1 * exp(-shift * log(disc1) / mat1),
                 disc2 * exp(-shift * log(disc2) / mat2), down_up, PREMIUM) -
             premium;
    break;

  default:
    answer = UNKNOWN_GREEK;
  }

  return (answer);
}
