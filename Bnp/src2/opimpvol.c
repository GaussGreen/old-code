/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optimpvol()
 *
 * PURPOSE      	: Calculates the volatility implied from a premium
 *		  	computed using the Black-Scholes formula
 *
 * DESCRIPTION  	: XX
 *
 * CALLS		: srt_f_optblksch()
 *
 * PARAMETERS
 *	INPUT	: premium 	- premium on European option
 *              	: fwd_price   	- forward price of underlying
 *              	: strike        - strike price
 *              	: vol_init      - initial volatility of underlying
 *              	: mat           - maturity  , in years
 *              	: disc        	- discount factor to expiry
 *              	: call_put      - type of option: 0 call  , 1 put
 *               : lognormal_normal  - diffusion process
 *
 * RETURNS      	: vol_guess     - guess of implied volatility
 *
 *******************************************************************************/
#define OPT_EPS 1.0e-15
#define VOL_MULT 5.0
#define VOL_MULT_ITER 3

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

static double optval_accurate(double fwd, double strike, double vol, double mat,
                              double disc_fact, SrtCallPutType call_put,
                              SrtGreekType greek,
                              SrtDiffusionType log_or_norm) {
  double premium;

  if (log_or_norm == SRT_LOGNORMAL) {
    premium = srt_f_optblksch_accurate(fwd, strike, vol, mat, disc_fact,
                                       call_put, greek);
  } else {
    premium = srt_f_optblknrm_accurate(fwd, strike, vol, mat, disc_fact,
                                       call_put, greek);
  }

  return premium;
}

Err srt_f_optimpvol(double premium, double fwd_price, double strike, double mat,
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
  vol_up = 10;
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
      /* bissection if Newton is out of range  , or not decreasing fast enough
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
      /* the change in newton is acceptable  , take it */
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

Err srt_f_optimpvol_accurate(double premium, double fwd_price, double strike,
                             double mat, double disc, SrtCallPutType call_put,
                             SrtDiffusionType lognormal_normal,
                             double *implied_vol) {
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
  intrinsic = optval_accurate(fwd_price, strike, NULL_VOL, mat, disc, call_put,
                              PREMIUM, lognormal_normal);

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
  vol_up = 10;
  vol_down = 0.000000001;
  vol_shift = 0.0000000001;
  vol_acc = 0.000001;

  prem_up = optval_accurate(fwd_price, strike, vol_up, mat, disc, call_put,
                            PREMIUM, lognormal_normal);

  prem_down = optval_accurate(fwd_price, strike, vol_down, mat, disc, call_put,
                              PREMIUM, lognormal_normal);

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
        prem_up = optval_accurate(fwd_price, strike, vol_up, mat, disc,
                                  call_put, PREMIUM, lognormal_normal);
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

  prem_new = optval_accurate(fwd_price, strike, vol_middle, mat, disc, call_put,
                             PREMIUM, lognormal_normal);

  deriv = (optval_accurate(fwd_price, strike, vol_middle + vol_shift, mat, disc,
                           call_put, PREMIUM, lognormal_normal) -
           prem_new) /
          vol_shift;

  for (i = 0; i <= MAX_ITER; i++) {
    if ((((vol_middle - vol_up) * deriv - (prem_new - premium)) *
             ((vol_middle - vol_down) * deriv - (prem_new - premium)) >=
         0.0) ||
        (fabs(2.0 * (prem_new - premium)) > fabs(vol_diffold * deriv))) {
      /* bissection if Newton is out of range  , or not decreasing fast enough
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
      /* the change in newton is acceptable  , take it */
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

    prem_new = optval_accurate(fwd_price, strike, vol_middle, mat, disc,
                               call_put, PREMIUM, lognormal_normal);

    deriv = (optval_accurate(fwd_price, strike, vol_middle + vol_shift, mat,
                             disc, call_put, PREMIUM, lognormal_normal) -
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

} /* srt_f_optimpvol_accurate() */

/******************************************************************************/

#undef OPT_EPS