/*******************************************************************************
**
**	dbarrier_fct.cxx
**
*******************************************************************************/
/*******************************************************************************
        DBARRIER - Double Barrier
        Calculates the density        , probability and option value for
        staying within two barrier levels between time T1 (i.e. forward start)
        and maturity        , for either risk neutral (measure_flag = 0) or
measure associated with the underlying (measure_flag = 1) .

        add test for mat=0 case; K L Chau 16/10/95
*******************************************************************************/

/* ==========================================================================
   include files
   ========================================================================== */

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/* ==========================================================================
   dbarrier_dens(...)
        - Gaussian Method:
                term used in the infinite sum
   ========================================================================== */

static double dbarrier_dens(double bound, double x1_n, double x2_n, double mu,
                            double vol, double mat) {
  double dens, coeff1, d1, coeff2, d2;
  double sig_sqrt_t = vol * sqrt(mat);
  double sig_sqr = vol * vol;

  coeff1 = exp(mu * x1_n / sig_sqr);
  d1 = (bound - x1_n - mu * mat) / sig_sqrt_t;

  coeff2 = exp(mu * x2_n / sig_sqr);
  d2 = (bound - x2_n - mu * mat) / sig_sqrt_t;

  dens = coeff1 * norm(d1) - coeff2 * norm(d2);

  return (dens);
}

/* ==========================================================================
   dbarrier_prob(...)
        - Gaussian Method
   This function computes the following probabilities for a drifted
   normal Brownian motion:
                if (call_put == SRT_CALL):
                                Prob[ X(T) > k; min X(t) > b_do ; max X(t) <
   b_up ] if (call_put == SRT_PUT): Prob[ X(T) < k; min X(t) > b_do ; max X(t) <
   b_up ]
                                                                        ==
                                Prob[ -X(T) > -k; max -X(t) < -b_do ; min -X(t)
   > -b_up ]
   ========================================================================== */

double dbarrier_prob(double k, double b_up, double b_do, double mu, double vol,
                     double mat, SrtCallPutType call_put, int nb_terms) {
  double x1_n, x2_n, prob, sum;
  double bound_do, bound_up;
  double rubbish;
  int i;

  if (mat <= 0) {
    if (b_up * b_do > 0.0)
      return 0.0;
    else
      return 1.0;
  }

  if (call_put == SRT_PUT) {
    k = -k;
    mu = -mu;
    rubbish = b_up;
    b_up = -b_do;
    b_do = -rubbish;
  }

  if (k >= b_up)
    return 0.0;
  else if (k > b_do) {
    bound_up = b_up;
    bound_do = k;
  } else {
    bound_up = b_up;
    bound_do = b_do;
  }

  prob = 0.0;
  for (i = -nb_terms; i <= nb_terms; i++) {
    x1_n = 2 * i * (b_up - b_do);
    x2_n = 2 * b_up - x1_n;

    sum = dbarrier_dens(bound_up, x1_n, x2_n, mu, vol, mat);
    sum -= dbarrier_dens(bound_do, x1_n, x2_n, mu, vol, mat);

    prob += sum;
  }

  return prob;
}

/* ==========================================================================
   srt_f_optdblbar(...)
        - Gaussian Method
   Price of a double barrier (extinguish) option
   ========================================================================== */

double srt_f_optdblbar(double spot, double forward, double strike,
                       double bar_do, double bar_up, double vol, double mat,
                       double disc, SrtCallPutType call_put, int nb_term,
                       SrtGreekType greek) {
  double b_up, b_do, k, mu_k, mu_s, premium, answer;
  double drift, proba_k, proba_s, cp, shift;

  if (mat < 0)
    return 0;

  cp = (call_put == SRT_CALL) ? 1 : -1;

  /* the spot should be between the two barriers */
  if ((spot <= bar_do) || (spot >= bar_up)) {
    premium = 0.0;
  } else
      /* To get a Payoff        , the spot has to cross the barrier */
      if (((strike >= bar_up) && (call_put == SRT_CALL)) ||
          ((strike <= bar_do) && (call_put == SRT_PUT))) {
    premium = 0.0;
  } else /* opt maturity is 0        , give intrinsic value */
      if ((mat == 0.0) || (vol == 0.0)) {
    /* the case at which spot is outside range is already treated
       above        , so it is only necessary to calculate intrinsic for the
       simple cases here */

    if (call_put == SRT_CALL) {
      premium = (spot > strike ? spot - strike : 0.0);
    } else {
      premium = (strike > spot ? strike - spot : 0.0);
    }
  } else {
    b_up = log(bar_up / spot);
    b_do = log(bar_do / spot);
    k = log(strike / spot);

    drift = log(forward / spot) / mat;
    mu_k = drift - vol * vol / 2;
    mu_s = drift + vol * vol / 2;

    proba_k = dbarrier_prob(k, b_up, b_do, mu_k, vol, mat, call_put, nb_term);
    proba_s = dbarrier_prob(k, b_up, b_do, mu_s, vol, mat, call_put, nb_term);

    premium = cp * disc * (forward * proba_s - strike * proba_k);
  }

  if (mat != 0) {
    switch (greek) {
    case PREMIUM: /*** PREMIUM ***/
      answer = premium;
      break;

    case DELTA_FWD: /*** DELTA FWD ***/
      shift = forward / 10000;
      answer = (srt_f_optdblbar(spot, forward + shift, strike, bar_do, bar_up,
                                vol, mat, disc, call_put, nb_term, PREMIUM) -
                premium) /
               shift;
      break;

    case DELTA: /*** DELTA SPOT + FWD ***/
      shift = spot / 10000;
      answer = (srt_f_optdblbar(spot + shift, forward * (1 + shift / spot),
                                strike, bar_do, bar_up, vol, mat, disc,
                                call_put, nb_term, PREMIUM) -
                premium) /
               shift;
      break;

    case GAMMA: /*** GAMMA ***/
      shift = spot / 10000;
      answer = srt_f_optdblbar(spot + shift, forward * (1 + shift / spot),
                               strike, bar_do, bar_up, vol, mat, disc, call_put,
                               nb_term, PREMIUM);
      answer += srt_f_optdblbar(spot - shift, forward * (1 - shift / spot),
                                strike, bar_do, bar_up, vol, mat, disc,
                                call_put, nb_term, PREMIUM);
      answer -= 2 * premium;
      answer /= shift * shift;
      break;

    case VEGA: /*** VEGA ***/
      shift = GVOPT.vol_add;
      answer =
          (srt_f_optdblbar(spot, forward, strike, bar_do, bar_up, vol + shift,
                           mat, disc, call_put, nb_term, PREMIUM) -
           premium) /
          shift;
      break;

    case THETA: /*** THETA  ***/
      shift = YEARS_IN_DAY;
      answer =
          srt_f_optdblbar(spot, forward, strike, bar_do, bar_up, vol,
                          mat - shift, disc * exp(-shift * log(disc) / mat),
                          call_put, nb_term, PREMIUM) -
          premium;
      break;

    default:
      answer = UNKNOWN_GREEK;
      break;
    }
  } else /* maturity = 0 case */
  {
    switch (greek) {
    case PREMIUM:
      answer = premium;
      break;

    case DELTA:
    case DELTA_FWD:
      /* the delta when spot = strike is assumed to be 0        , instead of
       * infinity
       */
      if (premium == 0)
        answer = 0.0;
      else
        answer = 1.0;
      break;

    case GAMMA:
      /* again        , ignore discontinuous delta at spot = strike and assume
         gamma is 0 throughout */
    case VEGA:
    case THETA:
      answer = 0.0;
      break;

    default:
      answer = UNKNOWN_GREEK;
      break;
    }
  }

  return (answer);

} /* END of srt_f_optdblbar */
