
#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/*******************************************************************************
 *			Private Static Declarations
 *******************************************************************************/

static double NI_a;
static double NI_b;
static double NI_d;
static double NI_vol;
static double NI_volsqr;
static double NI_mu;
static int NI_nterms;

/******************************************************************************/

/*******************************************************************************
 *			Private Function Prototypes
 *******************************************************************************/

/******************************************************************************/

static double life_out_time_prob_function(double mat) {
  int i;
  double alpha, beta;
  double d_alpha, d_beta;
  double sum;
  double prob;

  sum = 0.0;
  prob = 0.0;

  for (i = -NI_nterms; i < 0; i++) {
    alpha = exp((2 * NI_mu * (NI_a - i * NI_d)) / NI_volsqr);
    beta = exp((-2 * NI_mu * i * NI_d) / NI_volsqr);

    d_alpha =
        (NI_b - 2 * (NI_a - i * NI_d) - NI_mu * mat) / (NI_vol * sqrt(mat));

    d_beta = (NI_b + 2 * i * NI_d - NI_mu * mat) / (NI_vol * sqrt(mat));

    sum = alpha * norm(d_alpha) - beta * norm(d_beta);

    prob += sum;
  }

  for (i = 0; i <= NI_nterms; i++) {
    alpha = exp((2 * NI_mu * (NI_a - i * NI_d)) / NI_volsqr);
    beta = exp((-2 * NI_mu * i * NI_d) / NI_volsqr);

    d_alpha =
        (NI_b - 2 * (NI_a - i * NI_d) - NI_mu * mat) / (NI_vol * sqrt(mat));

    d_beta = (NI_b + 2 * i * NI_d - NI_mu * mat) / (NI_vol * sqrt(mat));

    sum = alpha * (norm(d_alpha) - 1) - beta * (norm(d_beta) - 1);

    prob += sum;
  }
  /* To renormalize */
  prob /= mat;

  return (prob);
}

/* ======================================================================== */

double srt_f_optlifout(double fwd, double spot, double b_life, double b_out,
                       double vol, double mat, double disc, int nterms,
                       SrtGreekType greek) {
  double mu;
  double integr;
  double premium;
  double answer;
  double shift;
  double a;
  double b;
  double (*f)(double);

  if (mat <= 0) {
    premium = 0.0;
  } else if (((spot > b_out) && (b_out > b_life)) ||
             ((spot < b_out) && (b_out < b_life))) {
    premium = 0.0;
  } else if (((spot > b_life) && (b_life > b_out)) ||
             ((spot < b_life) && (b_life < b_out))) {
    premium = mat * disc;
  } else /* To avoid numerical errors */
      if (fabs(spot - b_out) / spot < 1.0e-04) {
    premium = 0.0;
  } else if (fabs(spot - b_life) / spot < 1.0e-04) {
    premium = mat * disc;
  } else {
    if (b_out < b_life) {
      mu = 1 / mat * log(fwd / spot) - vol * vol / 2;
      a = log(b_out / spot);
      b = log(b_life / spot);
    } else /* Take the inverse Brownian motion */
    {
      mu = -1 / mat * log(fwd / spot) + vol * vol / 2;
      a = -log(b_out / spot);
      b = -log(b_life / spot);
    }

    NI_a = a;
    NI_b = b;
    NI_d = b - a;
    NI_vol = vol;
    NI_volsqr = vol * vol;
    NI_mu = mu;
    NI_nterms = nterms;

    f = life_out_time_prob_function;

    integr = mat * sm_qsimp(f, 0.0000001, mat, 1.0e-06);

    premium = integr * disc;
  }

  switch (greek) {
  case PREMIUM:
    answer = premium;
    break;

  case DELTA_FWD:
    shift = fwd / 10000;
    answer = (srt_f_optlifout(fwd + shift, spot, b_life, b_out, vol, mat, disc,
                              nterms, PREMIUM) -
              premium) /
             shift;
    break;

  case DELTA:
    shift = spot / 10000;
    answer = (srt_f_optlifout(fwd * (1 + shift / spot), spot + shift, b_life,
                              b_out, vol, mat, disc, nterms, PREMIUM) -
              premium) /
             shift;
    break;

  case GAMMA:
    shift = spot / 10000;
    answer = (srt_f_optlifout(fwd * (1 + shift / spot), spot + shift, b_life,
                              b_out, vol, mat, disc, nterms, PREMIUM) +
              srt_f_optlifout(fwd * (1 - shift / spot), spot - shift, b_life,
                              b_out, vol, mat, disc, nterms, PREMIUM) -
              2 * premium) /
             (shift * shift);
    break;

  case VEGA:
    shift = GVOPT.vol_add;
    answer = (srt_f_optlifout(fwd, spot, b_life, b_out, vol + shift, mat, disc,
                              nterms, PREMIUM) -
              premium) /
             shift;
    break;

  case THETA:
    shift = YEARS_IN_DAY;
    answer = (srt_f_optlifout(fwd, spot, b_life, b_out, vol, mat - shift,
                              disc * exp(-shift * log(disc) / mat), nterms,
                              PREMIUM) -
              premium);
    break;

  default:
    answer = UNKNOWN_GREEK;
  }

  return answer;

} /* END of srt_f_optlifout */
