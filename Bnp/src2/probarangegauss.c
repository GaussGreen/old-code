/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <math.h"
#include <num_h_proba.h"

/******************************************************************************/

static double dens_gauss_n(double var, double x1_n, double x2_n, double mu,
                           double vol, double mat) {
  double result, coeff1, d1, coeff2, d2;

  coeff1 = exp(mu * x1_n / (vol * vol));
  d1 = (var - x1_n - mu * mat) / (vol * sqrt(mat));

  coeff2 = exp(mu * x2_n / (vol * vol));
  d2 = (var - x2_n - mu * mat) / (vol * sqrt(mat));

  result = coeff1 * norm(d1) - coeff2 * norm(d2);

  return (result);
}

/*******************************************************************************
 *
 * FUNCTION     	: proba_gauss_fct(...)
 *
 * DESCRIPTION  	: Returns the probability to get 1 unit of currency at
 *maturity if the spot stays under barrier_up and above barrier_do during mat.
 *
 * CALLS			: gauss(...) in gen_math.c
 *				  norm(...)  in gen_math.c
 *
 * PARAMETERS
 *	INPUT		: spot			- spot underlying price
 *				: barrier_up    - barrier level to stay under
 *				: barrier_do    - barrier level to stay above
 *				: vol         	- annual volatility
 *				: mat         	- maturity  , in years
 *				: disc        	- discount factor for the domestic
 *currency : num_terms		- number of terms to compute in the series
 *
 * RETURNS		: proba			- probability
 *
 *******************************************************************************/

double proba_gauss_fct(double spot, double b_up, double b_do, double vol,
                       double mat, double disc, int nb_term) {
  double x1_n, x2_n, a, b, mu, bound_do, bound_up, price = 0.0;
  double price_k;
  int i;

  if ((spot >= b_up) || (spot <= b_do))
    return (0);

  a = log(b_up / spot);
  b = log(spot / b_do);

  mu = -log(disc) / mat - vol * vol / 2;

  bound_up = a;
  bound_do = -b;

  for (i = -nb_term; i <= nb_term; i++) {
    x1_n = 2 * i * (a + b);
    x2_n = 2 * a - x1_n;

    price_k = dens_gauss_n(bound_up, x1_n, x2_n, mu, vol, mat) -
              dens_gauss_n(bound_do, x1_n, x2_n, mu, vol, mat);

    price += disc * price_k;
  }

  return (price);
}

#undef GAUSS_SIMP_ESP