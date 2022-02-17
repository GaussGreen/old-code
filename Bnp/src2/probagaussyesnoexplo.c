/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include <math.h"
#include <num_h_allhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: probag_yes_no_explo_fct(...)
 *
 * DESCRIPTION  	: Returns the probability to get 1 unit of currency as
 *soon as the spot hits barrier_yes without having hit barrier_no before that.
 *
 * CALLS			: gauss(...) in gen_math.c
 *				  norm(...)  in gen_math.c
 *
 * PARAMETERS
 *	INPUT		: spot			- spot underlying price
 *				: barrier_yes   - barrier level to hit
 *				: barrier_no    - barrier level not to hit
 *				: vol         	- annual volatility
 *				: mat         	- maturity in years
 *				: disc        	- discount factor for the domestic
 *currency : num_terms		- number of terms to compute in the series
 *
 * RETURNS		: proba			- probability
 *
 *******************************************************************************/

double probag_yes_no_explo_fct(double spot, double b_yes, double b_no,
                               double vol, double mat, double disc,
                               int nb_term) {
  double a, b, price = 0.0, yn, coeff_tot, coeff;
  double mu, term, d1, d2, mu2;
  int i, mult, mult2;

  if (b_yes < b_no) {
    if (b_yes >= spot)
      return (1);

    if (b_no <= spot)
      return (0);

    mult2 = 1;
  } else {
    if (b_yes <= spot)
      return (1);

    if (b_no >= spot)
      return (0);

    mult2 = -1;
  }

  a = log(b_yes / spot);
  b = log(spot / b_no);

  mu = -log(disc) / mat - vol * vol / 2;
  mu2 = -sqrt(mu * mu - 2 * vol * vol * log(disc) / mat);

  for (i = -nb_term; i <= nb_term; i++) {
    yn = -(2 * i + 1) * (a + b) + b;
    if (yn > 0)
      mult = 1;
    else
      mult = -1;

    coeff = exp(-2 * mu2 * yn / vol / vol);

    d1 = mult * (-yn - mu2 * mat) / vol / sqrt(mat);
    d2 = mult * (-yn + mu2 * mat) / vol / sqrt(mat);

    coeff_tot = exp(-mu * 2 * i * (a + b) / vol / vol) *
                exp(yn * (mu2 - mu) / vol / vol);

    if (yn == 0)
      term = -0.5;
    else
      term = mult * (norm(d1) + coeff * norm(d2));

    price += term * coeff_tot;
  }

  return (mult2 * price);
}