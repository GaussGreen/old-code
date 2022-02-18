/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "math.h"
#include "num_h_simpso.h"

#include "utallhdr.h"

#define GAUSS_SIMP_EPS 1e-06

static double NI_var;
static double NI_vol;
static double NI_mu;
static double NI_mu2;

/******************************************************************************/

static double int_g(double t)
{
    double expo;
    double result;

    expo   = exp(-pow(NI_var - NI_mu2 * t, 2) / (2 * NI_vol * NI_vol * t));
    result = (NI_var + NI_mu * t) / (2 * NI_vol * t * sqrt(2 * SRT_PI * t)) * expo;

    return (result);
}

static double dens_gauss_n_explo(
    double var, double x1_n, double x2_n, double mu, double mu2, double vol, double mat)
{
    double result, coeff1, integ1, coeff2, integ2;
    double (*f)(double);

    coeff1 = exp(mu * x1_n / (vol * vol)) * exp((mu - mu2) * (var - x1_n) / (vol * vol));

    NI_var = var - x1_n;
    f      = int_g;
    integ1 = sm_qsimp(f, GAUSS_SIMP_EPS, mat, GAUSS_SIMP_EPS);

    coeff2 = exp(mu * x2_n / (vol * vol)) * exp((mu - mu2) * (var - x2_n) / (vol * vol));

    NI_var = var - x2_n;
    f      = int_g;
    integ2 = sm_qsimp(f, GAUSS_SIMP_EPS, mat, GAUSS_SIMP_EPS);

    result = coeff1 * integ1 - coeff2 * integ2;

    return (result);
}

/*******************************************************************************
 *
 * FUNCTION     	: proba_gauss_explo_fct(...)
 *
 * DESCRIPTION  	: Returns the probability to get 1 unit of domestic currency (if
 *				  dom/for=0) or 1 unit of foreign currency(if dom/for=1) as soon as
 *                 the spot hits barrier_up or barrier_do during mat.
 *
 * CALLS			: gauss(...) in gen_math.c
 *				  norm(...)  in gen_math.c
 *
 * PARAMETERS
 *	INPUT		: spot			- spot underlying price
 *				: barrier_up    - barrier level to stay under
 *				: barrier_do    - barrier level to stay above
 *				: vol         	- annual volatility
 *				: mat         	- maturity, in years
 *				: disc        	- discount factor for the domestic currency
 *				: num_terms		- number of terms to compute in the series
 *
 * RETURNS		: proba			- probability
 *
 *******************************************************************************/

double proba_gauss_explo_fct(
    double spot, double b_up, double b_do, double vol, double mat, double disc, int nb_term)
{
    double x1_n, x2_n, a, b, mu, bound_do, bound_up, price = 0.0;
    double price_k, mu2;
    int    i;

    if ((spot >= b_up) || (spot <= b_do))
        return (1);

    a = log(b_up / spot);
    b = log(spot / b_do);

    mu  = -log(disc) / mat - vol * vol / 2;
    mu2 = sqrt(mu * mu - 2 * vol * vol * log(disc) / mat);

    NI_vol = vol;
    NI_mu  = mu;
    NI_mu2 = mu2;

    bound_up = a;
    bound_do = -b;

    for (i = -nb_term; i <= nb_term; i++)
    {
        x1_n = 2 * i * (a + b);
        x2_n = 2 * a - x1_n;

        price_k = dens_gauss_n_explo(bound_up, x1_n, x2_n, mu, mu2, vol, mat) -
                  dens_gauss_n_explo(bound_do, x1_n, x2_n, mu, mu2, vol, mat);

        price += price_k;
    }

    return (price);
}

#undef GAUSS_SIMP_ESP