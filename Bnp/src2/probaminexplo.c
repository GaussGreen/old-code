/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "math.h"
#include "num_h_proba.h"

#include "utallhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: proba_min_explo_fct(...)
 *
 * DESCRIPTION  	: Returns the probability to get 1 unit of currency as soon as
 *                 the spot becomes lower than the barrier during mat.
 *
 * CALLS			: gauss(...) in gen_math.c
 *				  norm(...)  in gen_math.c
 *
 * PARAMETERS
 *	INPUT		: spot			- spot underlying price
 *				: barrier      	- barrier level to stay under
 *				: vol         	- annual volatility
 *				: mat         	- maturity, in years
 *				: disc        	- discount factor for the domestic currency
 *
 * RETURNS		: proba			- probability
 *
 *******************************************************************************/

double proba_min_explo_fct(double spot, double barrier, double vol, double mat, double disc)
{
    double mu, mu1, mu2, l, d1, d2, price = 0.0, coeff;

    if (spot <= barrier)
        return (1);

    mu1 = -log(disc) / mat - vol * vol / 2;
    mu2 = sqrt(mu1 * mu1 - 2 * vol * vol * log(disc) / mat);
    mu  = mu1 - mu2;

    l = log(barrier / spot);

    coeff = exp(2 * mu2 * l / vol / vol);

    d1 = (-l + mu2 * mat) / vol / sqrt(mat);
    d2 = (l + mu2 * mat) / vol / sqrt(mat);

    price = 1 - (norm(d1) - coeff * norm(d2));
    price *= exp(l * mu / vol / vol);

    return (price);
}
