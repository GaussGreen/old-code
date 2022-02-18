/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"

#include "utallhdr.h"

/******************************************************************************/

static double proba_joint_spot_min_max(
    double mu,
    double sig,
    double mat,
    double k,
    double b,
    double cp, /*1 if call:S>K  -1 if put:S<K */
    double du  /*1 if down:min>B  -1 if up:Max<B*/
);

static double proba_min_max(
    double mu, double sig, double mat, double b, double du /*1 if down:min>B  -1 if up:Max<B*/
);

/******************************************************************************
        computes the following probabilities:
                Prob ( x(T) > k ; m(T) > b ): 	cp = 1; du = 1
                Prob ( x(T) < k ; m(T) > b ): 	cp = -1; du = 1
                Prob ( x(T) > k ; M(T) < b ): 	cp = 1; du = -1
                Prob ( x(T) < k ; M(T) < b ): 	cp = -1; du = -1
******************************************************************************/
static double proba_joint_spot_min_max(
    double mu,
    double sig,
    double mat,
    double k,
    double b,
    double cp, /*1 if call:S>K  -1 if put:S<K */
    double du  /*1 if down:min>B  -1 if up:Max<B*/
)
{
    double prob;
    double sig_sqrt = sig * sqrt(mat);
    double d1, d2, alpha;

    /* Call (S>K) Down&Out with B>K: integration starts at b */
    /* Put  (S<K) Up&Out   with K>B: integration starts at b */
    if ((cp * du > 0.0) && (cp * (k - b) < 0))
        k = b;

    d1    = du * (-k + mu * mat) / sig_sqrt;
    alpha = exp(2 * mu * b / (sig * sig));
    d2    = du * (-k + 2 * b + mu * mat) / sig_sqrt;

    prob = norm(d1) - alpha * norm(d2);

    if (cp * du < 0.0) /* Call Up&Out or Put Down&Out */
    {
        prob *= -1;
        prob += proba_min_max(mu, sig, mat, b, du);
    }

    return (prob);
}

/******************************************************************************
        computes the following probabilities:
                        Prob ( m(T) > b ):	du = 1
                        Prob ( M(T) < b ):	du = -1
*******************************************************************************/
static double proba_min_max(
    double mu, double sig, double mat, double b, double du /*1 if down:min>B  -1 if up:Max<B*/
)
{
    double prob;
    double sig_sqrt = sig * sqrt(mat);

    prob = norm(du * (-b + mu * mat) / sig_sqrt);
    prob -= exp(2 * mu * b / (sig * sig)) * norm(du * (b + mu * mat) / sig_sqrt);

    return (prob);
}

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optexting()
 *
 * PURPOSE     	: extinguishable options for stocks, bonds & swaps
 *
 * DESCRIPTION  	: XX
 *
 * CALLS		: norm() in gen_math.c
 *
 * PARAMETERS  	: fwd_price 	- forward price of underlying
 *             	: spot     	- spot price of underlying
 *             	: strike   	- strike price
 *             	: barrier   	- barrier
 *             	: vol         	- volatility of underlying
 *             	: mat        	- maturity of underlying
 *             	: disc       	- discount factor
 *             	: call_put      - type of option: 0 call, 1 put
 *             	: down_up       - type of extinguish : 0 down, 1 up
 *             	: greek      	- result wanted (premium or greeks)
 *
 * RETURNS      	: premium    	- price of the option or greeks
 *
 *******************************************************************************/

double srt_f_optexting(
    double         fwd,
    double         spot,
    double         strike,
    double         barrier,
    double         vol,
    double         mat,
    double         disc,
    SrtCallPutType call_put,
    SrtBarrierType down_up,
    SrtGreekType   greek)
{
    double premium, shift, cp, du, answer;
    double mu, mu_1, mu_2, k, b;

    cp = (call_put == SRT_CALL) ? 1 : -1;
    du = (down_up == SRT_DOWN) ? 1 : -1;

    /* Spot should better be above the barrier for a Down&Out*/
    /* and better be below the barrier for an Up&Out */
    if ((spot - barrier) * du < 0)
        premium = 0.0;
    else
        /* The price for an option at maturity is MAX(S-K,0) for a call */
        /* and MAX(K-S,0) for a put */
        if (mat <= 0.0)
    {
        if ((spot - strike) * cp > 0)
            premium = cp * (spot - strike);
        else
            premium = 0.0;
    }
    else
        /* if the vol is equal to 0  then the spot depends only from the constant rate */
        /* at maturity the spot is equal to the forward */
        if (vol == 0.0)
    {
        if (((fwd - barrier) * du > 0) && ((fwd - strike) * cp > 0))
            premium = cp * (fwd - strike) * disc;
        else
            premium = 0.0;
    }
    else
        /* Whatever happens the spot will always cross the barrier */
        /* or give a payoff nul*/
        if ((cp * (strike - barrier) > 0) && (du * cp < 0))
        premium = 0.0;
    else
    /* Usual case */
    {
        mu = log(fwd / spot) / mat;

        k = log(strike / spot);
        b = log(barrier / spot);

        mu_2 = mu + 0.5 * vol * vol;
        mu_1 = mu - 0.5 * vol * vol;

        premium = fwd * proba_joint_spot_min_max(mu_2, vol, mat, k, b, cp, du);

        premium -= strike * proba_joint_spot_min_max(mu_1, vol, mat, k, b, cp, du);
        premium *= cp * disc;
    }

    switch (greek)
    {
    case PREMIUM: /*** PREMIUM ***/
        answer = premium;
        break;

    case DELTA_FWD: /*** DELTA FWD ***/
        shift = fwd / 10000;
        answer =
            (srt_f_optexting(
                 fwd + shift, spot, strike, barrier, vol, mat, disc, call_put, down_up, PREMIUM) -
             premium) /
            shift;
        break;

    case DELTA: /*** DELTA SPOT + FWD ***/
        shift  = spot / 10000;
        answer = (srt_f_optexting(
                      fwd * (1 + shift / spot),
                      spot + shift,
                      strike,
                      barrier,
                      vol,
                      mat,
                      disc,
                      call_put,
                      down_up,
                      PREMIUM) -
                  premium) /
                 shift;
        break;

    case GAMMA: /*** GAMMA ***/
        shift  = spot / 10000;
        answer = srt_f_optexting(
            fwd * (1 + shift / spot),
            spot + shift,
            strike,
            barrier,
            vol,
            mat,
            disc,
            call_put,
            down_up,
            PREMIUM);
        answer += srt_f_optexting(
            fwd * (1 - shift / spot),
            spot - shift,
            strike,
            barrier,
            vol,
            mat,
            disc,
            call_put,
            down_up,
            PREMIUM);
        answer -= 2 * premium;
        answer /= shift * shift;
        break;

    case VEGA: /*** VEGA ***/
        shift = GVOPT.vol_add;
        answer =
            (srt_f_optexting(
                 fwd, spot, strike, barrier, vol + shift, mat, disc, call_put, down_up, PREMIUM) -
             premium) /
            shift;
        break;

    case THETA: /*** THETA  ***/
        shift  = YEARS_IN_DAY;
        answer = srt_f_optexting(
                     fwd,
                     spot,
                     strike,
                     barrier,
                     vol,
                     mat - shift,
                     disc * exp(-shift * log(disc) / mat),
                     call_put,
                     down_up,
                     PREMIUM) -
                 premium;
        break;

    default:
        answer = UNKNOWN_GREEK;
        break;
    }

    return (answer);

} /* END srt_f_optexting() */

/******************************************************************************/
