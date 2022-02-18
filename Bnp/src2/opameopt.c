/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"

#include "utallhdr.h"

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_optameopt(...)
 *
 * PURPOSE      	: pricing an american option
 *
 * DESCRIPTION  	: ??
 *
 * CALLS			: SELF
 *
 * PARAMETERS   	: spot     	- spot price
 *              	: forward   	- forward price of underlying
 *              	: strike   	- strike price
 *				: disc_fact	- discount factor
 *              	: vol         	- volatility of underlying
 *              	: opt_mat 	- maturity of option
 *              	: call_put 	- type of option: 0 call, 1 put
 *              	: step_num 	- number of steps for ??
 *              	: greek     	- ??
 *
 * RETURNS      	: premium    	- ??
 *
 *******************************************************************************/

double srt_f_optameopt(
    double         spot,
    double         forward,
    double         strike,
    double         disc_fact,
    double         vol,
    double         opt_mat,
    SrtCallPutType call_put,
    int            step_num,
    SrtGreekType   greek)
{
    int    i, n;
    double p_u, p_d, delta_t;
    double a, bb, uu, dd;
    double drift;
    double asset;
    double amer_disc_premium;
    double euro_disc_premium;
    double factor_disc_premium;
    double intrinsic;
    double BSprice;
    double BSdelta;
    double BSgamma;
    double BStheta;
    int    bPayatTheExFlag = 0;
    double answer;

    double  diff_up;
    double  diff_do;
    double  amer_price;
    double  amer_delta;
    double  amer_gamma;
    double  amer_theta;
    double  amer_vega;
    double  euro_price;
    double  euro_delta;
    double  euro_gamma;
    double  euro_theta;
    double* euro_premium;
    double* amer_premium;
    double* disc_factor = NULL;
    int     old_step_num;
    double  df_n;
    double  r_new;
    double  theta_1d;
    double  cp;

    old_step_num = step_num;

    if (step_num < 0)
    {
        bPayatTheExFlag = 1;
        step_num        = -step_num;
    }

    if (call_put == SRT_STRADDLE)
    {
        answer = srt_f_optameopt(
                     spot, forward, strike, disc_fact, vol, opt_mat, SRT_CALL, step_num, greek) +
                 srt_f_optameopt(
                     spot, forward, strike, disc_fact, vol, opt_mat, SRT_PUT, step_num, greek);
        return (answer);
    }

    cp = (call_put == SRT_CALL) ? 1 : -1;

    if (opt_mat < 0)
    {
        return 0;
    }
    else if ((opt_mat <= 1.0e-10) || (vol <= 1.0e-9))
    {
        if (cp * (spot - strike) > 0)
        {
            if (greek == PREMIUM)
            {
                answer = cp * (spot - strike);
            }
            else if (greek == DELTA)
            {
                answer = cp;
            }
            else
            {
                answer = 0;
            }
        }
        else if (spot == strike)
        {
            if (greek == DELTA)
            {
                answer = cp * 0.5;
            }
            else
            {
                answer = 0;
            }
        }
        else
        {
            answer = 0;
        }

        return (answer);
    }
    else
    {
        /* Control variate on a European option: theroetical values */
        if (bPayatTheExFlag)
        {
            BSprice = srt_f_optblksch(forward, strike, vol, opt_mat, 1, call_put, PREMIUM);
            BSdelta = srt_f_optblksch(forward, strike, vol, opt_mat, 1, call_put, DELTA);
            BSgamma = srt_f_optblksch(forward, strike, vol, opt_mat, 1, call_put, GAMMA);
            BStheta = srt_f_optblksch(forward, strike, vol, opt_mat, 1, call_put, THETA);
        }
        else
        {
            BSprice = srt_f_optblksch(forward, strike, vol, opt_mat, disc_fact, call_put, PREMIUM);
            BSdelta = srt_f_optblksch(forward, strike, vol, opt_mat, disc_fact, call_put, DELTA);
            BSgamma = srt_f_optblksch(forward, strike, vol, opt_mat, disc_fact, call_put, GAMMA);
            BStheta = srt_f_optblksch(forward, strike, vol, opt_mat, disc_fact, call_put, THETA);
        }

        /* Memory Allocation for the full tree grid */
        if ((amer_premium = dvector(0, step_num + 3)) == NULL)
            return MEMORY_ERR;
        if ((euro_premium = dvector(0, step_num + 3)) == NULL)
            return MEMORY_ERR;
        if (bPayatTheExFlag)
        {
            if ((disc_factor = dvector(0, step_num + 3)) == NULL)
                return MEMORY_ERR;
        }

        /* Computation of drift, variance,... */
        delta_t = opt_mat / (double)step_num;
        drift   = log(forward / spot) / opt_mat;
        a       = exp(drift * delta_t);

        bb = a * a * (exp(vol * vol * delta_t) - 1);
        uu = ((a * a + bb + 1) + sqrt((a * a + bb + 1) * (a * a + bb + 1) - 4 * a * a)) / (2 * a);

        /* Probabilities */
        dd  = 1 / uu;
        p_u = (a - dd) / (uu - dd);
        p_d = 1 - p_u;

        /* Use for Greek calculation */
        diff_up = spot * (pow(uu, 2) - 1);
        diff_do = spot * (1 - pow(dd, 2));

        /* Initialisation of price at top of the tree */
        asset = spot * pow(uu, (double)(step_num + 2));

        /* Payoff at option expiry */
        for (i = 0; i < step_num + 3; i++)
        {
            intrinsic = cp * (asset - strike);

            if (intrinsic > 0)
            {
                euro_premium[i] = intrinsic;
                amer_premium[i] = intrinsic;
            }
            else
            {
                euro_premium[i] = 0;
                amer_premium[i] = 0;
            }
            asset /= uu * uu;
            if (bPayatTheExFlag)
                disc_factor[i] = 1.0;
        }

        /* Discount factor from step t to step t+delta_t */
        df_n = pow(disc_fact, 1 / (double)step_num);

        /* Backward induction from step T-delta_t to 0 */
        for (n = step_num + 1; n >= 1; n--)
        {
            asset = spot * pow(uu, (double)n);

            /* Loop on all nodes at the step */
            for (i = 0; i <= n; i++)
            {
                amer_disc_premium = (p_u * amer_premium[i] + p_d * amer_premium[i + 1]) * df_n;
                euro_disc_premium = (p_u * euro_premium[i] + p_d * euro_premium[i + 1]) * df_n;
                if (bPayatTheExFlag)
                    factor_disc_premium = (p_u * disc_factor[i] + p_d * disc_factor[i + 1]) * df_n;

                intrinsic = cp * (asset - strike);

                if (amer_disc_premium > intrinsic)
                {
                    amer_premium[i] = amer_disc_premium;
                    if (bPayatTheExFlag)
                        disc_factor[i] = factor_disc_premium;
                }
                else
                {
                    amer_premium[i] = intrinsic;
                    if (bPayatTheExFlag)
                        disc_factor[i] = 1.0;
                }
                euro_premium[i] = euro_disc_premium;

                asset /= uu * uu;
            }

            if (n == 4)
            {
                amer_theta = amer_premium[2];
                euro_theta = euro_premium[2];
            }

            if (n == 2)
            {
                if (bPayatTheExFlag)
                {
                    amer_premium[0] /= disc_factor[0];
                    amer_premium[1] /= disc_factor[1];
                    amer_premium[2] /= disc_factor[2];

                    euro_premium[0] /= disc_factor[0];
                    euro_premium[1] /= disc_factor[1];
                    euro_premium[2] /= disc_factor[2];
                }

                euro_price = euro_premium[1];
                amer_price = amer_premium[1];

                amer_delta = (amer_premium[0] - amer_premium[2]) / (diff_up + diff_do);
                euro_delta = (euro_premium[0] - euro_premium[2]) / (diff_up + diff_do);

                amer_gamma = (amer_premium[0] - amer_premium[1]) / diff_up;
                amer_gamma -= (amer_premium[1] - amer_premium[2]) / diff_do;
                amer_gamma /= 0.5 * (diff_up + diff_do);

                euro_gamma = (euro_premium[0] - euro_premium[1]) / diff_up;
                euro_gamma -= (euro_premium[1] - euro_premium[2]) / diff_do;
                euro_gamma /= 0.5 * (diff_up + diff_do);
            }
        }

        switch (greek)
        {
        case PREMIUM:
            answer = amer_price + BSprice - euro_price;
            break;

        case DELTA:
            answer = amer_delta + BSdelta - euro_delta;
            break;

        case GAMMA:
            answer = amer_gamma + BSgamma - euro_gamma;
            break;

        case THETA:
            euro_theta = (euro_price - euro_theta) / (2 * delta_t);
            amer_theta = (amer_price - amer_theta) / (2 * delta_t);
            answer     = amer_theta + BStheta * 365.0 - euro_theta;
            break;

        case THETA_1D:
            r_new     = (1.0 / disc_fact - 1.0) / opt_mat;
            disc_fact = 1.0 / (1.0 + r_new * (opt_mat - YEARS_IN_DAY));
            theta_1d  = srt_f_optameopt(
                spot,
                forward,
                strike,
                disc_fact,
                vol,
                opt_mat - YEARS_IN_DAY,
                call_put,
                old_step_num,
                PREMIUM);
            theta_1d -= amer_price;
            answer = theta_1d;
            break;

        case VEGA:

            amer_vega = srt_f_optameopt(
                spot,
                forward,
                strike,
                disc_fact,
                vol + GVOPT.vol_add,
                opt_mat,
                call_put,
                old_step_num,
                PREMIUM);
            amer_vega = amer_vega - amer_price + euro_price - BSprice;
            answer    = amer_vega;
            break;

        default:
            answer = UNKNOWN_GREEK;
            break;
        }

        if (disc_factor)
        {
            free_dvector(disc_factor, 0, step_num + 3);
            disc_factor = NULL;
        }
        free_dvector(amer_premium, 0, step_num + 3);
        amer_premium = NULL;
        free_dvector(euro_premium, 0, step_num + 3);
        euro_premium = NULL;
    }

    return (answer);

} /* END srt_f_optameopt(...) */

/******************************************************************************/
