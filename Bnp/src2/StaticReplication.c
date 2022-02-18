//================================================================
//	Replication of any payoff by giving:
//		a function to evaluate the payoff at any point
//		a function to evaluate the derivative of the payoff at any point
//		a list of tangent points
//
//	Creation of a robust algorithm to find the way to spread the tangent points if the function
//	has a constant convexity.
//
//	11.11.2004 J. Dinh
//================================================================

#include "StaticReplication.h"

#include "CPDVol.h"
#include "math.h"
#include "opfnctns.h"

//================================================================
//
// Replication of any positive payoff
//		by using a constant + list of calls + Put or Call for the extrapolation
//
//================================================================
Err Static_payoff_replication(
    double* tangent_points,
    int     n_tangent_points,
    double  cvx_strike,
    int     is_left_flat,
    int     is_right_flat,
    Err (*payoff)(double S, double K, double* out),
    Err (*deriv)(double S, double K, double* out),
    double* constant,
    double* strikes,
    double* notionals,
    int*    is_call)
{
    int     i;
    double *payoff_val = NULL, *deriv_val = NULL;
    Err     err = NULL;

    payoff_val = calloc(n_tangent_points, sizeof(double));
    deriv_val  = calloc(n_tangent_points, sizeof(double));
    if (!deriv_val || !payoff_val)
    {
        err = "Static_payoff_replication: memory allocation failure";
        goto FREE_RETURN;
    }

    // fill the infos
    for (i = 0; i < n_tangent_points; i++)
    {
        err = payoff(tangent_points[i], cvx_strike, &(payoff_val[i]));
        err = deriv(tangent_points[i], cvx_strike, &(deriv_val[i]));
    }

    // get the constant
    *constant = payoff_val[0];

    // get the strikes (n+1) given by the intersection of the tangents
    strikes[0]                = tangent_points[0];
    strikes[n_tangent_points] = tangent_points[n_tangent_points - 1];
    for (i = 0; i < n_tangent_points - 1; i++)
    {
        strikes[i + 1] =
            ((deriv_val[i] * tangent_points[i] - payoff_val[i]) -
             (deriv_val[i + 1] * tangent_points[i + 1] - payoff_val[i + 1]));
        strikes[i + 1] /= deriv_val[i] - deriv_val[i + 1];
    }

    // get the notionals to be flat left & right
    notionals[0] = deriv_val[0];
    for (i = 0; i < n_tangent_points - 1; i++)
        notionals[i + 1] = deriv_val[i + 1] - deriv_val[i];

    notionals[n_tangent_points] = -deriv_val[n_tangent_points - 1];

    // call or put replication
    for (i = 0; i < n_tangent_points + 1; i++)
        is_call[i] = 1;

    // Handels the left extrapolation
    if (!is_left_flat)
    {
        is_call[n_tangent_points + 1]   = 0;
        strikes[n_tangent_points + 1]   = strikes[0];
        notionals[n_tangent_points + 1] = -notionals[0];
    }
    else
    {
        is_call[n_tangent_points + 1]   = 0;
        strikes[n_tangent_points + 1]   = strikes[0];
        notionals[n_tangent_points + 1] = 0.0;
    }

    // Handels the left extrapolation
    if (!is_right_flat && payoff_val[n_tangent_points - 1] > 1.0e-10 &&
        fabs(deriv_val[n_tangent_points - 1]) > 1.0e-10)
    {
        notionals[n_tangent_points] = 0.0;

        is_call[n_tangent_points + 2] = 1;
        strikes[n_tangent_points + 2] =
            strikes[n_tangent_points] -
            payoff_val[n_tangent_points - 1] / deriv_val[n_tangent_points - 1];
        notionals[n_tangent_points + 2] = -deriv_val[n_tangent_points - 1];
    }
    else
    {
        is_call[n_tangent_points + 2]   = 1;
        strikes[n_tangent_points + 2]   = strikes[n_tangent_points];
        notionals[n_tangent_points + 2] = 0.0;
    }

FREE_RETURN:
    if (payoff_val)
        free(payoff_val);
    if (deriv_val)
        free(deriv_val);

    return err;
}

//================================================================
//
//	Get the Strikes for a smooth replication in the case of a monotonic derivative
//
//================================================================
Err SR_find_tangents(
    double start_point,
    double end_point,
    int    index,
    int    nb_points,
    Err (*deriv)(double S, double K, double* out),
    double* final_points,
    double* middle_point,
    int*    left_points,
    int*    right_points)
{
    double left_deriv, right_deriv, middle_deriv, coef, temp_middle;
    int    temp_left_points, temp_right_points;
    Err    err = NULL;

    *middle_point = 0.5 * (start_point + end_point);

    if (nb_points > 1)
    {
        err = deriv(start_point, 0.0, &left_deriv);
        if (err)
            goto FREE_RETURN;
        err = deriv(*middle_point, 0.0, &middle_deriv);
        if (err)
            goto FREE_RETURN;
        err = deriv(end_point, 0.0, &right_deriv);
        if (err)
            goto FREE_RETURN;

        coef = (left_deriv - middle_deriv) / (left_deriv - right_deriv);
        if (coef > 0.9)
            coef = 0.9;
        if (coef < 0.1)
            coef = 0.1;

        if (coef > 0.5)
        {
            *left_points  = (int)(nb_points * coef);
            *right_points = nb_points - *left_points;
        }
        else
        {
            *right_points = (int)(nb_points * (1.0 - coef));
            *left_points  = nb_points - *right_points;
        }

        err = SR_find_tangents(
            start_point,
            *middle_point,
            index,
            *left_points,
            deriv,
            final_points,
            &temp_middle,
            &temp_left_points,
            &temp_right_points);
        if (err)
            goto FREE_RETURN;
        err = SR_find_tangents(
            *middle_point,
            end_point,
            index + *left_points,
            *right_points,
            deriv,
            final_points,
            &temp_middle,
            &temp_left_points,
            &temp_right_points);
        if (err)
            goto FREE_RETURN;
    }
    else
        final_points[index] = *middle_point;

FREE_RETURN:
    return err;
}

//================================================================
//
// Convex Payoff max(1/S - 1/K,0)
//
//================================================================
Err convex_payoff(double S, double K, double* out)
{
    Err err = NULL;

    if (S == 0)
    {
        err = "function undifined";
        goto FREE_RETURN;
    }

    if (fabs(K) < 1.0e-15)
        *out = 1.0 / S;
    else
    {
        *out = 1.0 / S - 1.0 / K;
        if (*out < 0.0)
        {
            *out = 0.0;
        }
    }

FREE_RETURN:
    return err;
}

Err convex_deriv(double S, double K, double* out)
{
    Err err = NULL;

    if (S == 0)
    {
        err = "derivative undifined";
        goto FREE_RETURN;
    }

    if (fabs(K) < 1.0e-13)
        *out = -1.0 / S / S;
    else
    {
        if (S < K + 1.0e-13)
        {
            *out = -1.0 / S / S;
        }
        else
        {
            *out = 0.0;
        }
    }

FREE_RETURN:
    return err;
}

//================================================================
//
//	PRICING OF A CONVEX FORWARD
//
//================================================================
Err Convex_replication(  // Market
    long    Today,
    long    fx_spot_date,
    double  Spot,
    double* fx_mkt_vol,
    long*   fx_mkt_vol_date,
    double* fx_mkt_smile_alpha,
    double* fx_mkt_smile_beta,
    double* fx_mkt_smile_rho,
    double* fx_mkt_smile_pi,
    int     num_fx_mkt_vol,
    int     smile_spec_type,  //	0: lognormal vol + SABR params, 1: sigma-beta + SABR params, 2:
                              //lognormal vol + BMM, 3: lognormal vol + BMM
    char* for_yc,
    char* dom_yc,
    int   is_reversed_fx,
    // Convex Data
    double Notional,
    long   Fix_date,
    long   Settle_date,
    double Kmin,
    double Kmax,
    double Convex_strike,
    // algo data
    int nb_strikes,
    // Output
    double** Replication,
    double*  Price)
{
    double Forward, Fix_time;

    // VOL MKT STRUCTURE
    SMILE_VOL_MARKET smile_mkt        = NULL;
    SMILE_PARAMETERS smile_params     = NULL;
    double*          fwd_at_vol_dates = NULL;

    Err err = NULL;

    //=========================================
    //		CREATE A VOL MARKET STRUCTURE
    //=========================================
    smile_mkt    = calloc(1, sizeof(smile_vol_market));
    smile_params = calloc(1, sizeof(smile_parameters));

    if (!smile_mkt || !smile_params)
    {
        err = "Convex_replication: Memory allocation failure";
        if (err)
            goto FREE_RETURN;
    }

    Forward = Spot * swp_f_df(fx_spot_date, Settle_date, for_yc) /
              swp_f_df(fx_spot_date, Settle_date, dom_yc);
    //=========================================
    // memory allocation
    //=========================================
    err = cpd_alloc_smile_vol_market(num_fx_mkt_vol, smile_mkt);
    if (err)
        goto FREE_RETURN;

    //=========================================
    // fill the mkt structure
    //=========================================
    err = cpd_fill_smile_vol_market(
        Today,
        smile_spec_type,
        num_fx_mkt_vol,
        NULL,
        fx_mkt_vol,
        fx_mkt_smile_alpha,
        fx_mkt_smile_beta,
        fx_mkt_smile_rho,
        fx_mkt_smile_pi,
        NULL,
        fx_mkt_vol_date,
        smile_mkt);

    if (err)
        goto FREE_RETURN;

    //=========================================
    //		GET THE SMILE PARAMS & CALL THE REPLICATION FCT
    //=========================================
    Fix_time = (Fix_date - Today) * YEARS_IN_DAY;

    err = cpd_vol_get_smile_params(0, Fix_date, 0.0, smile_mkt, smile_params);
    if (err)
        goto FREE_RETURN;
    err = cpd_vol_get_linterpVol(0, Fix_date, 0.0, smile_mkt, smile_params);
    if (err)
        goto FREE_RETURN;

    err = Convex_replication_simple(
        Today,
        Forward,
        smile_params->sigma,
        smile_params->alpha,
        smile_params->beta,
        smile_params->rho,
        smile_params->pi,
        smile_params->smile_spec_type,
        is_reversed_fx,
        Notional,
        Fix_date,
        Settle_date,
        Kmin,
        Kmax,
        Convex_strike,
        nb_strikes,
        Replication,
        Price);

    if (err)
        goto FREE_RETURN;

FREE_RETURN:
    if (smile_mkt)
    {
        cpd_free_smile_vol_market(smile_mkt);
        free(smile_mkt);
    }

    if (smile_params)
        free(smile_params);
    return err;
}

Err Convex_replication_simple(  // Market
    long   Today,
    double Forward,
    double fx_mkt_vol,
    double fx_mkt_smile_alpha,
    double fx_mkt_smile_beta,
    double fx_mkt_smile_rho,
    double fx_mkt_smile_pi,
    int    smile_spec_type,  //	0: lognormal vol + SABR params, 1: sigma-beta + SABR params, 2:
                             //lognormal vol + BMM, 3: lognormal vol + BMM
    int is_reversed_fx,
    // Convex Data
    double Notional,
    long   Fix_date,
    long   Settle_date,
    double Kmin,
    double Kmax,
    double Convex_strike,
    // algo data
    int nb_strikes,
    // Output
    double** Replication,
    double*  Price)
{
    double dtemp, Fix_time;
    int    i, itemp1, itemp2, nb_tangents = nb_strikes - 1;

    // Tangents and vols
    double *pStrike = NULL, *pNotional = NULL, *tangent_points = NULL, *pVols = NULL;
    int*    pIsCall = NULL;
    double  dConstant;

    SMILE_PARAMETERS smile_params = NULL;

    Err err = NULL;

    //=========================================
    //		FILL THE SMILE PARAMS
    //=========================================
    smile_params = calloc(1, sizeof(smile_parameters));
    if (!smile_params)
    {
        err = "Convex_replication: Memory allocation failure";
        if (err)
            goto FREE_RETURN;
    }
    smile_params->smile_spec_type = smile_spec_type;
    smile_params->sigma           = fx_mkt_vol;
    smile_params->alpha           = fx_mkt_smile_alpha;
    smile_params->beta            = fx_mkt_smile_beta;
    smile_params->rho             = fx_mkt_smile_rho;
    smile_params->pi              = fx_mkt_smile_pi;

    //=========================================
    //		GET THE REPLICATION TANGENTS
    //=========================================
    tangent_points = calloc(nb_tangents, sizeof(double));

    if (!tangent_points)
    {
        err = "Convex_replication: Memory allocation failure";
        if (err)
            goto FREE_RETURN;
    }

    tangent_points[0] = Kmin;
    if (fabs(Convex_strike) < 1.0e-15)
        tangent_points[nb_tangents - 1] = Kmax;
    else
        tangent_points[nb_tangents - 1] = Convex_strike;

    if (nb_tangents > 2)
    {
        err = SR_find_tangents(
            Kmin,
            fabs(Convex_strike) < 1.0e-15 ? Kmax : Convex_strike,
            1,
            nb_tangents - 2,
            convex_deriv,
            tangent_points,
            &dtemp,
            &itemp1,
            &itemp2);
        if (err)
            goto FREE_RETURN;
    }

    //=========================================
    //		GET THE REPLICATION STRIKES
    //=========================================
    pIsCall =
        (int*)calloc(nb_strikes + 2, sizeof(int));  //+2 for the interpolation on the left & right
    pStrike = (double*)calloc(
        nb_strikes + 2, sizeof(double));  //+2 for the interpolation on the left & right
    pNotional = (double*)calloc(
        nb_strikes + 2, sizeof(double));  //+2 for the interpolation on the left & right
    pVols = (double*)calloc(
        nb_strikes + 2, sizeof(double));  //+2 for the interpolation on the left & right

    if (!pIsCall || !pStrike || !pNotional || !pVols)
    {
        err = "Convex_replication: Memory allocation failure";
        goto FREE_RETURN;
    }

    err = Static_payoff_replication(
        tangent_points,
        nb_tangents,
        Convex_strike,
        0,
        0,
        convex_payoff,
        convex_deriv,
        &dConstant,
        pStrike,
        pNotional,
        pIsCall);

    if (err)
        goto FREE_RETURN;

    //=========================================
    //		GET THE VOL AT THE STRIKES
    //=========================================
    Fix_time = (Fix_date - Today) * YEARS_IN_DAY;

    for (i = 0; i < nb_strikes + 2; i++)
    {
        if (!is_reversed_fx)
        {
            err = cpd_vol_get_vol(
                Forward, Fix_time, pStrike[i], smile_params, SRT_LOGNORMAL, &(pVols[i]));
            if (err)
                goto FREE_RETURN;
        }
        else
        {
            err = cpd_vol_get_vol(
                1 / Forward, Fix_time, 1 / pStrike[i], smile_params, SRT_LOGNORMAL, &(pVols[i]));
            if (err)
                goto FREE_RETURN;
        }
    }

    //=========================================
    //		GET THE PRICE
    //=========================================
    *Price = dConstant;

    for (i = 0; i < nb_strikes + 2; i++)
        *Price += pNotional[i] * srt_f_optblksch(
                                     Forward,
                                     pStrike[i],
                                     pVols[i],
                                     Fix_time,
                                     1.0,
                                     pIsCall[i] ? SRT_CALL : SRT_PUT,
                                     PREMIUM);

    if (Replication)
    {
        for (i = 0; i < nb_strikes + 2; i++)
        {
            Replication[i][0] = pNotional[i];
            Replication[i][1] = pStrike[i];
            Replication[i][2] = pVols[i];
            Replication[i][3] = (double)pIsCall[i];
        }
    }
    Replication[0][4] = dConstant;

FREE_RETURN:
    if (smile_params)
        free(smile_params);
    if (tangent_points)
        free(tangent_points);
    if (pIsCall)
        free(pIsCall);
    if (pStrike)
        free(pStrike);
    if (pNotional)
        free(pNotional);
    if (pVols)
        free(pVols);

    return err;
}