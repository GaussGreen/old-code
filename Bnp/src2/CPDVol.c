/* ==========================================================
        FILENAME :			CPDVol.c

        PURPOSE:			define a vol market

        AUTHOR:				J. Dinh

        DATE:				22-OCT-2004
   ========================================================== */

#include "CPDVol.h"

#include "OTCutils.h"
#include "math.h"
#include "opfnctns.h"
#include "opsabrgenericinterp.h"
#include "srt_h_all.h"

/* =========================================================================
 Macro
========================================================================= */
#define POS_VAL(X) ((X) > 0 ? (X) : 0)

#define CALL_VAL(FWD, STRIKE, D, S) ((FWD)*norm((D) + (S)) - (STRIKE)*norm((D)))

#define PUT_VAL(FWD, STRIKE, D, S) (-(FWD)*norm(-(D) - (S)) + (STRIKE)*norm(-(D)))

#define OPT_VAL_MACRO(TYPE, FWD, STRIKE, STD, HALF_STD)                                         \
    ((TYPE) == 0                                                                                \
         ? 0.0                                                                                  \
         : ((TYPE) == 1                                                                         \
                ? POS_VAL((FWD) - (STRIKE))                                                     \
                : ((TYPE) == 2 ? POS_VAL((STRIKE) - (FWD))                                      \
                               : ((TYPE) == 3 ? CALL_VAL(                                       \
                                                    (FWD),                                      \
                                                    (STRIKE),                                   \
                                                    log((FWD) / (STRIKE)) / (STD) - (HALF_STD), \
                                                    (STD))                                      \
                                              : PUT_VAL(                                        \
                                                    (FWD),                                      \
                                                    (STRIKE),                                   \
                                                    log((FWD) / (STRIKE)) / (STD) - (HALF_STD), \
                                                    (STD))))))

Err cpd_alloc_smile_vol_market(int num_vols, SMILE_VOL_MARKET smile_mkt)
{
    Err err = NULL;

    smile_mkt->forward = calloc(num_vols, sizeof(double));
    smile_mkt->sigma   = calloc(num_vols, sizeof(double));
    smile_mkt->alpha   = calloc(num_vols, sizeof(double));
    smile_mkt->beta    = calloc(num_vols, sizeof(double));
    smile_mkt->rho     = calloc(num_vols, sizeof(double));
    smile_mkt->pi      = calloc(num_vols, sizeof(double));
    smile_mkt->times   = calloc(num_vols, sizeof(double));
    smile_mkt->dates   = calloc(num_vols, sizeof(long));

    if (!smile_mkt->forward || !smile_mkt->sigma || !smile_mkt->alpha || !smile_mkt->beta ||
        !smile_mkt->rho || !smile_mkt->pi || !smile_mkt->times || !smile_mkt->dates)
    {
        err = "cpd_alloc_smile_vol_market: Memory allocation failure";
    }

    return err;
}

Err cpd_fill_smile_vol_market(
    long             today,
    int              smile_spec_type,  // 0 SABR with ATMLOG, 1 SABR with ATMBETA, 2 BMM ...
    int              num_vols,
    double*          forward,
    double*          sigma,
    double*          alpha,
    double*          beta,
    double*          rho,
    double*          pi,
    double*          times,
    long*            dates,
    SMILE_VOL_MARKET smile_mkt)
{
    int i;
    Err err = NULL;

    smile_mkt->num_vols        = num_vols;
    smile_mkt->smile_spec_type = smile_spec_type;

    for (i = 0; i < smile_mkt->num_vols; i++)
    {
        if (times)
            smile_mkt->times[i] = times[i];
        if (dates)
            smile_mkt->dates[i] = dates[i];
        if (forward)
            smile_mkt->forward[i] = forward[i];
        if (sigma)
            smile_mkt->sigma[i] = sigma[i];
        if (alpha)
            smile_mkt->alpha[i] = alpha[i];
        if (beta)
            smile_mkt->beta[i] = beta[i];
        if (rho)
            smile_mkt->rho[i] = rho[i];
        if (pi)
            smile_mkt->pi[i] = pi[i];

        if (!dates && times)
            smile_mkt->dates[i] = today + (long)(smile_mkt->times[i] * DAYS_IN_YEAR);
        if (!times && dates)
            smile_mkt->times[i] = (smile_mkt->dates[i] - today) * YEARS_IN_DAY;

        if (!times && !dates)
        {
            err = "smile_mkt: dates or times have to be passed";
            goto FREE_RETURN;
        }
    }

FREE_RETURN:
    return err;
}

Err cpd_check_smile_vol_market(SMILE_VOL_MARKET smile_mkt)
{
    int i, IsTime = 0;
    Err err = NULL;

    if (smile_mkt->times[smile_mkt->num_vols - 1] > 0.0)
        IsTime++;

    if (IsTime && smile_mkt->dates[smile_mkt->num_vols - 1] > 0.0)
        IsTime++;

    // times or dates have to be non zero
    for (i = 0; i < smile_mkt->num_vols - 1; i++)
    {
        if (smile_mkt->dates[i] == 0 && fabs(smile_mkt->times[i]) < 1e-15)
        {
            err = "cpd_check_smile_vol_market: dates or times have to be non zero";
            goto FREE_RETURN;
        }
    }

    // times have to be increasing
    if (IsTime > 0)
    {
        for (i = 0; i < smile_mkt->num_vols - 2; i++)
        {
            if (smile_mkt->times[i + 1] <= smile_mkt->times[i])
            {
                err = "cpd_check_smile_vol_market: times have to be increasing";
                goto FREE_RETURN;
            }
        }
    }

    // dates have to be increasing
    if (IsTime != 1)
    {
        for (i = 0; i < smile_mkt->num_vols; i++)
        {
            if (smile_mkt->dates[i + 1] <= smile_mkt->dates[i])
            {
                err = "cpd_check_smile_vol_market: dates have to be increasing";
                goto FREE_RETURN;
            }
        }
    }

FREE_RETURN:
    return err;
}

Err cpd_free_smile_vol_market(SMILE_VOL_MARKET smile_mkt)
{
    Err err = NULL;

    if (smile_mkt->forward)
        free(smile_mkt->forward);
    if (smile_mkt->sigma)
        free(smile_mkt->sigma);
    if (smile_mkt->alpha)
        free(smile_mkt->alpha);
    if (smile_mkt->beta)
        free(smile_mkt->beta);
    if (smile_mkt->rho)
        free(smile_mkt->rho);
    if (smile_mkt->pi)
        free(smile_mkt->pi);
    if (smile_mkt->times)
        free(smile_mkt->times);
    if (smile_mkt->dates)
        free(smile_mkt->dates);

    return err;
}

Err cpd_vol_get_vol(
    double           forward,
    double           fix_time,
    double           strike,
    SMILE_PARAMETERS smile_params,
    SrtDiffusionType output_vol,
    double*          smile_std)
{
    Err err = NULL;

    if (smile_params->smile_spec_type < 2)
    {
        // SABR
        err = srt_f_optsarbvol(
            forward,
            strike,
            fix_time,
            smile_params->sigma,
            smile_params->alpha,
            smile_params->beta,
            smile_params->rho,
            smile_params->smile_spec_type ? SRT_BETAVOL : SRT_LOGNORMAL,
            output_vol,
            smile_std);
    }
    else if (smile_params->smile_spec_type == 2 || smile_params->smile_spec_type == 3)
    {
        // BMM
        err = srt_f_optbmmvol(
            forward,
            strike,
            fix_time,
            smile_params->sigma,
            smile_params->alpha,
            smile_params->beta,
            smile_params->rho,
            smile_params->pi,
            smile_params->smile_spec_type == 3 ? SRT_BETAVOL : SRT_LOGNORMAL,
            output_vol,
            smile_std);
    }
    else
    {
        err = "smile_spec_type has to be less or equal to 3)";
    }

    return err;
}

Err cpd_vol_get_price(
    int              type,
    double           Forward,
    double           Maturity,
    double           Strike,
    SMILE_PARAMETERS smile_params,
    double*          Value)
{
    double smile_std = 0.0, smile_half_std = 0.0;
    Err    err = NULL;

    if (smile_params->smile_spec_type < 3)  // SABR Log, SABR Beta or BMM Log
    {
        if (type >= 3)
        {
            err =
                cpd_vol_get_vol(Forward, Maturity, Strike, smile_params, SRT_LOGNORMAL, &smile_std);
            if (err)
                goto FREE_RETURN;

            smile_std *= sqrt(Maturity);
            smile_half_std = 0.5 * smile_std;
        }

        if ((Forward < 0.0 && type >= 3) || (smile_std < 1.0e-6 && type >= 3))
            type -= 2;

        *Value = OPT_VAL_MACRO(type, Forward, Strike, smile_std, smile_half_std);
    }
    else if (smile_params->smile_spec_type == 3)  // BMM beta
    {
        if (type >= 3)
        {
            if (!smile_params->Use_BetaQuick_in_MC)
            {
                err = BMM_Option_Price(
                    Forward,
                    Strike,
                    Maturity,
                    smile_params->sigma,
                    smile_params->alpha,
                    smile_params->beta,
                    smile_params->rho,
                    smile_params->pi,
                    type == 3 ? SRT_CALL : SRT_PUT,
                    Value);
                if (err)
                    goto FREE_RETURN;
            }
            else
            {
                err = BMM_Option_Price_quick(
                    Forward,
                    Strike,
                    Maturity,
                    smile_params->sigma,
                    smile_params->alpha,
                    smile_params->beta,
                    smile_params->rho,
                    smile_params->pi,
                    type == 3 ? SRT_CALL : SRT_PUT,
                    Value);

                if (err)
                    goto FREE_RETURN;
            }
        }
        else
        {
            *Value = OPT_VAL_MACRO(type, Forward, Strike, smile_std, smile_half_std);
        }
    }

FREE_RETURN:
    return err;
}

Err cpd_vol_get_smile_params(
    int              IsTime,  // 0 will use date, 1 will use time
    long             date,
    double           time,
    SMILE_VOL_MARKET smile_mkt,
    SMILE_PARAMETERS smile_params)
{
    int     j;
    double* dates = NULL;

    Err err = NULL;

    //====================================
    // LINEAR INTERP OF THE PARAMETERS
    //====================================
    smile_params->smile_spec_type = smile_mkt->smile_spec_type;

    if (IsTime)
    {
        smile_params->alpha =
            interp(smile_mkt->times, smile_mkt->alpha, smile_mkt->num_vols, time, 0, NULL);
        smile_params->beta =
            interp(smile_mkt->times, smile_mkt->beta, smile_mkt->num_vols, time, 0, NULL);
        smile_params->rho =
            interp(smile_mkt->times, smile_mkt->rho, smile_mkt->num_vols, time, 0, NULL);

        if (smile_mkt->smile_spec_type == 2 || smile_mkt->smile_spec_type == 3)
            smile_params->pi =
                interp(smile_mkt->times, smile_mkt->pi, smile_mkt->num_vols, time, 0, NULL);
    }
    else
    {
        dates = calloc(smile_mkt->num_vols, sizeof(double));
        for (j = 0; j < smile_mkt->num_vols; j++)
            dates[j] = (double)smile_mkt->dates[j];

        smile_params->alpha = interp(dates, smile_mkt->alpha, smile_mkt->num_vols, date, 0, NULL);
        smile_params->beta  = interp(dates, smile_mkt->beta, smile_mkt->num_vols, date, 0, NULL);
        smile_params->rho   = interp(dates, smile_mkt->rho, smile_mkt->num_vols, date, 0, NULL);

        if (smile_mkt->smile_spec_type == 2 || smile_mkt->smile_spec_type == 3)
            smile_params->pi = interp(dates, smile_mkt->pi, smile_mkt->num_vols, date, 0, NULL);
    }

    if (dates)
        free(dates);

    return err;
}

Err cpd_vol_get_linterpVol(
    int              IsTime,  // 0 will use date, 1 will use time
    long             date,
    double           time,
    SMILE_VOL_MARKET smile_mkt,
    SMILE_PARAMETERS smile_params)
{
    int     j;
    double* dates = NULL;

    Err err = NULL;

    //====================================
    // LINEAR INTERP OF THE PARAMETERS
    //====================================
    smile_params->smile_spec_type = smile_mkt->smile_spec_type;

    if (IsTime)
        smile_params->sigma =
            interp(smile_mkt->times, smile_mkt->sigma, smile_mkt->num_vols, time, 0, NULL);
    else
    {
        dates = calloc(smile_mkt->num_vols, sizeof(double));
        for (j = 0; j < smile_mkt->num_vols; j++)
            dates[j] = (double)smile_mkt->dates[j];

        smile_params->sigma = interp(dates, smile_mkt->sigma, smile_mkt->num_vols, date, 0, NULL);
    }

    if (dates)
        free(dates);

    return err;
}

Err cpd_vol_get_smile_params_from_otc(
    int              smile_spec_type,
    int              row_idx,
    int              col_idx,
    OTCPDPRECALC     otc_precalc,
    SMILE_PARAMETERS smile_params)
{
    Err err = NULL;

    smile_params->smile_spec_type = smile_spec_type;

    smile_params->sigma = otc_precalc->pd_fwdsmile_vol[row_idx][col_idx];
    smile_params->alpha = otc_precalc->pd_fwdsmile_alpha[row_idx][col_idx];
    smile_params->beta  = otc_precalc->pd_fwdsmile_beta[row_idx][col_idx];
    smile_params->rho   = otc_precalc->pd_fwdsmile_rho[row_idx][col_idx];

    if (smile_params->smile_spec_type == 2 || smile_params->smile_spec_type == 3)
        smile_params->pi = otc_precalc->pd_fwdsmile_pi[row_idx][col_idx];  // To do

    return err;
}
