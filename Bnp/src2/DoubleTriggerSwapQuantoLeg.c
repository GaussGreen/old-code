/* ===================================================================================
   FILENAME:      DoubleTriggerSwapQuantoLeg.c / LY J-M

   PURPOSE:       Computes the PV of a DoubleTriggerSwapQuanto
   =================================================================================== */

#pragma warning(disable : 4786)  // Disable long name warnings

#include "DoubleTriggerSwapQuantoLeg.h"

#include "DigitalQuantoFloatLeg.h"
#include "math.h"
#include "num_h_allhdr.h"
#include "opHeston.h"
#include "opfnctns.h"
#include "swp_h_all.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "swp_h_vol.h"

// Prices a double put with payoff E[max(strikeX - X_fixed_at_timeX,0) * max(strikeY -
// Y_fixed_at_timeY,0) where X and Y are both lognormals with fwds fwdX and fwdY and vols sigmaX and
// sigmaY and instantaneous correlation rhoXY X is fixed at timeX and Y is fixed at time Y.
Err double_put_lognormal(
    double  x0,
    double  B,
    double  sigma_x,
    double  T_x,
    double  y0,
    double  K,
    double  sigma_y,
    double  T_y,
    double  correl,
    double* result)
{
    Err    err;
    double e_x_1, e_y_1, e_x_2, e_y_2, e_x_3, e_y_3, e_x_4, e_y_4;
    double bivar1, bivar2, bivar3, bivar4, new_correl;

    if (T_x <= 0.0 || T_y <= 0.0 || B <= 0.0 || K <= 0.0 || x0 <= 0.0 || y0 <= 0.0)
    {
        err = "double_put_lognormal: ERROR (maturities <= 0.0 OR strikes <= 0.0 OR fwds <= 0.0)";
        goto FREE_RETURN;
    }

    new_correl = correl * min(T_x, T_y) / (sqrt(T_x) * sqrt(T_y));
    e_x_1      = 1 / (sigma_x * sqrt(T_x)) * (log(B / x0) + 0.5 * sigma_x * sigma_x * T_x);
    e_y_1      = 1 / (sigma_y * sqrt(T_y)) * (log(K / y0) + 0.5 * sigma_y * sigma_y * T_y);
    e_x_2 =
        1 / (sigma_x * sqrt(T_x)) *
        (log(B / x0) - correl * sigma_x * sigma_y * min(T_x, T_y) + 0.5 * sigma_x * sigma_x * T_x);
    e_y_2 = 1 / (sigma_y * sqrt(T_y)) * (log(K / y0) - 0.5 * sigma_y * sigma_y * T_y);
    e_x_3 = 1 / (sigma_x * sqrt(T_x)) * (log(B / x0) - 0.5 * sigma_x * sigma_x * T_x);
    e_y_3 =
        1 / (sigma_y * sqrt(T_y)) *
        (log(K / y0) - correl * sigma_x * sigma_y * min(T_x, T_y) + 0.5 * sigma_y * sigma_y * T_y);
    e_x_4 =
        1 / (sigma_x * sqrt(T_x)) *
        (log(B / x0) - correl * sigma_x * sigma_y * min(T_x, T_y) - 0.5 * sigma_x * sigma_x * T_x);
    e_y_4 =
        1 / (sigma_y * sqrt(T_y)) *
        (log(K / y0) - correl * sigma_x * sigma_y * min(T_x, T_y) - 0.5 * sigma_y * sigma_y * T_y);

    bivar1 = bivar(e_x_1, e_y_1, new_correl);
    bivar2 = bivar(e_x_2, e_y_2, new_correl);
    bivar3 = bivar(e_x_3, e_y_3, new_correl);
    bivar4 = bivar(e_x_4, e_y_4, new_correl);

    *result = B * K * bivar1 - B * y0 * bivar2 - K * x0 * bivar3 +
              x0 * y0 * exp(correl * sigma_x * sigma_y * min(T_x, T_y)) * bivar4;

FREE_RETURN:
    return NULL;
}

// Prices a double inf barrier with payoff E[1{X<B}1{Y<K}]
// where X and Y are both lognormals with fwds fwdX and fwdY
// and vols sigmaX and sigmaY and instantaneous correlation rhoXY
// X is fixed at timeX and Y is fixed at time Y.
// The call spreads of X are B1 < B2 and of Y are K1 < K2
Err double_inf_barrier_lognormal(
    double  fwdX,
    double  B1,
    double  volX_B1,
    double  B2,
    double  volX_B2,
    double  timeX,
    double  fwdY,
    double  K1,
    double  volY_K1,
    double  K2,
    double  volY_K2,
    double  timeY,
    double  rhoXY,
    double* result)
{
    Err    err = NULL;
    double c1, c2, c3, c4, dX, dY;

    dX = B2 - B1;
    dY = K2 - K1;

    if (dX < 0 || dY < 0)
    {
        err = "Strikes in the wrong order B1<B2 & K1<K2";
        goto FREE_RETURN;
    }

    err = double_put_lognormal(fwdX, B2, volX_B2, timeX, fwdY, K2, volY_K2, timeY, rhoXY, &c1);
    if (err)
        goto FREE_RETURN;

    err = double_put_lognormal(fwdX, B2, volX_B2, timeX, fwdY, K1, volY_K1, timeY, rhoXY, &c2);
    if (err)
        goto FREE_RETURN;

    err = double_put_lognormal(fwdX, B1, volX_B1, timeX, fwdY, K2, volY_K2, timeY, rhoXY, &c3);
    if (err)
        goto FREE_RETURN;

    err = double_put_lognormal(fwdX, B1, volX_B1, timeX, fwdY, K1, volY_K1, timeY, rhoXY, &c4);
    if (err)
        goto FREE_RETURN;

    *result = (c1 - c2 - c3 + c4) / (dX * dY);

FREE_RETURN:
    return err;
}

/*
//From frequency "A" and "BB" get the numerical compounding 1 and the dRate
Err   frequency_and_basis_to_compounding_and_dRateConv(char* frequency,
                                                                                                           char* basis,
                                                                                                           SrtCompounding* compounding,
                                                                                                           double* dRateConv)
{
        //Cf CTSProdStruct ref_conv
        Err err;
        SrtBasisCode srt_basis_code;

        err = interp_compounding(frequency, compounding);
        if(err) goto FREE_RETURN;

        err = interp_basis(basis, &srt_basis_code);
        if(err) goto FREE_RETURN;

        if(srt_basis_code == BASIS_ACT_360)
                *dRateConv = 365.0 / 360.0;
        else
                *dRateConv = 1.0;

FREE_RETURN:
        return err;
}



//Calculates the CMS end date using the tenor and the start date
Err   cms_get_tenor_num_and_unit_from_tenor(char* cms_tenor,
                                                                                        int*
tenor_num, char* tenor_unit)
{
        int string_len;
        char tenor[256]; // "5Y"

        strcpy(tenor, cms_tenor);
        string_len = strlen(tenor);
        *tenor_unit = tenor[string_len - 1];
        tenor[string_len - 1] = '\0';
        *tenor_num = atoi(tenor);

        return NULL;
}


//Gets the lognormal vol (to be modified for MAD compatibility by using the MAD function
//used in the swp_f_truncvol function used in swp_f_cms_option
Err digital_quanto_float_leg_get_lognormal_vol(char* vc_name,
                                                                                           double
fwd, double fixing_time, long start_date, long end_date, double strike, double* result_volatility)
{
        Err err;
        double vol, power, temp_vol;

        err = swp_f_vol(vc_name,
                                        start_date,
                                        end_date,
                                        strike,
                                        &vol,
                                        &power);
        if (err) goto FREE_RETURN;

        //Conversion from normal to lognormal
        if (power == 0.0)
        {
                err = srt_f_optsarbvol(fwd,
                                                           strike,
                                                           fixing_time,
                                                           vol,
                                                           0,
                                                           0,
                                                           0,
                                                           SRT_BETAVOL,
                                                           SRT_LOGNORMAL,
                                                           &temp_vol);

                if (err) goto FREE_RETURN;

                vol = temp_vol;
        }

        FREE_RETURN:
        *result_volatility = vol;
        return err;
}
*/

// Calculates a CMS rate quanto or non quanto
Err double_trigger_swap_quanto_cms_rate_quanto(
    double fwd,
    double cms_fixing_time,
    double cms_start_time,
    double cms_nfp,

    SrtCompounding cms_compounding,
    double         cms_delay,
    double         cms_dRateConv,

    SrtDiffusionType cms_vol_type,
    int              cms_flat_vol_or_smile_approx,
    long             cms_start_date,
    long             cms_theo_end_date,
    int              cms_cash_vol,
    double           fwd_spread,

    char*   cms_vc_name,
    int     n_cms_strikes_in_vol,
    double* cms_strikes_in_vol,
    int     cms_dom_for,

    double* cms_quanto_corr_times,
    double* cms_quanto_corr_values,
    int     n_cms_quanto_corr,

    double* fwd_fx_vols_times,
    double* fwd_fx_vols_values,
    int     n_fwd_fx_vols,

    int use_CMSRate_CMSTEC,

    double* result_cms_rate_quanto)
{
    Err err = NULL;
    // double cms, cms_atm_log_vol, cms_quanto_correlation, fwd_fx_vol, flat_vol, temp_vol;
    double cms = 0.0;

    if (use_CMSRate_CMSTEC == 0)
    {
        // Call the old functions to calculate the CMS Quanto
        digital_quanto_float_leg_cms_rate_quanto(
            fwd,
            cms_fixing_time,
            cms_start_time,
            cms_nfp,

            cms_compounding,
            cms_delay,
            cms_dRateConv,

            cms_vol_type,
            cms_flat_vol_or_smile_approx,
            cms_start_date,
            cms_theo_end_date,
            cms_cash_vol,
            fwd_spread,

            cms_vc_name,
            n_cms_strikes_in_vol,
            cms_strikes_in_vol,
            cms_dom_for,

            cms_quanto_corr_times,
            cms_quanto_corr_values,
            n_cms_quanto_corr,

            fwd_fx_vols_times,
            fwd_fx_vols_values,
            n_fwd_fx_vols,

            &cms);
    }
    else
    {
        err = "CMSTEC not yet implemented in this code";
        goto FREE_RETURN;
    }

FREE_RETURN:
    *result_cms_rate_quanto = cms;
    return err;
}

// Calculates the PV of a double trigger swap quanto leg
Err double_trigger_swap_quanto_leg(
    char*            dom_yc_name,
    char*            dom_vc_name,
    int              dom_spot_lag,
    double*          dom_strikes_in_vol,
    int              n_dom_strikes_in_vol,
    SrtDiffusionType dom_vol_type,
    int              dom_cash_vol,

    char*            for_yc_name,
    char*            for_vc_name,
    int              for_spot_lag,
    double*          for_strikes_in_vol,
    int              n_for_strikes_in_vol,
    SrtDiffusionType for_vol_type,
    int              for_cash_vol,

    long today,

    char* cms1_tenor,
    char* cms1_frequency,
    char* cms1_basis,
    char* cms1_refrate,
    int   cms1_dom_for,

    char* cms2_tenor,
    char* cms2_frequency,
    char* cms2_basis,
    char* cms2_refrate,
    int   cms2_dom_for,

    long*   exo_start_dates,
    long*   exo_end_dates,
    long*   exo_pay_dates,
    char**  exo_basis,
    double* exo_notionals,
    long*   exo_cms1_fixing_dates,
    double* exo_cms1_barriers,
    long*   exo_cms2_fixing_dates,
    double* exo_cms2_barriers,
    double* exo_m12,
    double* exo_g2,
    double* exo_m2,
    double* exo_g3,
    double* exo_m3,
    double* exo_g4_X,
    double* exo_g4_Y,
    double* exo_m4,
    double* exo_cms1_past_fixings,
    double* exo_cms2_past_fixings,
    int     n_exo_coupons,

    double* cms1_cms2_corr_times,
    double* cms1_cms2_corr_values_bid,
    double* cms1_cms2_corr_values_ask,
    int     n_cms1_cms2_corr,

    double* cms1_quanto_corr_times,
    double* cms1_quanto_corr_values,
    int     n_cms1_quanto_corr,

    double* cms2_quanto_corr_times,
    double* cms2_quanto_corr_values,
    int     n_cms2_quanto_corr,

    double* fwd_fx_vols_times,
    double* fwd_fx_vols_values,
    int     n_fwd_fx_vols,

    int    pay_rec,
    double call_spread,
    int    use_CMSRate_CMSTEC,
    int    eod_fix_flag,
    int    eod_pay_flag,

    double* cpn1,
    double* cpn2,
    double* cpn3,
    double* cpn4,
    double* DoubleTriggerSwapQuantoLeg)
{
    Err            err = NULL;
    double         fwd1, fwd2, cms1, cms2, fwd_spread1, fwd_spread2;
    double         cms1_fixing_time, cms1_delay, cms1_start_time, cms1_nfp, cms1_dRateConv;
    double         cms2_fixing_time, cms2_delay, cms2_start_time, cms2_nfp, cms2_dRateConv;
    long           cms1_start_date, cms2_start_date, cms1_theo_end_date, cms2_theo_end_date;
    char           cms1_tenor_unit, cms2_tenor_unit;
    SrtCompounding cms1_compounding, cms2_compounding;
    SrtBasisCode   srt_exo_basis_code;
    double cms1_ats_log_vol, cms1_B1_log_vol, cms1_B2_log_vol, cms2_ats_log_vol, cms2_K1_log_vol,
        cms2_K2_log_vol;
    double cms1_numeraire_adjustment = 1.0, cms2_numeraire_adjustment = 1.0, cms1_cms2_correlation;
    double put1, put2, B1, B2, K1, K2, df_pay, cvg, double_inf_barrier;
    int    i, cms1_tenor_num, cms2_tenor_num;

    // CMS1 market details
    char cms1_yc_name[256];
    char cms1_vc_name[256];
    int  cms1_spot_lag;

    double*          cms1_strikes_in_vol = NULL;
    int              n_cms1_strikes_in_vol;
    SrtDiffusionType cms1_vol_type;
    int              cms1_cash_vol;

    // CMS2 market details
    char cms2_yc_name[256];
    char cms2_vc_name[256];
    int  cms2_spot_lag;

    double*          cms2_strikes_in_vol = NULL;
    int              n_cms2_strikes_in_vol;
    SrtDiffusionType cms2_vol_type;
    int              cms2_cash_vol;

    // Results
    double exotic_leg1_PV = 0.0, exotic_leg2_PV = 0.0, exotic_leg3_PV = 0.0, exotic_leg4_PV = 0.0;
    double exotic_coupon1_PV = 0.0, exotic_coupon2_PV = 0.0, exotic_coupon3_PV = 0.0,
           exotic_coupon4_PV = 0.0;

    // Checks that fixing & pay dates are increasing
    i = 0;
    while (i < n_exo_coupons - 1)
    {
        if ((exo_cms1_fixing_dates[i] >= exo_cms1_fixing_dates[i + 1]) ||
            (exo_cms2_fixing_dates[i] >= exo_cms2_fixing_dates[i + 1]))
        {
            err = "cms_fixing_dates are not in increasing order";
            goto FREE_RETURN;
        }
        i++;
    }

    // Checks that fixing dates <= pay dates
    i = 0;
    while (i < n_exo_coupons)
    {
        if ((exo_cms1_fixing_dates[i] > exo_pay_dates[i]) ||
            (exo_cms2_fixing_dates[i] > exo_pay_dates[i]))
        {
            err = "cms is paid before the before the fixing date";
            goto FREE_RETURN;
        }
        i++;
    }

    // Reading CMS1 market parameters
    if (cms1_dom_for == 0)
    {
        // CMS1 is domestic
        strcpy(cms1_yc_name, dom_yc_name);
        strcpy(cms1_vc_name, dom_vc_name);
        cms1_spot_lag = dom_spot_lag;

        cms1_strikes_in_vol   = dom_strikes_in_vol;
        n_cms1_strikes_in_vol = n_dom_strikes_in_vol;
        cms1_vol_type         = dom_vol_type;
        cms1_cash_vol         = dom_cash_vol;
    }
    else
    {
        // CMS1 is foreign
        strcpy(cms1_yc_name, for_yc_name);
        strcpy(cms1_vc_name, for_vc_name);
        cms1_spot_lag = for_spot_lag;

        cms1_strikes_in_vol   = for_strikes_in_vol;
        n_cms1_strikes_in_vol = n_for_strikes_in_vol;
        cms1_vol_type         = for_vol_type;
        cms1_cash_vol         = for_cash_vol;
    }

    // Reading CMS2 market parameters
    if (cms2_dom_for == 0)
    {
        // CMS2 is domestic
        strcpy(cms2_yc_name, dom_yc_name);
        strcpy(cms2_vc_name, dom_vc_name);
        cms2_spot_lag = dom_spot_lag;

        cms2_strikes_in_vol   = dom_strikes_in_vol;
        n_cms2_strikes_in_vol = n_dom_strikes_in_vol;
        cms2_vol_type         = dom_vol_type;
        cms2_cash_vol         = dom_cash_vol;
    }
    else
    {
        // CMS2 is foreign
        strcpy(cms2_yc_name, for_yc_name);
        strcpy(cms2_vc_name, for_vc_name);
        cms2_spot_lag = for_spot_lag;

        cms2_strikes_in_vol   = for_strikes_in_vol;
        n_cms2_strikes_in_vol = n_for_strikes_in_vol;
        cms2_vol_type         = for_vol_type;
        cms2_cash_vol         = for_cash_vol;
    }

    // Computes the compounding and the day rate convention
    frequency_and_basis_to_compounding_and_dRateConv(
        cms1_frequency, cms1_basis, &cms1_compounding, &cms1_dRateConv);

    frequency_and_basis_to_compounding_and_dRateConv(
        cms2_frequency, cms2_basis, &cms2_compounding, &cms2_dRateConv);

    // Parsing the tenor num and tenor unit of the CMS assets
    cms_get_tenor_num_and_unit_from_tenor(cms1_tenor, &cms1_tenor_num, &cms1_tenor_unit);

    cms_get_tenor_num_and_unit_from_tenor(cms2_tenor, &cms2_tenor_num, &cms2_tenor_unit);

    // Initializes the coupon index
    i = 0;

    // Remove all periods which have already been paid
    while (exo_pay_dates[i] < today + eod_pay_flag)
        i++;

    // Pricing of coupons where at least one CMS has been fixed TO BE CHANGED
    while (min(exo_cms1_fixing_dates[i], exo_cms2_fixing_dates[i]) < today + eod_fix_flag)
    {
        // Both fixed
        if ((exo_cms1_fixing_dates[i] < today + eod_fix_flag) &&
            (exo_cms2_fixing_dates[i] < today + eod_fix_flag))
        {
            if ((exo_cms1_past_fixings[i] < exo_cms1_barriers[i]) &&
                (exo_cms2_past_fixings[i] < exo_cms2_barriers[i]))
            {
                // If CMS1<B && CMS2<K then m12
                exotic_coupon1_PV = exo_m12[i];
            }
            else if (
                (exo_cms1_past_fixings[i] >= exo_cms1_barriers[i]) &&
                (exo_cms2_past_fixings[i] < exo_cms2_barriers[i]))
            {
                // If CMS1>=B && CMS2<K then g2 * CMS1 + m2
                exotic_coupon2_PV = exo_g2[i] * exo_cms1_past_fixings[i] + exo_m2[i];
            }
            else if (
                (exo_cms1_past_fixings[i] < exo_cms1_barriers[i]) &&
                (exo_cms2_past_fixings[i] >= exo_cms2_barriers[i]))
            {
                // If CMS1<B && CMS2>=K then g3 * CMS2 + m3
                exotic_coupon3_PV = exo_g3[i] * exo_cms2_past_fixings[i] + exo_m3[i];
            }
            else if (
                (exo_cms1_past_fixings[i] >= exo_cms1_barriers[i]) &&
                (exo_cms2_past_fixings[i] >= exo_cms2_barriers[i]))
            {
                // If CMS1<B && CMS2>=K then 0.5 * (g4_X * CMS1 + g4_Y * CMS2 + m4)
                exotic_coupon4_PV = 0.5 * (exo_g4_X[i] * exo_cms1_past_fixings[i] +
                                           exo_g4_Y[i] * exo_cms2_past_fixings[i] + exo_m4[i]);
            }
        }
        else
        {
            err =
                "Handling of past fixings where one of the CMS is fixed in advance and the other "
                "one in arrears is not yet available";
            goto FREE_RETURN;
        }

        /*
                //One asset fixed in the case of a Digital Quanto Float TO BE REWRITTEN if
           DoubleTipTop with fixing in advance and in arrears else
                {
                        //CMS1 has fixed
                        if (exo_cms1_fixing_dates[i] < today + eod_fix_flag)
                        {
                                //Payoff is a forward CMS2 + margin (quanto corrected if needed)
                                if (((exo_cms1_past_fixings[i] > exo_barriers[i]) && (above_below ==
           0)) ||
                                        ((exo_cms1_past_fixings[i] < exo_barriers[i])  &&
           (above_below == 1)) ||
                                        ((exo_cms1_past_fixings[i] == exo_barriers[i]) &&
           (accrue_on_barrier == 1)))
                                {
                                        //Calculates the start and theo end date of the CMS2
           underlying cms2_start_date = add_unit(exo_cms2_fixing_dates[i],cms2_spot_lag, SRT_BDAY,
           MODIFIED_SUCCEEDING); if(cms2_tenor_unit == 'M')
                                        {
                                                cms2_theo_end_date = add_unit(cms2_start_date,
           cms2_tenor_num, SRT_MONTH, NO_BUSDAY_CONVENTION);
                                        }
                                        else
                                        {
                                                cms2_theo_end_date = add_unit(cms2_start_date,
           cms2_tenor_num, SRT_YEAR, NO_BUSDAY_CONVENTION);
                                        }

                                        //Calculates number of full periods for the underlying of
           CMS2 cms2_fixing_time = (exo_cms2_fixing_dates[i] - today) * YEARS_IN_DAY;
                                        cms2_start_time = (cms2_start_date - today) * YEARS_IN_DAY;
                                        cms2_delay = (exo_pay_dates[i] - cms2_start_date) *
           YEARS_IN_DAY; cms2_nfp = ((int) ((cms2_theo_end_date - cms2_start_date) * YEARS_IN_DAY *
           cms2_compounding + 0.1)) * 1.0; fwd_spread2 = swp_f_spread(cms2_start_date,
           cms2_theo_end_date, cms2_refrate);

                                        //Calculates the forward
                                        err = swp_f_ForwardRate(cms2_start_date,
                                                                                        cms2_theo_end_date,
                                                                                        cms2_frequency,
                                                                                        cms2_basis,
                                                                                        cms2_yc_name,
                                                                                        cms2_refrate,
                                                                                        &fwd2);
                                        if(err) goto FREE_RETURN;

                                        //Calculates the CMS Rate
                                        err = digital_quanto_float_leg_cms_rate_quanto(fwd2,
                                                                                                                                   cms2_fixing_time,
                                                                                                                                   cms2_start_time,
                                                                                                                                   cms2_nfp,
                                                                                                                                   cms2_compounding,
                                                                                                                                   cms2_delay,
                                                                                                                                   cms2_dRateConv,
                                                                                                                                   cms2_vol_type,
                                                                                                                                   cms_flat_vol_or_smile_approx,
                                                                                                                                   cms2_start_date,
                                                                                                                                   cms2_theo_end_date,
                                                                                                                                   cms2_cash_vol,
                                                                                                                                   fwd_spread2,
                                                                                                                                   cms2_vc_name,
                                                                                                                                   n_cms2_strikes_in_vol,
                                                                                                                                   cms2_strikes_in_vol,
                                                                                                                                   cms2_dom_for,
                                                                                                                                   cms2_quanto_corr_times,
                                                                                                                                   cms2_quanto_corr_values,
                                                                                                                                   n_cms2_quanto_corr,
                                                                                                                                   fwd_fx_vols_times,
                                                                                                                                   fwd_fx_vols_values,
                                                                                                                                   n_fwd_fx_vols,
                                                                                                                                   &cms2);
                                        if(err) goto FREE_RETURN;
                                        exotic_coupon_PV = exo_gearings[i] * cms2 + exo_margins[i];
                                }
                                else
                                        exotic_coupon_PV = 0.0;
                        }


                         //CMS2 has fixed we need to price a digital on CMS1
                        else
                        {
                                //Calculates the start and theo end date of the CMS1 underlying
                                cms1_start_date = add_unit(exo_cms1_fixing_dates[i],cms1_spot_lag,
           SRT_BDAY, MODIFIED_SUCCEEDING); if(cms1_tenor_unit == 'M')
                                {
                                        cms1_theo_end_date = add_unit(cms1_start_date,
           cms1_tenor_num, SRT_MONTH, NO_BUSDAY_CONVENTION);
                                }
                                else
                                {
                                        cms1_theo_end_date = add_unit(cms1_start_date,
           cms1_tenor_num, SRT_YEAR, NO_BUSDAY_CONVENTION);
                                }

                                //Calculates number of full periods for the underlying of CMS1
                                cms1_fixing_time = (exo_cms1_fixing_dates[i] - today) *
           YEARS_IN_DAY; cms1_start_time = (cms1_start_date - today) * YEARS_IN_DAY; cms1_delay =
           (exo_pay_dates[i] - cms1_start_date) * YEARS_IN_DAY; cms1_nfp = ((int)
           ((cms1_theo_end_date - cms1_start_date) * YEARS_IN_DAY * cms1_compounding + 0.1)) * 1.0;
                                fwd_spread1 = swp_f_spread(cms1_start_date, cms1_theo_end_date,
           cms1_refrate);

                                //Calculates the forward
                                err = swp_f_ForwardRate(cms1_start_date,
                                                                                cms1_theo_end_date,
                                                                                cms1_frequency,
                                                                                cms1_basis,
                                                                                cms1_yc_name,
                                                                                cms1_refrate,
                                                                                &fwd1);
                                if(err) goto FREE_RETURN;

                                //Calculates the CMS Rate
                                err = digital_quanto_float_leg_cms_rate_quanto(fwd1,
                                                                                                                           cms1_fixing_time,
                                                                                                                           cms1_start_time,
                                                                                                                           cms1_nfp,
                                                                                                                           cms1_compounding,
                                                                                                                           cms1_delay,
                                                                                                                           cms1_dRateConv,
                                                                                                                           cms1_vol_type,
                                                                                                                           cms_flat_vol_or_smile_approx,
                                                                                                                           cms1_start_date,
                                                                                                                           cms1_theo_end_date,
                                                                                                                           cms1_cash_vol,
                                                                                                                           fwd_spread1,
                                                                                                                           cms1_vc_name,
                                                                                                                           n_cms1_strikes_in_vol,
                                                                                                                           cms1_strikes_in_vol,
                                                                                                                           cms1_dom_for,
                                                                                                                           cms1_quanto_corr_times,
                                                                                                                           cms1_quanto_corr_values,
                                                                                                                           n_cms1_quanto_corr,
                                                                                                                           fwd_fx_vols_times,
                                                                                                                           fwd_fx_vols_values,
                                                                                                                           n_fwd_fx_vols,
                                                                                                                           &cms1);
                                if (err) goto FREE_RETURN;


                                //CMS1 ATS B log vol
                                err = digital_quanto_float_leg_get_lognormal_vol(cms1_vc_name,
                                                                                                                                 fwd1,
                                                                                                                                 cms1_fixing_time,
                                                                                                                                 cms1_start_date,
                                                                                                                                 cms1_theo_end_date,
                                                                                                                                 exo_barriers[i],
                                                                                                                                 &cms1_ats_log_vol);
                                if (err) goto FREE_RETURN;

                                put_margin_B = srt_f_optblksch(cms1,
                                                                                           exo_barriers[i],
                                                                                           cms1_ats_log_vol,
                                                                                           cms1_fixing_time,
                                                                                           1.0,
                                                                                           SRT_PUT,
                                                                                           SRT_PREMIUM);

                                //Calculation of the strike for the shifted option (B+dB or B-dB)
                                if(above_below == 1)
                                {
                                        if(pay_rec == 0)
                                                shifted_strike = exo_barriers[i] + call_spread;
                                        else
                                                shifted_strike = exo_barriers[i] - call_spread;
                                }
                                else
                                {
                                        if(pay_rec == 0)
                                                shifted_strike = exo_barriers[i] - call_spread;
                                        else
                                                shifted_strike = exo_barriers[i] + call_spread;
                                }

                                //CMS1 ATS shifted B log vol
                                err = digital_quanto_float_leg_get_lognormal_vol(cms1_vc_name,
                                                                                                                                 fwd1,
                                                                                                                                 cms1_fixing_time,
                                                                                                                                 cms1_start_date,
                                                                                                                                 cms1_theo_end_date,
                                                                                                                                 shifted_strike,
                                                                                                                                 &cms1_ats_log_vol);
                                if (err) goto FREE_RETURN;

                                //Pricing of the put at shifted_strike
                                put_margin_B_shift = srt_f_optblksch(cms1,
                                                                                                         shifted_strike,
                                                                                                         cms1_ats_log_vol,
                                                                                                         cms1_fixing_time,
                                                                                                         1.0,
                                                                                                         SRT_PUT,
                                                                                                         SRT_PREMIUM);

                                //Price the digital quanto float coupon
                                exotic_coupon_PV = exo_gearings[i] * exo_cms2_past_fixings[i] +
           exo_margins[i]; if(above_below == 1) //CMS1<B
                                {
                                        if(pay_rec == 0) //PAY
                                                exotic_coupon_PV *= - (put_margin_B_shift -
           put_margin_B) / call_spread; else //REC exotic_coupon_PV *= (put_margin_B -
           put_margin_B_shift) / call_spread;
                                }
                                else //CMS1>B
                                {
                                        if(pay_rec == 0) //PAY
                                                exotic_coupon_PV *= - 1 + (put_B + put_margin_B -
           put_B_shift - put_margin_B_shift) / call_spread; else //REC exotic_coupon_PV *= 1 -
           (put_B_shift + put_margin_B_shift - put_B - put_margin_B) / call_spread;
                                }
                        }
                }*/

        // Parsing the exotic leg basis
        err = interp_basis(exo_basis[i], &srt_exo_basis_code);
        if (err)
            goto FREE_RETURN;

        cvg    = coverage(exo_start_dates[i], exo_end_dates[i], srt_exo_basis_code);
        df_pay = swp_f_df(today, exo_pay_dates[i], dom_yc_name);

        exotic_coupon1_PV *= df_pay * cvg * exo_notionals[i];
        exotic_coupon2_PV *= df_pay * cvg * exo_notionals[i];
        exotic_coupon3_PV *= df_pay * cvg * exo_notionals[i];
        exotic_coupon4_PV *= df_pay * cvg * exo_notionals[i];

        // Adds the coupon to the sum
        exotic_leg1_PV += exotic_coupon1_PV;
        exotic_leg2_PV += exotic_coupon2_PV;
        exotic_leg3_PV += exotic_coupon3_PV;
        exotic_leg4_PV += exotic_coupon4_PV;

        i++;
    }

    // Pricing of coupons where both CMS have not yet fixed and paid
    while (i < n_exo_coupons)
    {
        // Calculates the start and theo end date of the CMS1 underlying
        cms1_start_date =
            add_unit(exo_cms1_fixing_dates[i], cms1_spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
        if (cms1_tenor_unit == 'M')
        {
            cms1_theo_end_date =
                add_unit(cms1_start_date, cms1_tenor_num, SRT_MONTH, NO_BUSDAY_CONVENTION);
        }
        else
        {
            cms1_theo_end_date =
                add_unit(cms1_start_date, cms1_tenor_num, SRT_YEAR, NO_BUSDAY_CONVENTION);
        }

        // Calculates the start and theo end date of the CMS2 underlying
        cms2_start_date =
            add_unit(exo_cms2_fixing_dates[i], cms2_spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
        if (cms2_tenor_unit == 'M')
        {
            cms2_theo_end_date =
                add_unit(cms2_start_date, cms2_tenor_num, SRT_MONTH, NO_BUSDAY_CONVENTION);
        }
        else
        {
            cms2_theo_end_date =
                add_unit(cms2_start_date, cms2_tenor_num, SRT_YEAR, NO_BUSDAY_CONVENTION);
        }

        // Calculates the fixing times & delays for CMS1 and CMS2
        cms1_fixing_time = (exo_cms1_fixing_dates[i] - today) * YEARS_IN_DAY;
        cms2_fixing_time = (exo_cms2_fixing_dates[i] - today) * YEARS_IN_DAY;
        cms1_start_time  = (cms1_start_date - today) * YEARS_IN_DAY;
        cms2_start_time  = (cms2_start_date - today) * YEARS_IN_DAY;
        cms1_delay       = (exo_pay_dates[i] - cms1_start_date) * YEARS_IN_DAY;
        cms2_delay       = (exo_pay_dates[i] - cms2_start_date) * YEARS_IN_DAY;
        cms1_nfp =
            ((int)((cms1_theo_end_date - cms1_start_date) * YEARS_IN_DAY * cms1_compounding + 0.1)) *
            1.0;
        cms2_nfp =
            ((int)((cms2_theo_end_date - cms2_start_date) * YEARS_IN_DAY * cms2_compounding + 0.1)) *
            1.0;
        fwd_spread1 = swp_f_spread(cms1_start_date, cms1_theo_end_date, cms1_refrate);
        fwd_spread2 = swp_f_spread(cms2_start_date, cms2_theo_end_date, cms2_refrate);

        // Calculates the forward
        err = swp_f_ForwardRate(
            cms1_start_date,
            cms1_theo_end_date,
            cms1_frequency,
            cms1_basis,
            cms1_yc_name,
            cms1_refrate,
            &fwd1);
        if (err)
            goto FREE_RETURN;

        // Calculates the forward
        err = swp_f_ForwardRate(
            cms2_start_date,
            cms2_theo_end_date,
            cms2_frequency,
            cms2_basis,
            cms2_yc_name,
            cms2_refrate,
            &fwd2);
        if (err)
            goto FREE_RETURN;

        err = double_trigger_swap_quanto_cms_rate_quanto(
            fwd1,
            cms1_fixing_time,
            cms1_start_time,
            cms1_nfp,
            cms1_compounding,
            cms1_delay,
            cms1_dRateConv,
            cms1_vol_type,
            1,
            cms1_start_date,
            cms1_theo_end_date,
            cms1_cash_vol,
            fwd_spread1,
            cms1_vc_name,
            n_cms1_strikes_in_vol,
            cms1_strikes_in_vol,
            cms1_dom_for,
            cms1_quanto_corr_times,
            cms1_quanto_corr_values,
            n_cms1_quanto_corr,
            fwd_fx_vols_times,
            fwd_fx_vols_values,
            n_fwd_fx_vols,
            use_CMSRate_CMSTEC,
            &cms1);
        if (err)
            goto FREE_RETURN;

        err = double_trigger_swap_quanto_cms_rate_quanto(
            fwd2,
            cms2_fixing_time,
            cms2_start_time,
            cms2_nfp,
            cms2_compounding,
            cms2_delay,
            cms2_dRateConv,
            cms2_vol_type,
            1,
            cms2_start_date,
            cms2_theo_end_date,
            cms2_cash_vol,
            fwd_spread2,
            cms2_vc_name,
            n_cms2_strikes_in_vol,
            cms2_strikes_in_vol,
            cms2_dom_for,
            cms2_quanto_corr_times,
            cms2_quanto_corr_values,
            n_cms2_quanto_corr,
            fwd_fx_vols_times,
            fwd_fx_vols_values,
            n_fwd_fx_vols,
            use_CMSRate_CMSTEC,
            &cms2);
        if (err)
            goto FREE_RETURN;

        // ---------------------------------------------
        // Prices the first coupon m12 * E[1{X<B}1{Y<K}]
        // ---------------------------------------------
        if (pay_rec == 0)
        {
            B1 = exo_cms1_barriers[i];
            B2 = exo_cms1_barriers[i] + call_spread;
            K1 = exo_cms2_barriers[i];
            K2 = exo_cms2_barriers[i] + call_spread;
        }
        else
        {
            B1 = exo_cms1_barriers[i] - call_spread;
            B2 = exo_cms1_barriers[i];
            K1 = exo_cms2_barriers[i] - call_spread;
            K2 = exo_cms2_barriers[i];
        }

        // Get all the ATS lognormal vols
        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name,
            fwd1,
            fwd_spread1,
            cms1_fixing_time,
            cms1_start_date,
            cms1_theo_end_date,
            B1,
            &cms1_B1_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name,
            fwd1,
            fwd_spread1,
            cms1_fixing_time,
            cms1_start_date,
            cms1_theo_end_date,
            B2,
            &cms1_B2_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name,
            fwd2,
            fwd_spread2,
            cms2_fixing_time,
            cms2_start_date,
            cms2_theo_end_date,
            K1,
            &cms2_K1_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name,
            fwd2,
            fwd_spread2,
            cms2_fixing_time,
            cms2_start_date,
            cms2_theo_end_date,
            K2,
            &cms2_K2_log_vol);
        if (err)
            goto FREE_RETURN;

        // Correlation
        err = constant_or_linear_interpolation(
            cms1_cms2_corr_times,
            cms1_cms2_corr_values_bid,
            n_cms1_cms2_corr,
            LINEAR_METHOD,
            min(cms1_fixing_time, cms2_fixing_time),
            &cms1_cms2_correlation);

        // Price the first coupon
        err = double_inf_barrier_lognormal(
            cms1,
            B1,
            cms1_B1_log_vol,
            B2,
            cms1_B2_log_vol,
            cms1_fixing_time,
            cms2,
            K1,
            cms2_K1_log_vol,
            K2,
            cms2_K2_log_vol,
            cms2_fixing_time,
            cms1_cms2_correlation,
            &exotic_coupon1_PV);
        if (err)
            goto FREE_RETURN;

        exotic_coupon1_PV *= exo_m12[i];

        // --------------------------------------------------
        // Prices the second coupon E[(g2*X+m2)1{X>=B}1{Y<K}]
        // --------------------------------------------------

        if (pay_rec == 0)
        {
            B1 = exo_cms1_barriers[i] - call_spread;
            B2 = exo_cms1_barriers[i];
            K1 = exo_cms2_barriers[i];
            K2 = exo_cms2_barriers[i] + call_spread;
        }
        else
        {
            B1 = exo_cms1_barriers[i];
            B2 = exo_cms1_barriers[i] + call_spread;
            K1 = exo_cms2_barriers[i] - call_spread;
            K2 = exo_cms2_barriers[i];
        }

        // Get all the ATS lognormal vols
        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name,
            fwd1,
            fwd_spread1,
            cms1_fixing_time,
            cms1_start_date,
            cms1_theo_end_date,
            B1,
            &cms1_B1_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name,
            fwd1,
            fwd_spread1,
            cms1_fixing_time,
            cms1_start_date,
            cms1_theo_end_date,
            B2,
            &cms1_B2_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name,
            fwd2,
            fwd_spread2,
            cms2_fixing_time,
            cms2_start_date,
            cms2_theo_end_date,
            K1,
            &cms2_K1_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name,
            fwd2,
            fwd_spread2,
            cms2_fixing_time,
            cms2_start_date,
            cms2_theo_end_date,
            K2,
            &cms2_K2_log_vol);
        if (err)
            goto FREE_RETURN;

        // Correlation
        err = constant_or_linear_interpolation(
            cms1_cms2_corr_times,
            cms1_cms2_corr_values_bid,
            n_cms1_cms2_corr,
            LINEAR_METHOD,
            min(cms1_fixing_time, cms2_fixing_time),
            &cms1_cms2_correlation);

        // For the numeraire adjustments, the vols are taken at exo_cms1_barriers and
        // exo_cms2_barriers
        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name,
            fwd1,
            fwd_spread1,
            cms1_fixing_time,
            cms1_start_date,
            cms1_theo_end_date,
            exo_cms1_barriers[i],
            &cms1_ats_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name,
            fwd2,
            fwd_spread2,
            cms2_fixing_time,
            cms2_start_date,
            cms2_theo_end_date,
            exo_cms2_barriers[i],
            &cms2_ats_log_vol);
        if (err)
            goto FREE_RETURN;

        // Prices the payoff with g2
        cms1_numeraire_adjustment =
            cms1 * exp(cms1_ats_log_vol * cms1_ats_log_vol * cms1_fixing_time);
        cms2_numeraire_adjustment =
            cms2 * exp(cms1_cms2_correlation * cms1_ats_log_vol * cms2_ats_log_vol *
                       min(cms1_fixing_time, cms2_fixing_time));

        put2 = srt_f_optblksch(
            cms2_numeraire_adjustment,
            K2,
            cms2_K2_log_vol,
            cms2_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        put1 = srt_f_optblksch(
            cms2_numeraire_adjustment,
            K1,
            cms2_K1_log_vol,
            cms2_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        exotic_coupon2_PV = exo_g2[i] * cms1 * (put2 - put1) / call_spread;

        err = double_inf_barrier_lognormal(
            cms1_numeraire_adjustment,
            B1,
            cms1_B1_log_vol,
            B2,
            cms1_B2_log_vol,
            cms1_fixing_time,
            cms2_numeraire_adjustment,
            K1,
            cms2_K1_log_vol,
            K2,
            cms2_K2_log_vol,
            cms2_fixing_time,
            cms1_cms2_correlation,
            &double_inf_barrier);
        if (err)
            goto FREE_RETURN;

        exotic_coupon2_PV -= exo_g2[i] * cms1 * double_inf_barrier;

        // Prices the payoff with m2
        put2 =
            srt_f_optblksch(cms2, K2, cms2_K2_log_vol, cms2_fixing_time, 1.0, SRT_PUT, SRT_PREMIUM);

        put1 =
            srt_f_optblksch(cms2, K1, cms2_K1_log_vol, cms2_fixing_time, 1.0, SRT_PUT, SRT_PREMIUM);

        exotic_coupon2_PV += exo_m2[i] * (put2 - put1) / call_spread;

        err = double_inf_barrier_lognormal(
            cms1,
            B1,
            cms1_B1_log_vol,
            B2,
            cms1_B2_log_vol,
            cms1_fixing_time,
            cms2,
            K1,
            cms2_K1_log_vol,
            K2,
            cms2_K2_log_vol,
            cms2_fixing_time,
            cms1_cms2_correlation,
            &double_inf_barrier);
        if (err)
            goto FREE_RETURN;

        exotic_coupon2_PV -= exo_m2[i] * double_inf_barrier;

        // --------------------------------------------------
        // Prices the third coupon E[(g3*Y+m3)1{X<B}1{Y>=K}]
        // --------------------------------------------------

        if (pay_rec == 0)
        {
            B1 = exo_cms1_barriers[i];
            B2 = exo_cms1_barriers[i] + call_spread;
            K1 = exo_cms2_barriers[i] - call_spread;
            K2 = exo_cms2_barriers[i];
        }
        else
        {
            B1 = exo_cms1_barriers[i] - call_spread;
            B2 = exo_cms1_barriers[i];
            K1 = exo_cms2_barriers[i];
            K2 = exo_cms2_barriers[i] + call_spread;
        }

        // Get all the ATS lognormal vols
        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name,
            fwd1,
            fwd_spread1,
            cms1_fixing_time,
            cms1_start_date,
            cms1_theo_end_date,
            B1,
            &cms1_B1_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name,
            fwd1,
            fwd_spread1,
            cms1_fixing_time,
            cms1_start_date,
            cms1_theo_end_date,
            B2,
            &cms1_B2_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name,
            fwd2,
            fwd_spread2,
            cms2_fixing_time,
            cms2_start_date,
            cms2_theo_end_date,
            K1,
            &cms2_K1_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name,
            fwd2,
            fwd_spread2,
            cms2_fixing_time,
            cms2_start_date,
            cms2_theo_end_date,
            K2,
            &cms2_K2_log_vol);
        if (err)
            goto FREE_RETURN;

        // Correlation
        err = constant_or_linear_interpolation(
            cms1_cms2_corr_times,
            cms1_cms2_corr_values_bid,
            n_cms1_cms2_corr,
            LINEAR_METHOD,
            min(cms1_fixing_time, cms2_fixing_time),
            &cms1_cms2_correlation);

        // For the numeraire adjustments, the vols are taken at exo_cms1_barriers and
        // exo_cms2_barriers
        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name,
            fwd1,
            fwd_spread1,
            cms1_fixing_time,
            cms1_start_date,
            cms1_theo_end_date,
            exo_cms1_barriers[i],
            &cms1_ats_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name,
            fwd2,
            fwd_spread2,
            cms2_fixing_time,
            cms2_start_date,
            cms2_theo_end_date,
            exo_cms2_barriers[i],
            &cms2_ats_log_vol);
        if (err)
            goto FREE_RETURN;

        // Prices the payoff with g3
        cms1_numeraire_adjustment =
            cms1 * exp(cms1_cms2_correlation * cms1_ats_log_vol * cms2_ats_log_vol *
                       min(cms1_fixing_time, cms2_fixing_time));
        cms2_numeraire_adjustment =
            cms2 * exp(cms2_ats_log_vol * cms2_ats_log_vol * cms2_fixing_time);

        put2 = srt_f_optblksch(
            cms1_numeraire_adjustment,
            B2,
            cms1_B2_log_vol,
            cms1_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        put1 = srt_f_optblksch(
            cms1_numeraire_adjustment,
            B1,
            cms1_B1_log_vol,
            cms1_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        exotic_coupon3_PV = exo_g3[i] * cms2 * (put2 - put1) / call_spread;

        err = double_inf_barrier_lognormal(
            cms1_numeraire_adjustment,
            B1,
            cms1_B1_log_vol,
            B2,
            cms1_B2_log_vol,
            cms1_fixing_time,
            cms2_numeraire_adjustment,
            K1,
            cms2_K1_log_vol,
            K2,
            cms2_K2_log_vol,
            cms2_fixing_time,
            cms1_cms2_correlation,
            &double_inf_barrier);
        if (err)
            goto FREE_RETURN;

        exotic_coupon3_PV -= exo_g3[i] * cms2 * double_inf_barrier;

        // Prices the payoff with m3
        put2 =
            srt_f_optblksch(cms1, B2, cms1_B2_log_vol, cms1_fixing_time, 1.0, SRT_PUT, SRT_PREMIUM);

        put1 =
            srt_f_optblksch(cms1, B1, cms1_B1_log_vol, cms1_fixing_time, 1.0, SRT_PUT, SRT_PREMIUM);

        exotic_coupon3_PV += exo_m3[i] * (put2 - put1) / call_spread;

        err = double_inf_barrier_lognormal(
            cms1,
            B1,
            cms1_B1_log_vol,
            B2,
            cms1_B2_log_vol,
            cms1_fixing_time,
            cms2,
            K1,
            cms2_K1_log_vol,
            K2,
            cms2_K2_log_vol,
            cms2_fixing_time,
            cms1_cms2_correlation,
            &double_inf_barrier);
        if (err)
            goto FREE_RETURN;

        exotic_coupon3_PV -= exo_m3[i] * double_inf_barrier;

        // ----------------------------------------------------------------------
        // Prices the fourth coupon 0.5 * E[(g4_X*X + g4_Y*Y + m4)1{X>=B}1{Y>=K}]
        // ----------------------------------------------------------------------

        if (pay_rec == 0)
        {
            B1 = exo_cms1_barriers[i] - call_spread;
            B2 = exo_cms1_barriers[i];
            K1 = exo_cms2_barriers[i] - call_spread;
            K2 = exo_cms2_barriers[i];
        }
        else
        {
            B1 = exo_cms1_barriers[i];
            B2 = exo_cms1_barriers[i] + call_spread;
            K1 = exo_cms2_barriers[i];
            K2 = exo_cms2_barriers[i] + call_spread;
        }

        // Get all the ATS lognormal vols
        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name,
            fwd1,
            fwd_spread1,
            cms1_fixing_time,
            cms1_start_date,
            cms1_theo_end_date,
            B1,
            &cms1_B1_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name,
            fwd1,
            fwd_spread1,
            cms1_fixing_time,
            cms1_start_date,
            cms1_theo_end_date,
            B2,
            &cms1_B2_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name,
            fwd2,
            fwd_spread2,
            cms2_fixing_time,
            cms2_start_date,
            cms2_theo_end_date,
            K1,
            &cms2_K1_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name,
            fwd2,
            fwd_spread2,
            cms2_fixing_time,
            cms2_start_date,
            cms2_theo_end_date,
            K2,
            &cms2_K2_log_vol);
        if (err)
            goto FREE_RETURN;

        // Correlation
        err = constant_or_linear_interpolation(
            cms1_cms2_corr_times,
            cms1_cms2_corr_values_bid,
            n_cms1_cms2_corr,
            LINEAR_METHOD,
            min(cms1_fixing_time, cms2_fixing_time),
            &cms1_cms2_correlation);

        // For the numeraire adjustments, the vols are taken at exo_cms1_barriers and
        // exo_cms2_barriers
        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name,
            fwd1,
            fwd_spread1,
            cms1_fixing_time,
            cms1_start_date,
            cms1_theo_end_date,
            exo_cms1_barriers[i],
            &cms1_ats_log_vol);
        if (err)
            goto FREE_RETURN;

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name,
            fwd2,
            fwd_spread2,
            cms2_fixing_time,
            cms2_start_date,
            cms2_theo_end_date,
            exo_cms2_barriers[i],
            &cms2_ats_log_vol);
        if (err)
            goto FREE_RETURN;

        // Prices the payoff with g4_X
        cms1_numeraire_adjustment =
            cms1 * exp(cms1_ats_log_vol * cms1_ats_log_vol * cms1_fixing_time);
        cms2_numeraire_adjustment =
            cms2 * exp(cms1_cms2_correlation * cms1_ats_log_vol * cms2_ats_log_vol *
                       min(cms1_fixing_time, cms2_fixing_time));

        exotic_coupon4_PV = exo_g4_X[i] * cms1;

        put2 = srt_f_optblksch(
            cms2_numeraire_adjustment,
            K2,
            cms2_K2_log_vol,
            cms2_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        put1 = srt_f_optblksch(
            cms2_numeraire_adjustment,
            K1,
            cms2_K1_log_vol,
            cms2_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        exotic_coupon4_PV -= exo_g4_X[i] * cms1 * (put2 - put1) / call_spread;

        put2 = srt_f_optblksch(
            cms1_numeraire_adjustment,
            B2,
            cms1_B2_log_vol,
            cms1_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        put1 = srt_f_optblksch(
            cms1_numeraire_adjustment,
            B1,
            cms1_B1_log_vol,
            cms1_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        exotic_coupon4_PV -= exo_g4_X[i] * cms1 * (put2 - put1) / call_spread;

        err = double_inf_barrier_lognormal(
            cms1_numeraire_adjustment,
            B1,
            cms1_B1_log_vol,
            B2,
            cms1_B2_log_vol,
            cms1_fixing_time,
            cms2_numeraire_adjustment,
            K1,
            cms2_K1_log_vol,
            K2,
            cms2_K2_log_vol,
            cms2_fixing_time,
            cms1_cms2_correlation,
            &double_inf_barrier);
        if (err)
            goto FREE_RETURN;

        exotic_coupon4_PV += exo_g4_X[i] * cms1 * double_inf_barrier;

        // Prices the payoff with g4_Y
        cms1_numeraire_adjustment =
            cms1 * exp(cms1_cms2_correlation * cms1_ats_log_vol * cms2_ats_log_vol *
                       min(cms1_fixing_time, cms2_fixing_time));
        cms2_numeraire_adjustment =
            cms2 * exp(cms2_ats_log_vol * cms2_ats_log_vol * cms2_fixing_time);

        exotic_coupon4_PV += exo_g4_Y[i] * cms2;

        put2 = srt_f_optblksch(
            cms2_numeraire_adjustment,
            K2,
            cms2_K2_log_vol,
            cms2_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        put1 = srt_f_optblksch(
            cms2_numeraire_adjustment,
            K1,
            cms2_K1_log_vol,
            cms2_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        exotic_coupon4_PV -= exo_g4_Y[i] * cms2 * (put2 - put1) / call_spread;

        put2 = srt_f_optblksch(
            cms1_numeraire_adjustment,
            B2,
            cms1_B2_log_vol,
            cms1_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        put1 = srt_f_optblksch(
            cms1_numeraire_adjustment,
            B1,
            cms1_B1_log_vol,
            cms1_fixing_time,
            1.0,
            SRT_PUT,
            SRT_PREMIUM);

        exotic_coupon4_PV -= exo_g4_Y[i] * cms2 * (put2 - put1) / call_spread;

        err = double_inf_barrier_lognormal(
            cms1_numeraire_adjustment,
            B1,
            cms1_B1_log_vol,
            B2,
            cms1_B2_log_vol,
            cms1_fixing_time,
            cms2_numeraire_adjustment,
            K1,
            cms2_K1_log_vol,
            K2,
            cms2_K2_log_vol,
            cms2_fixing_time,
            cms1_cms2_correlation,
            &double_inf_barrier);
        if (err)
            goto FREE_RETURN;

        exotic_coupon4_PV += exo_g4_Y[i] * cms2 * double_inf_barrier;

        // Prices the payoff with m4

        exotic_coupon4_PV += exo_m4[i];

        put2 =
            srt_f_optblksch(cms2, K2, cms2_K2_log_vol, cms2_fixing_time, 1.0, SRT_PUT, SRT_PREMIUM);

        put1 =
            srt_f_optblksch(cms2, K1, cms2_K1_log_vol, cms2_fixing_time, 1.0, SRT_PUT, SRT_PREMIUM);

        exotic_coupon4_PV -= exo_m4[i] * (put2 - put1) / call_spread;

        put2 =
            srt_f_optblksch(cms1, B2, cms1_B2_log_vol, cms1_fixing_time, 1.0, SRT_PUT, SRT_PREMIUM);

        put1 =
            srt_f_optblksch(cms1, B1, cms1_B1_log_vol, cms1_fixing_time, 1.0, SRT_PUT, SRT_PREMIUM);

        exotic_coupon4_PV -= exo_m4[i] * (put2 - put1) / call_spread;

        err = double_inf_barrier_lognormal(
            cms1,
            B1,
            cms1_B1_log_vol,
            B2,
            cms1_B2_log_vol,
            cms1_fixing_time,
            cms2,
            K1,
            cms2_K1_log_vol,
            K2,
            cms2_K2_log_vol,
            cms2_fixing_time,
            cms1_cms2_correlation,
            &double_inf_barrier);
        if (err)
            goto FREE_RETURN;

        exotic_coupon4_PV += exo_m4[i] * double_inf_barrier;
        exotic_coupon4_PV *= 0.5;

        // ----------------------------------------------------------------
        // End of the pricing of the four coupons
        // ----------------------------------------------------------------

        // Parsing the exotic leg basis
        err = interp_basis(exo_basis[i], &srt_exo_basis_code);
        if (err)
            goto FREE_RETURN;

        cvg    = coverage(exo_start_dates[i], exo_end_dates[i], srt_exo_basis_code);
        df_pay = swp_f_df(today, exo_pay_dates[i], dom_yc_name);
        exotic_coupon1_PV *= df_pay * cvg * exo_notionals[i];
        exotic_coupon2_PV *= df_pay * cvg * exo_notionals[i];
        exotic_coupon3_PV *= df_pay * cvg * exo_notionals[i];
        exotic_coupon4_PV *= df_pay * cvg * exo_notionals[i];

        // Adds the coupon to the sum
        exotic_leg1_PV += exotic_coupon1_PV;
        exotic_leg2_PV += exotic_coupon2_PV;
        exotic_leg3_PV += exotic_coupon3_PV;
        exotic_leg4_PV += exotic_coupon4_PV;

        // Increment the coupon and continue
        i++;
    }

// Returns the final result
FREE_RETURN:
    *cpn1                       = exotic_leg1_PV;
    *cpn2                       = exotic_leg2_PV;
    *cpn3                       = exotic_leg3_PV;
    *cpn4                       = exotic_leg4_PV;
    *DoubleTriggerSwapQuantoLeg = (*cpn1) + (*cpn2) + (*cpn3) + (*cpn4);

    if (pay_rec == 0)
    {
        (*cpn1) *= -1.0;
        (*cpn2) *= -1.0;
        (*cpn3) *= -1.0;
        (*cpn4) *= -1.0;
        (*DoubleTriggerSwapQuantoLeg) *= -1.0;
    }

    return err;
}
