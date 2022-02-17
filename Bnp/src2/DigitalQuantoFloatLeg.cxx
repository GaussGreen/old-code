/* ===================================================================================
   FILENAME:      DigitalQuantoFloatLeg.cxx / LY J-M

   PURPOSE:       Computes the PV of a DigitalQuantoFloatLeg
   ===================================================================================
 */

#pragma warning(disable : 4786) // Disable long name warnings

#include "DigitalQuantoFloatLeg.h"
#include "math.h"
#include "num_h_allhdr.h"
#include "opHeston.h"
#include "opfnctns.h"
#include "swp_h_all.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "swp_h_vol.h"

// From frequency "A" and "BB" get the numerical compounding 1 and the dRate
Err frequency_and_basis_to_compounding_and_dRateConv(
    char *frequency, char *basis, SrtCompounding *compounding,
    double *dRateConv) {
  // Cf CTSProdStruct ref_conv
  Err err;
  SrtBasisCode srt_basis_code;

  err = interp_compounding(frequency, compounding);
  if (err)
    goto FREE_RETURN;

  err = interp_basis(basis, &srt_basis_code);
  if (err)
    goto FREE_RETURN;

  if (srt_basis_code == BASIS_ACT_360)
    *dRateConv = 365.0 / 360.0;
  else
    *dRateConv = 1.0;

FREE_RETURN:
  return err;
}

// Interpolate linearly or piecewise constant
Err constant_or_linear_interpolation(double *dates, double *values, int n_dates,
                                     int method, double date_wanted,
                                     double *value_wanted) {
  int i = 0;
  double result;

  // Flat intrapolation if before
  if (date_wanted <= dates[0])
    result = values[0];
  else {
    // Checks if the index is below the limit before comparing the dates
    // to avoid release compilation problems
    while (i < n_dates && date_wanted > dates[i])
      i++;
    i--;

    if (i < n_dates - 1) {
      // Interpolate linearly
      result = values[i] + (values[i + 1] - values[i]) /
                               (dates[i + 1] - dates[i]) *
                               (date_wanted - dates[i]);
    } else {
      // Flat extrapolation
      result = values[n_dates - 1];
    }
  }
  *value_wanted = result;
  return NULL;
}

// Calculates the CMS end date using the tenor and the start date
Err cms_get_tenor_num_and_unit_from_tenor(char *cms_tenor, int *tenor_num,
                                          char *tenor_unit) {
  int string_len;
  char tenor[256]; // "5Y"

  strcpy(tenor, cms_tenor);
  string_len = strlen(tenor);
  *tenor_unit = tenor[string_len - 1];
  tenor[string_len - 1] = '\0';
  *tenor_num = atoi(tenor);

  return NULL;
}

// Gets the lognormal vol (to be modified for MAD compatibility by using the MAD
// function used in the swp_f_truncvol function used in swp_f_cms_option
Err digital_quanto_float_leg_get_lognormal_vol(
    char *vc_name, double fwd, double fwd_spread, double fixing_time,
    long start_date, long end_date, double strike, double *result_volatility) {
  Err err;
  double vol, power, temp_vol;

  err = swp_f_vol(vc_name, start_date, end_date, strike, &vol, &power);
  if (err)
    goto FREE_RETURN;

  // Conversion from normal to lognormal
  if (power == 0.0) {
    err = srt_f_optsarbvol(fwd, strike, fixing_time, vol, 0, 0, 0, SRT_BETAVOL,
                           SRT_LOGNORMAL, &temp_vol);

    if (err)
      goto FREE_RETURN;

    vol = temp_vol;
  }

FREE_RETURN:
  *result_volatility = vol;
  return err;
}

// Calculates a CMS rate quanto or non quanto
Err digital_quanto_float_leg_cms_rate_quanto(
    double fwd, double cms_fixing_time, double cms_start_time, double cms_nfp,

    SrtCompounding cms_compounding, double cms_delay, double cms_dRateConv,

    SrtDiffusionType cms_vol_type, int cms_flat_vol_or_smile_approx,
    long cms_start_date, long cms_theo_end_date, int cms_cash_vol,
    double fwd_spread, char *cms_vc_name, int n_cms_strikes_in_vol,
    double *cms_strikes_in_vol,

    int cms_dom_for,

    double *cms_quanto_corr_times, double *cms_quanto_corr_values,
    int n_cms_quanto_corr,

    double *fwd_fx_vols_times, double *fwd_fx_vols_values, int n_fwd_fx_vols,

    double *result_cms_rate_quanto) {
  Err err;
  double cms, cms_atm_log_vol, cms_quanto_correlation, fwd_fx_vol, flat_vol,
      temp_vol;

  if (cms_flat_vol_or_smile_approx == 0) {
    // Calculates the CMS Rate with Flat Vol
    err = digital_quanto_float_leg_get_lognormal_vol(
        cms_vc_name, fwd, fwd_spread, cms_fixing_time, cms_start_date,
        cms_theo_end_date, fwd, &flat_vol);
    if (err)
      goto FREE_RETURN;

    if (cms_vol_type == SRT_NORMAL) {
      err = srt_f_optsarbvol(fwd, fwd, cms_fixing_time, flat_vol, 0, 1, 0,
                             SRT_BETAVOL, SRT_NORMAL, &temp_vol);

      if (err)
        goto FREE_RETURN;

      flat_vol = temp_vol;
    }

    err = swp_f_Cms_Rate(fwd, cms_fixing_time, cms_nfp, cms_compounding,
                         cms_delay, cms_dRateConv, cms_vol_type, flat_vol,
                         0, // 0: Use flat vol
                         cms_start_date, cms_theo_end_date, cms_cash_vol,
                         fwd_spread, cms_vc_name, n_cms_strikes_in_vol,
                         cms_strikes_in_vol, &cms);

    if (err)
      goto FREE_RETURN;
  } else {
    // Calculates the CMS Rate with Smile Approx
    err = swp_f_Cms_Rate(fwd, cms_fixing_time, cms_nfp, cms_compounding,
                         cms_delay, cms_dRateConv, cms_vol_type, 0.0,
                         1, // 1: Smile linear interpolation
                         cms_start_date, cms_theo_end_date, cms_cash_vol,
                         fwd_spread, cms_vc_name, n_cms_strikes_in_vol,
                         cms_strikes_in_vol, &cms);

    if (err)
      goto FREE_RETURN;
  }

  if (cms_dom_for == 1) {
    // CMS ATM log vol
    err = digital_quanto_float_leg_get_lognormal_vol(
        cms_vc_name, fwd, fwd_spread, cms_fixing_time, cms_start_date,
        cms_theo_end_date, fwd, &cms_atm_log_vol);

    if (err)
      goto FREE_RETURN;

    err = constant_or_linear_interpolation(
        cms_quanto_corr_times, cms_quanto_corr_values, n_cms_quanto_corr,
        LINEAR_METHOD, cms_start_time, &cms_quanto_correlation);

    if (err)
      goto FREE_RETURN;

    err = constant_or_linear_interpolation(
        fwd_fx_vols_times, fwd_fx_vols_values, n_fwd_fx_vols, LINEAR_METHOD,
        cms_start_time, &fwd_fx_vol);

    if (err)
      goto FREE_RETURN;

    cms *= exp(-cms_quanto_correlation * cms_atm_log_vol * fwd_fx_vol *
               cms_fixing_time);
  }

FREE_RETURN:
  *result_cms_rate_quanto = cms;
  return err;
}

// Calculates the PV of a digital quanto float leg
Err digital_quanto_float_leg(
    char *dom_yc_name, char *dom_vc_name, int dom_spot_lag,

    double *dom_strikes_in_vol, int n_dom_strikes_in_vol,
    SrtDiffusionType dom_vol_type, int dom_cash_vol,

    char *for_yc_name, char *for_vc_name, int for_spot_lag,

    double *for_strikes_in_vol, int n_for_strikes_in_vol,
    SrtDiffusionType for_vol_type, int for_cash_vol,

    long today,

    char *cms1_tenor, char *cms1_frequency, char *cms1_basis,
    char *cms1_refrate, int cms1_dom_for,

    char *cms2_tenor, char *cms2_frequency, char *cms2_basis,
    char *cms2_refrate, int cms2_dom_for,

    long *exo_start_dates, long *exo_end_dates, long *exo_pay_dates,
    char **exo_basis, double *exo_notionals, long *exo_cms1_fixing_dates,
    long *exo_cms2_fixing_dates, double *exo_gearings, double *exo_barriers,
    double *exo_margins, double *exo_cms1_past_fixings,
    double *exo_cms2_past_fixings, int n_exo_coupons,

    double *cms1_cms2_corr_times, double *cms1_cms2_corr_values_bid,
    double *cms1_cms2_corr_values_ask, int n_cms1_cms2_corr,

    double *cms1_quanto_corr_times, double *cms1_quanto_corr_values,
    int n_cms1_quanto_corr,

    double *cms2_quanto_corr_times, double *cms2_quanto_corr_values,
    int n_cms2_quanto_corr,

    double *fwd_fx_vols_times, double *fwd_fx_vols_values, int n_fwd_fx_vols,

    int pay_rec, int above_below, double call_spread, int accrue_on_barrier,
    int atm_ats_std, int cms_flat_vol_or_smile_approx,

    int eod_fix_flag, int eod_pay_flag,

    double *DigitalQuantoFloatLeg) {
  Err err;
  double fwd1, fwd2, cms1, cms2, fwd_spread1, fwd_spread2;
  double cms1_fixing_time, cms1_delay, cms1_start_time, cms1_nfp,
      cms1_dRateConv;
  double cms2_fixing_time, cms2_delay, cms2_start_time, cms2_nfp,
      cms2_dRateConv;
  long cms1_start_date, cms2_start_date, cms1_theo_end_date, cms2_theo_end_date;
  char cms1_tenor_unit, cms2_tenor_unit;
  SrtCompounding cms1_compounding, cms2_compounding;
  SrtBasisCode srt_exo_basis_code;
  double cms1_ats_log_vol, cms1_atm_log_vol, cms2_atm_log_vol, cms2_log_vol,
      cms1_nb_std, cms2_strike;
  double numeraire_adjustment = 1.0, cms1_cms2_correlation;
  double shifted_strike, df_pay, cvg;
  int i, cms1_tenor_num, cms2_tenor_num;

  // CMS1 market details
  char cms1_yc_name[256];
  char cms1_vc_name[256];
  int cms1_spot_lag;

  double *cms1_strikes_in_vol = NULL;
  int n_cms1_strikes_in_vol;
  SrtDiffusionType cms1_vol_type;
  int cms1_cash_vol;

  // CMS2 market details
  char cms2_yc_name[256];
  char cms2_vc_name[256];
  int cms2_spot_lag;

  double *cms2_strikes_in_vol = NULL;
  int n_cms2_strikes_in_vol;
  SrtDiffusionType cms2_vol_type;
  int cms2_cash_vol;

  // Results
  double put_B = 0.0;
  double put_margin_B = 0.0;
  double put_B_shift = 0.0;
  double put_margin_B_shift = 0.0;
  double exotic_leg_PV = 0.0;
  double exotic_coupon_PV = 0.0;

  // Checks that fixing & pay dates are increasing
  i = 0;
  while (i < n_exo_coupons - 1) {
    if ((exo_cms1_fixing_dates[i] >= exo_cms1_fixing_dates[i + 1]) ||
        (exo_cms2_fixing_dates[i] >= exo_cms2_fixing_dates[i + 1])) {
      err = "cms_fixing_dates are not in increasing order";
      goto FREE_RETURN;
    }
    i++;
  }

  // Checks that fixing dates <= pay dates
  i = 0;
  while (i < n_exo_coupons) {
    if ((exo_cms1_fixing_dates[i] > exo_pay_dates[i]) ||
        (exo_cms2_fixing_dates[i] > exo_pay_dates[i])) {
      err = "cms is paid before the before the fixing date";
      goto FREE_RETURN;
    }
    i++;
  }

  // Reading CMS1 market parameters
  if (cms1_dom_for == 0) {
    // CMS1 is domestic
    strcpy(cms1_yc_name, dom_yc_name);
    strcpy(cms1_vc_name, dom_vc_name);
    cms1_spot_lag = dom_spot_lag;

    cms1_strikes_in_vol = dom_strikes_in_vol;
    n_cms1_strikes_in_vol = n_dom_strikes_in_vol;
    cms1_vol_type = dom_vol_type;
    cms1_cash_vol = dom_cash_vol;
  } else {
    // CMS1 is foreign
    strcpy(cms1_yc_name, for_yc_name);
    strcpy(cms1_vc_name, for_vc_name);
    cms1_spot_lag = for_spot_lag;

    cms1_strikes_in_vol = for_strikes_in_vol;
    n_cms1_strikes_in_vol = n_for_strikes_in_vol;
    cms1_vol_type = for_vol_type;
    cms1_cash_vol = for_cash_vol;
  }

  // Reading CMS2 market parameters
  if (cms2_dom_for == 0) {
    // CMS2 is domestic
    strcpy(cms2_yc_name, dom_yc_name);
    strcpy(cms2_vc_name, dom_vc_name);
    cms2_spot_lag = dom_spot_lag;

    cms2_strikes_in_vol = dom_strikes_in_vol;
    n_cms2_strikes_in_vol = n_dom_strikes_in_vol;
    cms2_vol_type = dom_vol_type;
    cms2_cash_vol = dom_cash_vol;
  } else {
    // CMS2 is foreign
    strcpy(cms2_yc_name, for_yc_name);
    strcpy(cms2_vc_name, for_vc_name);
    cms2_spot_lag = for_spot_lag;

    cms2_strikes_in_vol = for_strikes_in_vol;
    n_cms2_strikes_in_vol = n_for_strikes_in_vol;
    cms2_vol_type = for_vol_type;
    cms2_cash_vol = for_cash_vol;
  }

  // Computes the compounding and the day rate convention
  frequency_and_basis_to_compounding_and_dRateConv(
      cms1_frequency, cms1_basis, &cms1_compounding, &cms1_dRateConv);

  frequency_and_basis_to_compounding_and_dRateConv(
      cms2_frequency, cms2_basis, &cms2_compounding, &cms2_dRateConv);

  // Parsing the tenor num and tenor unit of the CMS assets
  cms_get_tenor_num_and_unit_from_tenor(cms1_tenor, &cms1_tenor_num,
                                        &cms1_tenor_unit);

  cms_get_tenor_num_and_unit_from_tenor(cms2_tenor, &cms2_tenor_num,
                                        &cms2_tenor_unit);

  // Initializes the coupon index
  i = 0;

  // Remove all periods which have already been paid
  while (exo_pay_dates[i] < today + eod_pay_flag)
    i++;

  // Pricing of coupons where at least one CMS has been fixed
  while (min(exo_cms1_fixing_dates[i], exo_cms2_fixing_dates[i]) <
         today + eod_fix_flag) {
    // Both fixed
    if ((exo_cms1_fixing_dates[i] < today + eod_fix_flag) &&
        (exo_cms2_fixing_dates[i] < today + eod_fix_flag)) {
      if (((exo_cms1_past_fixings[i] > exo_barriers[i]) &&
           (above_below == 0)) ||
          ((exo_cms1_past_fixings[i] < exo_barriers[i]) &&
           (above_below == 1)) ||
          ((exo_cms1_past_fixings[i] == exo_barriers[i]) &&
           (accrue_on_barrier == 1)))
        exotic_coupon_PV =
            exo_gearings[i] * exo_cms2_past_fixings[i] + exo_margins[i];
      else
        exotic_coupon_PV = 0.0;
    }

    // One asset fixed
    else {
      // CMS1 has fixed
      if (exo_cms1_fixing_dates[i] < today + eod_fix_flag) {
        // Payoff is a forward CMS2 + margin (quanto corrected if needed)
        if (((exo_cms1_past_fixings[i] > exo_barriers[i]) &&
             (above_below == 0)) ||
            ((exo_cms1_past_fixings[i] < exo_barriers[i]) &&
             (above_below == 1)) ||
            ((exo_cms1_past_fixings[i] == exo_barriers[i]) &&
             (accrue_on_barrier == 1))) {
          // Calculates the start and theo end date of the CMS2 underlying
          cms2_start_date = add_unit(exo_cms2_fixing_dates[i], cms2_spot_lag,
                                     SRT_BDAY, MODIFIED_SUCCEEDING);
          if (cms2_tenor_unit == 'M') {
            cms2_theo_end_date = add_unit(cms2_start_date, cms2_tenor_num,
                                          SRT_MONTH, NO_BUSDAY_CONVENTION);
          } else {
            cms2_theo_end_date = add_unit(cms2_start_date, cms2_tenor_num,
                                          SRT_YEAR, NO_BUSDAY_CONVENTION);
          }

          // Calculates number of full periods for the underlying of CMS2
          cms2_fixing_time = (exo_cms2_fixing_dates[i] - today) * YEARS_IN_DAY;
          cms2_start_time = (cms2_start_date - today) * YEARS_IN_DAY;
          cms2_delay = (exo_pay_dates[i] - cms2_start_date) * YEARS_IN_DAY;
          cms2_nfp = ((int)((cms2_theo_end_date - cms2_start_date) *
                                YEARS_IN_DAY * cms2_compounding +
                            0.1)) *
                     1.0;
          fwd_spread2 =
              swp_f_spread(cms2_start_date, cms2_theo_end_date, cms2_refrate);

          // Calculates the forward
          err = swp_f_ForwardRate(cms2_start_date, cms2_theo_end_date,
                                  cms2_frequency, cms2_basis, cms2_yc_name,
                                  cms2_refrate, &fwd2);
          if (err)
            goto FREE_RETURN;

          // Calculates the CMS Rate
          err = digital_quanto_float_leg_cms_rate_quanto(
              fwd2, cms2_fixing_time, cms2_start_time, cms2_nfp,
              cms2_compounding, cms2_delay, cms2_dRateConv, cms2_vol_type,
              cms_flat_vol_or_smile_approx, cms2_start_date, cms2_theo_end_date,
              cms2_cash_vol, fwd_spread2, cms2_vc_name, n_cms2_strikes_in_vol,
              cms2_strikes_in_vol, cms2_dom_for, cms2_quanto_corr_times,
              cms2_quanto_corr_values, n_cms2_quanto_corr, fwd_fx_vols_times,
              fwd_fx_vols_values, n_fwd_fx_vols, &cms2);
          if (err)
            goto FREE_RETURN;
          exotic_coupon_PV = exo_gearings[i] * cms2 + exo_margins[i];
        } else
          exotic_coupon_PV = 0.0;
      }

      // CMS2 has fixed we need to price a digital on CMS1
      else {
        // Calculates the start and theo end date of the CMS1 underlying
        cms1_start_date = add_unit(exo_cms1_fixing_dates[i], cms1_spot_lag,
                                   SRT_BDAY, MODIFIED_SUCCEEDING);
        if (cms1_tenor_unit == 'M') {
          cms1_theo_end_date = add_unit(cms1_start_date, cms1_tenor_num,
                                        SRT_MONTH, NO_BUSDAY_CONVENTION);
        } else {
          cms1_theo_end_date = add_unit(cms1_start_date, cms1_tenor_num,
                                        SRT_YEAR, NO_BUSDAY_CONVENTION);
        }

        // Calculates number of full periods for the underlying of CMS1
        cms1_fixing_time = (exo_cms1_fixing_dates[i] - today) * YEARS_IN_DAY;
        cms1_start_time = (cms1_start_date - today) * YEARS_IN_DAY;
        cms1_delay = (exo_pay_dates[i] - cms1_start_date) * YEARS_IN_DAY;
        cms1_nfp = ((int)((cms1_theo_end_date - cms1_start_date) *
                              YEARS_IN_DAY * cms1_compounding +
                          0.1)) *
                   1.0;
        fwd_spread1 =
            swp_f_spread(cms1_start_date, cms1_theo_end_date, cms1_refrate);

        // Calculates the forward
        err = swp_f_ForwardRate(cms1_start_date, cms1_theo_end_date,
                                cms1_frequency, cms1_basis, cms1_yc_name,
                                cms1_refrate, &fwd1);
        if (err)
          goto FREE_RETURN;

        // Calculates the CMS Rate
        err = digital_quanto_float_leg_cms_rate_quanto(
            fwd1, cms1_fixing_time, cms1_start_time, cms1_nfp, cms1_compounding,
            cms1_delay, cms1_dRateConv, cms1_vol_type,
            cms_flat_vol_or_smile_approx, cms1_start_date, cms1_theo_end_date,
            cms1_cash_vol, fwd_spread1, cms1_vc_name, n_cms1_strikes_in_vol,
            cms1_strikes_in_vol, cms1_dom_for, cms1_quanto_corr_times,
            cms1_quanto_corr_values, n_cms1_quanto_corr, fwd_fx_vols_times,
            fwd_fx_vols_values, n_fwd_fx_vols, &cms1);
        if (err)
          goto FREE_RETURN;

        // CMS1 ATS B log vol
        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name, fwd1, fwd_spread1, cms1_fixing_time, cms1_start_date,
            cms1_theo_end_date, exo_barriers[i], &cms1_ats_log_vol);
        if (err)
          goto FREE_RETURN;

        put_margin_B =
            srt_f_optblksch(cms1, exo_barriers[i], cms1_ats_log_vol,
                            cms1_fixing_time, 1.0, SRT_PUT, SRT_PREMIUM);

        // Calculation of the strike for the shifted option (B+dB or B-dB)
        if (above_below == 1) {
          if (pay_rec == 0)
            shifted_strike = exo_barriers[i] + call_spread;
          else
            shifted_strike = exo_barriers[i] - call_spread;
        } else {
          if (pay_rec == 0)
            shifted_strike = exo_barriers[i] - call_spread;
          else
            shifted_strike = exo_barriers[i] + call_spread;
        }

        // CMS1 ATS shifted B log vol
        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name, fwd1, fwd_spread1, cms1_fixing_time, cms1_start_date,
            cms1_theo_end_date, shifted_strike, &cms1_ats_log_vol);
        if (err)
          goto FREE_RETURN;

        // Pricing of the put at shifted_strike
        put_margin_B_shift =
            srt_f_optblksch(cms1, shifted_strike, cms1_ats_log_vol,
                            cms1_fixing_time, 1.0, SRT_PUT, SRT_PREMIUM);

        // Price the digital quanto float coupon
        exotic_coupon_PV =
            exo_gearings[i] * exo_cms2_past_fixings[i] + exo_margins[i];
        if (above_below == 1) // CMS1<B
        {
          if (pay_rec == 0) // PAY
            exotic_coupon_PV *=
                -(put_margin_B_shift - put_margin_B) / call_spread;
          else // REC
            exotic_coupon_PV *=
                (put_margin_B - put_margin_B_shift) / call_spread;
        } else // CMS1>B
        {
          if (pay_rec == 0) // PAY
            exotic_coupon_PV *=
                -1 + (put_B + put_margin_B - put_B_shift - put_margin_B_shift) /
                         call_spread;
          else // REC
            exotic_coupon_PV *=
                1 - (put_B_shift + put_margin_B_shift - put_B - put_margin_B) /
                        call_spread;
        }
      }
    }

    // Parsing the exotic leg basis
    err = interp_basis(exo_basis[i], &srt_exo_basis_code);
    if (err)
      goto FREE_RETURN;

    cvg = coverage(exo_start_dates[i], exo_end_dates[i], srt_exo_basis_code);
    df_pay = swp_f_df(today, exo_pay_dates[i], dom_yc_name);
    exotic_coupon_PV *= df_pay * cvg * exo_notionals[i];

    // Adds the coupon to the sum
    exotic_leg_PV += exotic_coupon_PV;
    i++;
  }

  // Pricing of coupons where both CMS have not yet fixed and paid
  while (i < n_exo_coupons) {

    // Calculates the start and theo end date of the CMS1 underlying
    cms1_start_date = add_unit(exo_cms1_fixing_dates[i], cms1_spot_lag,
                               SRT_BDAY, MODIFIED_SUCCEEDING);
    if (cms1_tenor_unit == 'M') {
      cms1_theo_end_date = add_unit(cms1_start_date, cms1_tenor_num, SRT_MONTH,
                                    NO_BUSDAY_CONVENTION);
    } else {
      cms1_theo_end_date = add_unit(cms1_start_date, cms1_tenor_num, SRT_YEAR,
                                    NO_BUSDAY_CONVENTION);
    }

    // Calculates the start and theo end date of the CMS2 underlying
    cms2_start_date = add_unit(exo_cms2_fixing_dates[i], cms2_spot_lag,
                               SRT_BDAY, MODIFIED_SUCCEEDING);
    if (cms2_tenor_unit == 'M') {
      cms2_theo_end_date = add_unit(cms2_start_date, cms2_tenor_num, SRT_MONTH,
                                    NO_BUSDAY_CONVENTION);
    } else {
      cms2_theo_end_date = add_unit(cms2_start_date, cms2_tenor_num, SRT_YEAR,
                                    NO_BUSDAY_CONVENTION);
    }

    // Calculates the fixing times & delays for CMS1 and CMS2
    cms1_fixing_time = (exo_cms1_fixing_dates[i] - today) * YEARS_IN_DAY;
    cms2_fixing_time = (exo_cms2_fixing_dates[i] - today) * YEARS_IN_DAY;
    cms1_start_time = (cms1_start_date - today) * YEARS_IN_DAY;
    cms2_start_time = (cms2_start_date - today) * YEARS_IN_DAY;
    cms1_delay = (exo_pay_dates[i] - cms1_start_date) * YEARS_IN_DAY;
    cms2_delay = (exo_pay_dates[i] - cms2_start_date) * YEARS_IN_DAY;
    cms1_nfp = ((int)((cms1_theo_end_date - cms1_start_date) * YEARS_IN_DAY *
                          cms1_compounding +
                      0.1)) *
               1.0;
    cms2_nfp = ((int)((cms2_theo_end_date - cms2_start_date) * YEARS_IN_DAY *
                          cms2_compounding +
                      0.1)) *
               1.0;
    fwd_spread1 =
        swp_f_spread(cms1_start_date, cms1_theo_end_date, cms1_refrate);
    fwd_spread2 =
        swp_f_spread(cms2_start_date, cms2_theo_end_date, cms2_refrate);

    // Calculates the forward
    err = swp_f_ForwardRate(cms1_start_date, cms1_theo_end_date, cms1_frequency,
                            cms1_basis, cms1_yc_name, cms1_refrate, &fwd1);
    if (err)
      goto FREE_RETURN;

    // Calculates the forward
    err = swp_f_ForwardRate(cms2_start_date, cms2_theo_end_date, cms2_frequency,
                            cms2_basis, cms2_yc_name, cms2_refrate, &fwd2);
    if (err)
      goto FREE_RETURN;

    err = digital_quanto_float_leg_cms_rate_quanto(
        fwd1, cms1_fixing_time, cms1_start_time, cms1_nfp, cms1_compounding,
        cms1_delay, cms1_dRateConv, cms1_vol_type, cms_flat_vol_or_smile_approx,
        cms1_start_date, cms1_theo_end_date, cms1_cash_vol, fwd_spread1,
        cms1_vc_name, n_cms1_strikes_in_vol, cms1_strikes_in_vol, cms1_dom_for,
        cms1_quanto_corr_times, cms1_quanto_corr_values, n_cms1_quanto_corr,
        fwd_fx_vols_times, fwd_fx_vols_values, n_fwd_fx_vols, &cms1);
    if (err)
      goto FREE_RETURN;

    err = digital_quanto_float_leg_cms_rate_quanto(
        fwd2, cms2_fixing_time, cms2_start_time, cms2_nfp, cms2_compounding,
        cms2_delay, cms2_dRateConv, cms2_vol_type, cms_flat_vol_or_smile_approx,
        cms2_start_date, cms2_theo_end_date, cms2_cash_vol, fwd_spread2,
        cms2_vc_name, n_cms2_strikes_in_vol, cms2_strikes_in_vol, cms2_dom_for,
        cms2_quanto_corr_times, cms2_quanto_corr_values, n_cms2_quanto_corr,
        fwd_fx_vols_times, fwd_fx_vols_values, n_fwd_fx_vols, &cms2);
    if (err)
      goto FREE_RETURN;

    // Prices the exact option (B-CMS1)+*(CMS2+m)

    // CMS1 ATS B log vol
    err = digital_quanto_float_leg_get_lognormal_vol(
        cms1_vc_name, fwd1, fwd_spread1, cms1_fixing_time, cms1_start_date,
        cms1_theo_end_date, exo_barriers[i], &cms1_ats_log_vol);
    if (err)
      goto FREE_RETURN;

    if (exo_gearings[i] != 0) {
      // Numeraire adjustment ne prend que corr_bid pour le moment

      err = constant_or_linear_interpolation(
          cms1_cms2_corr_times, cms1_cms2_corr_values_bid, n_cms1_cms2_corr,
          LINEAR_METHOD, min(cms1_start_time, cms2_start_time),
          &cms1_cms2_correlation);
      if (err)
        goto FREE_RETURN;

      if (atm_ats_std == 0) {
        // CMS2 vol is taken at the money
        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name, fwd2, fwd_spread2, cms2_fixing_time, cms2_start_date,
            cms2_theo_end_date, fwd2, &cms2_log_vol);
        if (err)
          goto FREE_RETURN;
      } else if (atm_ats_std == 1) {
        // CMS2 vol is taken at the strike
        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name, fwd2, fwd_spread2, cms2_fixing_time, cms2_start_date,
            cms2_theo_end_date, exo_barriers[i], &cms2_log_vol);
        if (err)
          goto FREE_RETURN;
      } else if (atm_ats_std == 2) {
        // CMS2 vol is taken + nCMS1*CMS2logvol*sqrt(CMS2_fixing_time)
        // where nCMS1 is the number of lognormal stdev between the strike
        // exo_barriers[i] and the fwd1

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms1_vc_name, fwd1, fwd_spread1, cms1_fixing_time, cms1_start_date,
            cms1_theo_end_date, fwd1, &cms1_atm_log_vol);
        if (err)
          goto FREE_RETURN;

        cms1_nb_std = log(exo_barriers[i] / fwd1) /
                      (cms1_atm_log_vol * sqrt(cms1_fixing_time));

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name, fwd2, fwd_spread2, cms2_fixing_time, cms2_start_date,
            cms2_theo_end_date, fwd2, &cms2_atm_log_vol);
        if (err)
          goto FREE_RETURN;

        cms2_strike = fwd2 * exp(cms1_cms2_correlation * cms2_atm_log_vol *
                                 cms1_nb_std * sqrt(cms2_fixing_time));

        err = digital_quanto_float_leg_get_lognormal_vol(
            cms2_vc_name, fwd2, fwd_spread2, cms2_fixing_time, cms2_start_date,
            cms2_theo_end_date, cms2_strike, &cms2_log_vol);
        if (err)
          goto FREE_RETURN;
      }

      numeraire_adjustment =
          exp(cms1_cms2_correlation * cms1_ats_log_vol * cms2_log_vol *
              min(cms1_fixing_time, cms2_fixing_time));
    }

    // Calculation of the exact option (strike B)
    put_B = exo_gearings[i] * cms2 *
            srt_f_optblksch(cms1 * numeraire_adjustment, exo_barriers[i],
                            cms1_ats_log_vol, cms1_fixing_time, 1.0, SRT_PUT,
                            SRT_PREMIUM);

    put_margin_B = exo_margins[i] *
                   srt_f_optblksch(cms1, exo_barriers[i], cms1_ats_log_vol,
                                   cms1_fixing_time, 1.0, SRT_PUT, SRT_PREMIUM);

    // Calculation of the strike for the shifted option (B+dB or B-dB)
    if (above_below == 1) {
      if (pay_rec == 0)
        shifted_strike = exo_barriers[i] + call_spread;
      else
        shifted_strike = exo_barriers[i] - call_spread;
    } else {
      if (pay_rec == 0)
        shifted_strike = exo_barriers[i] - call_spread;
      else
        shifted_strike = exo_barriers[i] + call_spread;
    }

    // CMS1 ATS shifted B log vol

    err = digital_quanto_float_leg_get_lognormal_vol(
        cms1_vc_name, fwd1, fwd_spread1, cms1_fixing_time, cms1_start_date,
        cms1_theo_end_date, shifted_strike, &cms1_ats_log_vol);
    if (err)
      goto FREE_RETURN;

    // Pricing of the puts at shifted_strike
    put_B_shift = exo_gearings[i] * cms2 *
                  srt_f_optblksch(cms1 * numeraire_adjustment, shifted_strike,
                                  cms1_ats_log_vol, cms1_fixing_time, 1.0,
                                  SRT_PUT, SRT_PREMIUM);

    put_margin_B_shift =
        exo_margins[i] * srt_f_optblksch(cms1, shifted_strike, cms1_ats_log_vol,
                                         cms1_fixing_time, 1.0, SRT_PUT,
                                         SRT_PREMIUM);

    // Price the digital quanto float coupon
    if (above_below == 1) // CMS1<B
    {
      if (pay_rec == 0) // PAY
        exotic_coupon_PV =
            -(put_B_shift + put_margin_B_shift - put_B - put_margin_B) /
            call_spread;
      else // REC
        exotic_coupon_PV =
            (put_B + put_margin_B - put_B_shift - put_margin_B_shift) /
            call_spread;
    } else // CMS1>B
    {
      if (pay_rec == 0) // PAY
        exotic_coupon_PV =
            -(exo_gearings[i] * cms2 + exo_margins[i]) +
            (put_B + put_margin_B - put_B_shift - put_margin_B_shift) /
                call_spread;
      else // REC
        exotic_coupon_PV =
            exo_gearings[i] * cms2 + exo_margins[i] -
            (put_B_shift + put_margin_B_shift - put_B - put_margin_B) /
                call_spread;
    }

    // Parsing the exotic leg basis
    err = interp_basis(exo_basis[i], &srt_exo_basis_code);
    if (err)
      goto FREE_RETURN;

    cvg = coverage(exo_start_dates[i], exo_end_dates[i], srt_exo_basis_code);
    df_pay = swp_f_df(today, exo_pay_dates[i], dom_yc_name);
    exotic_coupon_PV *= df_pay * cvg * exo_notionals[i];

    // Adds the coupon to the sum
    exotic_leg_PV += exotic_coupon_PV;

    // Increment the coupon and continue
    i++;
  }

// Returns the final result
FREE_RETURN:
  *DigitalQuantoFloatLeg = exotic_leg_PV;
  return err;
}
