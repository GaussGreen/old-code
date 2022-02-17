#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srtaccess.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "swp_h_vol.h"

#include "CTSQTOProdStruct.h"
#include "FundLegProdStruct.h"
#include "LGMQuantoUnd.h"
#include "RangeAccrualProdStruct.h"

static int NUM_HERMITE = 10;

static double dens(double x) { return INV_SQRT_TWO_PI * exp(-x * x / 2.0); }

Err ctsqto_caller(
    /*	Today's date */
    long today,
    /*	The underlying */
    int use_calib, /*	0: use lgm2fund        , 1: calibrate */
    /*		if calib */
    char *dom_yc,         /*	dom yc */
    char *dom_vc,         /*	dom vc */
    char *dom_ref,        /*	dom ref rate */
    char *dom_swap_freq,  /*	dom swap freq */
    char *dom_swap_basis, /*	dom swap basis */
    double dom_lambda,    /*	dom lambda if unique */

    char *for_yc,          /*	for yc */
    char *for_vc,          /*	for vc */
    char *for_ref,         /*	for ref rate */
    char *for_instr_freq,  /*	for instr freq */
    char *for_instr_basis, /*	for instr basis */
    double for_lambda,     /*	for lambda if unique */
    int forcalib,          /*	0 : RA Und        , 1 : Diag */

    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*		if no calilb */
    char *fx3dund,
    /*	The structure */
    long start_date, /*	Date at which initial notional exchange occurs */
    /*		funding */
    int fund_ccy,    /*	0: domestic        , 1: other */
    double fund_not, /*	If different from domestic or foreign (fund_ccy = 1) */
    char *
        fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 1) */
    double fx_fund_dom, /*	If different from domestic or foreign (fund_ccy
                     = 1) 2 bd fwd */
    long fx_fund_dom_spot_date, int fund_ncpn, long *fund_fix, long *fund_start,
    long *fund_pay, char **fund_basis, double *fund_spr, double *fund_mrg,
    double *fund_fix_cpn, /*	Past coupon fixing if relevant        ,
                            includes spr        , but not mrg        , cvg and
                            notional */
    /*		ra */
    double ra_not, int ra_cpn_type, int ra_ncpn, double *ra_cpns,
    long *ra_start, long *ra_pay, char *ra_refrate, char *ra_basis,
    int *ra_nfixings, long **ra_fixingdates,
    double **ra_fixings, /*	Past coupon fixing if relevant */

    // RA floating coupons
    int ra_float_refrate_is_dom_for, char *ra_float_refrate,
    long *ra_float_startl, double *ra_float_past_fixings,
    double *ra_float_gearings,

    double *upper_barr, double *lower_barr, double c_spread, long obs_freq_iv,
    long obs_freq_op, double rho_df,

    // Params for numeraire adjustment if a floating coupon is paid
    double correl_start, double correl_end, int float_adj_strike,

    int n_fxvol_dates, long *fxvol_dates, double *fxvol, double *qtocorrel,

    int typeVol,

    double *corr_times, double *correl_dom_for, double *correl_dom_fx,
    double *correl_for_fx, int corr_n_times,

    /*		calls */
    int ncall, int pay_rec, /*	0: rec exo        , 1: pay exo */
    long *ex_date, long *set_date, double *fee,
    /*	Numerical params */
    int req_stp, int req_stpx,
    /*	Calib params */
    int dom_force_atm, int for_force_atm, double max_std_long,
    double max_std_short, int fix_lambda, /*	0: calib lambda to cap        ,
                                 1: fix lambda calib to diagonal */
    int one_f_equi,                       /*	1F equivalent flag:
                                                                    if set to 1        , then 2F
                                       lambda will calibrate                       to the cap priced
                                       within calibrated                       1F                       with the given
                                       lambda */
    int skip_last, /*	If 1        , the last option is disregarded
                                             and the forward volatility is
                flat from option n-1 */

    int calc_fwdiv, /* output the model and market fwdiv  */
    int adjust_fee, /* adjust the fee to price correctly the fwd iv */

    /*	EOD Flags */
    int eod_fix_flag, /*	0: I        , 1: E */
    int eod_pay_flag, /*	0: I        , 1: E */
    int eod_ex_flag,  /*	0: I        , 1: E */

    /*	Exercised flag */
    int exercised,    /*	Flag */
    long ex_date_ex,  /*	Date when exercised */
    long ex_date_set, /*	Corresponding settlement date */
    double ex_fee,    /*	Corresponding fee */
    /*	Results */
    double *fund_val, /*	Value of the funding leg */
    double *ra_val,   /*	Value of the Range Accrual leg */
    double *call_val, /*	Value of the callable feature */
    int export_ts,    /*	1: Export TS        , 0: don't */
    LGMQTO_UND und_exp) {
  ctsqto_str *ctsqto = NULL;
  lgmQto_und *und = NULL;
  lgmQto_adi_arg *adi_arg = NULL;

  FUNDING_LEG fund_leg;
  FUNDING_CPN fund_cpn;
  RangeAccrualStruct *ra_leg = NULL;
  RangeAccrualStruct *ra_leg_fixing = NULL;

  int call_feat;
  double fund_leg_pv, ra_leg_pv, ra_leg_pv_fixing;
  int nb_fixing;
  long *all_fixings;
  SrtBasisCode bas;
  char ra_recpay[256];

  double call;
  int i, j;
  int free_struct = 0;

  double temp;

  double ra_cpn_float;

  CTSQTO_CALL callt;
  double fact;

  Err err = NULL;

  long *dates = NULL;

  double *fwdPV = NULL;
  int nbFwdPVDates;
  long *FwdPVDates = NULL;

  int for_fund;
  long fund_start_date, fin_not_date;
  double eq_final_ex, eq_init_ex;
  int ra_spot_lag;
  int rabasis, rafreq;
  char *rafreqchar = NULL;
  int fund_spot_lag;
  int spotlag_f;
  long dates_fixing[2];

  // In case of broken period with a DRS floating coupon
  long d_cpn, d1_cpn, d2_cpn, fixing_cpn_float;
  SrtBasisCode basis_d, basis_float_cpn;
  SrtCompounding freq_d, freq_float_cpn;
  int tenor_float_cpn, spotlag, q;
  char *yc_cpn_float, *volc_cpn_float;
  double vol_fra_cpn, power, vol_fx_cpn, rho_ffx_cpn, vol_d_cpn, coef_cpn,
      quanto_corr_cpn, DRS_corr_cpn, fra_d_cpn, cvg_d_cpn;

  double weight_fixing;

  double *fund_mrg2 = NULL, *fund_spr2 = NULL, *fund_fix_cpn2 = NULL;

  if (pay_rec == 1) {
    strcpy(ra_recpay, "REC");
    fact = -1.0;
  } else {
    strcpy(ra_recpay, "PAY");
    fact = 1.0;
  }

  err = srt_f_get_spot_lag_from_refrate(dom_ref, &fund_spot_lag);

  err = swp_f_get_ref_rate_details(ra_refrate, &rabasis, &rafreq);
  if (err) {
    smessage("Error in swp_f_get_ref_rate_details");
    goto FREE_RETURN;
  }
  rafreqchar = (char *)calloc(256, sizeof(char));
  if (!rafreqchar) {
    smessage("Memory Allocation Failed in ctsqto_caller (1)");
    err = "Memory Allocation Failed in ctsqto_caller (1)";
    goto FREE_RETURN;
  }
  if (rafreq == 1) {
    strcpy(rafreqchar, "A");
  } else if (rafreq == 2) {
    strcpy(rafreqchar, "S");
  } else if (rafreq == 4) {
    strcpy(rafreqchar, "Q");
  } else if (rafreq == 12) {
    strcpy(rafreqchar, "M");
  }

  err = srt_f_get_spot_lag_from_refrate(ra_refrate, &ra_spot_lag);
  if (err) {
    smessage("Error in srt_f_get_spot_lag_from_refrate");
    goto FREE_RETURN;
  }

  /*	If exercised */
  if (exercised) {
    i = 0;
    while (i < ra_ncpn && ra_start[i] < ex_date_ex) {
      i++;
    }
    ra_ncpn = i;

    /*	Structure is called before start: return 0 */
    if (ra_ncpn == 0) {
      *fund_val = *ra_val = *call_val = 0.0;
      return NULL;
    }

    i = 0;
    while (i < fund_ncpn && fund_start[i] < ex_date_ex) {
      i++;
    }
    fund_ncpn = i;

    ncall = 0;
    exercised = 0;

    err = ctsqto_caller(
        today, use_calib,

        dom_yc, dom_vc, dom_ref, dom_swap_freq, dom_swap_basis, dom_lambda,

        for_yc, for_vc, for_ref, for_instr_freq, for_instr_basis, for_lambda,
        forcalib,

        get_cash_vol, fx3dund,

        start_date,

        fund_ccy, fund_not, fund_ccy_yc, fx_fund_dom, fx_fund_dom_spot_date,
        fund_ncpn, fund_fix, fund_start, fund_pay, fund_basis, fund_spr,
        fund_mrg, fund_fix_cpn,

        ra_not, ra_cpn_type, ra_ncpn, ra_cpns, ra_start, ra_pay, ra_refrate,
        ra_basis, ra_nfixings, ra_fixingdates, ra_fixings,

        // RA floating coupons
        ra_float_refrate_is_dom_for, ra_float_refrate, ra_float_startl,
        ra_float_past_fixings, ra_float_gearings,

        upper_barr, lower_barr, c_spread, obs_freq_iv, obs_freq_op, rho_df,

        correl_start, correl_end, float_adj_strike,

        n_fxvol_dates, fxvol_dates, fxvol, qtocorrel,

        typeVol,

        corr_times, correl_dom_for, correl_dom_fx, correl_for_fx, corr_n_times,

        ncall, pay_rec, ex_date, set_date, fee, req_stp, req_stpx,

        dom_force_atm, for_force_atm, max_std_long, max_std_short, fix_lambda,
        one_f_equi, skip_last, calc_fwdiv, adjust_fee, eod_fix_flag,
        eod_pay_flag, eod_ex_flag,

        exercised, ex_date_ex, ex_date_set, ex_fee, fund_val, ra_val, call_val,
        export_ts, und_exp);

    if (err) {
      smessage("Error in ctsqto_caller");
      goto FREE_RETURN;
    }

    if (ex_date_set >= today + eod_pay_flag) {
      *call_val = -ex_fee * swp_f_df(today, ex_date_set, dom_yc);
    }

    goto FREE_RETURN;
  }

  /* save the initial fund margins */
  fund_mrg2 = calloc(fund_ncpn, sizeof(double));
  fund_spr2 = calloc(fund_ncpn, sizeof(double));
  fund_fix_cpn2 = calloc(fund_ncpn, sizeof(double));

  if (!fund_mrg2 || !fund_spr2 || !fund_fix_cpn2) {
    err = "Memory allocation error in cts_QTO_caller";
    goto FREE_RETURN;
  }

  memcpy(fund_mrg2, fund_mrg, fund_ncpn * sizeof(double));
  memcpy(fund_spr2, fund_spr, fund_ncpn * sizeof(double));
  memcpy(fund_fix_cpn2, fund_fix_cpn, fund_ncpn * sizeof(double));

  if (fund_ccy == 1) {
    fund_ccy = 0;
    for_fund = 1;
    err = convert_funding_to_domestic(
        today, start_date, eod_fix_flag, eod_pay_flag, fx_fund_dom,
        fx_fund_dom_spot_date, ra_not, dom_yc, fund_ncpn, fund_fix, fund_start,
        fund_pay, fund_basis, fund_ccy_yc, &fund_not, fund_spr2, fund_mrg2,
        fund_fix_cpn2, &fund_start_date, &eq_final_ex, &eq_init_ex);
    if (err) {
      smessage("Error in convert_funding_to_domestic");
      return err;
    }
  } else {
    for_fund = 0;
  }

  /*	Initialise structures */
  dates = lvector(0, ra_ncpn);
  if (!dates) {
    smessage("Memory allocation dates");
    err = "Memory allocation failed in CTSQTO Caller";
    goto FREE_RETURN;
  }

  dates[0] = ra_start[0];
  for (i = 0; i < ra_ncpn; ++i) {
    dates[i + 1] = ra_pay[i];
  }

  ctsqto = calloc(1, sizeof(ctsqto_str));
  und = calloc(1, sizeof(lgmQto_und));
  adi_arg = calloc(1, sizeof(lgmQto_adi_arg));

  if (!ctsqto || !und || !adi_arg) {
    smessage("Memory allocation ctsqto        , und or adi_arg");
    err = "memory allocation failure in ctsqto_caller";
    goto FREE_RETURN;
  }

  free_struct = 0;

  err = ctsqto_fill_check_all_struct(
      today, use_calib, dom_yc, dom_vc, dom_ref,
      //					rafreqchar        ,
      dom_swap_freq, dom_swap_basis, dom_lambda,

      for_yc, for_vc, for_ref, rafreqchar, for_instr_basis, for_lambda,
      forcalib, get_cash_vol,

      fx3dund, fund_not, fund_ncpn, fund_fix, fund_start, fund_pay, fund_basis,
      fund_spr2, fund_mrg2,

      ra_not, ra_cpn_type, n_fxvol_dates, fxvol_dates, fxvol, qtocorrel,
      typeVol,

      ra_ncpn, dates, ra_cpns, ra_basis, ra_nfixings, ra_fixingdates,
      ra_fixings, /*	Past coupon fixing if relevant */

      // RA floating coupons
      ra_float_refrate_is_dom_for, ra_float_refrate, ra_float_startl,
      ra_float_past_fixings, ra_float_gearings,

      upper_barr, lower_barr, ra_recpay, ra_refrate, c_spread, obs_freq_op,
      rho_df,

      correl_start, correl_end, float_adj_strike,

      ncall, pay_rec, ex_date, set_date, fee,

      req_stp, req_stpx,

      dom_force_atm, for_force_atm, max_std_long, max_std_short, fix_lambda,
      one_f_equi, skip_last,

      fxvol_dates, fxvol, n_fxvol_dates, corr_times, correl_dom_for,
      correl_dom_fx, correl_for_fx, corr_n_times,

      eod_fix_flag, eod_ex_flag,

      fund_spot_lag, ra_spot_lag,

      ctsqto, und,

      &call_feat,

      adi_arg);

  smessage("ctsqto_fill_check_all_struct OK");

  if (err) {
    smessage("Error in ctsqto_fill_check_all_struct");
    goto FREE_RETURN;
  }
  free_struct = 1;

  /*  0) calculate the fwd iv in the model */

  if (calc_fwdiv && (ctsqto->num_calls > 0)) {
    und->has_fwd_iv = 1;
    und->nb_fwdiv = ctsqto->num_calls;

    und->exercise_date = (double *)calloc(und->nb_fwdiv, sizeof(double));
    und->market_fwdiv = (double *)calloc(und->nb_fwdiv, sizeof(double));
    und->model_fwdiv = (double *)calloc(und->nb_fwdiv, sizeof(double));
    und->extra_fees = (double *)calloc(und->nb_fwdiv, sizeof(double));

    if (!und->exercise_date || !und->market_fwdiv || !und->model_fwdiv ||
        !und->extra_fees) {
      smessage("Memory allocation und->exercise_date");
      err = "Memory allocation faillure in ctsqto caller";
      goto FREE_RETURN;
    }

    //		NUM_HERMITE = req_stpx;
    err = ctsqto_calc_mdl_iv_fwd(ctsqto, und, adi_arg, NUM_HERMITE,
                                 und->model_fwdiv);

    smessage("ctsqto_calc_mdl_iv_fwd OK");
    if (err) {
      smessage("Error in ctsqto_calc_mdl_iv_fwd");
      goto FREE_RETURN;
    }

    for (i = 0; i < ctsqto->num_calls; i++) {
      callt = ctsqto->call + i;
      und->exercise_date[i] = callt->ex_date;
    }
  }

  /*	1)	Value funding leg */

  fund_leg = ctsqto->fund_leg;
  fund_leg_pv = 0.0;

  smessage("funding leg");

  /*	Cash libor */
  if (fund_leg->num_cpn > 0) {
    if (for_fund) {
      fund_leg_pv +=
          swp_f_df(today, fund_leg->cpn[0].start_date, dom_yc) * eq_final_ex;

      if (calc_fwdiv) {
        // initialisation of the initial notional
        for (i = 0; i < ctsqto->num_calls; i++) {
          callt = ctsqto->call + i;
          und->market_fwdiv[i] =
              fact *
              swp_f_df(today,
                       (ctsqto->fund_leg->cpn + callt->fund_idx)->start_date,
                       dom_yc) *
              eq_final_ex;
        }
      }
    } else {
      fund_leg_pv += swp_f_df(today, fund_leg->cpn[0].start_date, dom_yc) *
                     fund_leg->notional;

      if (calc_fwdiv) {
        // initialisation of the initial notional
        for (i = 0; i < ctsqto->num_calls; i++) {
          callt = ctsqto->call + i;
          und->market_fwdiv[i] =
              fact *
              swp_f_df(today,
                       (ctsqto->fund_leg->cpn + callt->fund_idx)->start_date,
                       dom_yc) *
              fund_leg->notional;
        }
      }
    }
  } else {
    fin_not_date = fund_pay[fund_ncpn - 1];
    if (fin_not_date >= today + eod_pay_flag) {
      if (for_fund) {
        temp = swp_f_df(today, fin_not_date, dom_yc) * eq_final_ex;
        fund_leg_pv += temp;

        if (calc_fwdiv) {
          for (i = 0; i < ctsqto->num_calls; i++) {
            und->market_fwdiv[i] += fact * temp;
          }
        }
      }
    }
  }

  smessage("funding leg OK");

  smessage("spread + margin");
  /*	Coupons: spread + margin */
  for (i = 0; i < fund_leg->num_cpn; i++) {
    fund_cpn = fund_leg->cpn + i;
    temp = swp_f_df(today, fund_cpn->pay_date, dom_yc) * fund_cpn->cpn;
    fund_leg_pv += temp;
    //		smessage("for spread");

    if (calc_fwdiv) {
      j = 0;
      if (j < ctsqto->num_calls) {
        //				smessage("test spread");
        while ((j < ctsqto->num_calls) && ((ctsqto->call + j)->fund_idx <= i)) {
          und->market_fwdiv[j] += fact * temp;
          //					smessage("while spread");
          j++;
        }
      }
    }
  }

  smessage("spread + margin OK");

  smessage("Notional Exchage");
  /*	Notional exchange */
  if (for_fund) {
    if (start_date >= today + eod_pay_flag) {
      fund_leg_pv -= swp_f_df(today, start_date, dom_yc) * eq_init_ex;
    }

    if (calc_fwdiv) {
      for (i = 0; i < ctsqto->num_calls; i++) {
        callt = ctsqto->call + i;
        und->market_fwdiv[i] -=
            fact *
            swp_f_df(today,
                     (ctsqto->fund_leg->cpn + callt->fund_idx)->start_date,
                     dom_yc) *
            eq_final_ex;
      }
    }
  } else {
    if (fund_leg->num_cpn > 0) {
      temp = swp_f_df(today, fund_leg->cpn[fund_leg->num_cpn - 1].pay_date,
                      dom_yc) *
             fund_leg->notional;
      fund_leg_pv -= temp;
      if (calc_fwdiv) {
        for (i = 0; i < ctsqto->num_calls; i++) {
          und->market_fwdiv[i] -= fact * temp;
        }
      }
    }
  }

  smessage("Notional Exchage OK");

  smessage("Past Fixings");
  /*	PV of coupons fixed in the past and not yet paid */
  i = 0;
  //	while (i < fund_ncpn && add_unit( fund_start[i]        , -spot_lag ,
  //						SRT_BDAY        , MODIFIED_SUCCEEDING ) < today
  //+ eod_fix_flag)
  while (i < fund_ncpn && fund_fix[i] < today + eod_fix_flag) {
    if (fund_pay[i] >= today + eod_pay_flag) {
      err = interp_basis(fund_basis[i], &bas);
      if (err) {
        smessage("Error in interp_basis");
        goto FREE_RETURN;
      }

      fund_leg_pv += (fund_fix_cpn2[i] + fund_mrg2[i]) *
                     coverage(fund_start[i], fund_pay[i], bas) * fund_not *
                     swp_f_df(today, fund_pay[i], dom_yc);
    }

    i++;
  }

  smessage("Past Fixings OK");

  /*	2)	Value ra leg */

  smessage("RA alloc");
  ra_leg = (RangeAccrualStruct *)calloc(1, sizeof(RangeAccrualStruct));
  if (!ra_leg) {
    err = "Memory Allocation Failed ra_leg";
    smessage("RA alloc OK");
  }

  //	ra_leg = ctsqto->ra_leg;
  ra_leg_pv = 0.0;

  smessage("start ra_init_struct");
  err =
      ra_init_struct(today, dom_yc, for_yc, dom_vc, for_vc, dom_ref, for_ref,
                     ra_not, ra_cpn_type, n_fxvol_dates, fxvol_dates, fxvol,
                     qtocorrel, typeVol, ra_ncpn, dates, ra_cpns, ra_basis,
                     ra_nfixings, ra_fixingdates, ra_fixings,

                     // RA float coupons
                     ra_float_refrate_is_dom_for, ra_float_refrate,
                     ra_float_startl, ra_float_past_fixings, ra_float_gearings,

                     upper_barr, lower_barr, ra_recpay, ra_refrate, c_spread,
                     obs_freq_iv, rho_df,

                     // Params for the floating coupon
                     correl_start, correl_end, float_adj_strike,

                     eod_fix_flag, ra_leg);

  smessage("ra_init_struct OK");
  if (err) {
    smessage("Error in ra_init_struct");
    goto FREE_RETURN;
  }

  nbFwdPVDates = ctsqto->num_calls + 1;

  fwdPV = dvector(0, ctsqto->num_calls);
  if (!fwdPV) {
    smessage("Memory allocation error fwdPV");
    err = "Memory allocationfailed in CTSQTOCaller";
    goto FREE_RETURN;
  }

  FwdPVDates = calloc(nbFwdPVDates, sizeof(long));
  if (!FwdPVDates) {
    smessage("Memory allocation error fwdPVDates");
    err = "Memory allocationfailed in CTSQTOCaller";
    goto FREE_RETURN;
  }
  FwdPVDates[0] = today;
  for (i = 0; i < ctsqto->num_calls; ++i) {
    FwdPVDates[i + 1] = ctsqto->call[i].ex_date;
  }

  smessage("RA_FwdPV");
  err =
      RA_FwdPV(dom_yc, for_yc, ra_leg, today, nbFwdPVDates, FwdPVDates, fwdPV);
  if (err) {
    smessage("Error in RA_FwdPV");
    goto FREE_RETURN;
  }

  ra_leg_pv = -fwdPV[0];

  /*	PV of coupons fixed in the past and not yet paid */
  i = 0;
  err = interp_basis(ra_basis, &bas);
  if (err)
    goto FREE_RETURN;

  err = srt_f_get_spot_lag_from_refrate(for_ref, &spotlag_f);
  if (err)
    goto FREE_RETURN;

  // while (i < ra_ncpn && add_unit( ra_start[i]        , -ra_spot_lag        ,
  //					SRT_BDAY        , MODIFIED_SUCCEEDING ) < today
  //+ eod_fix_flag)
  // Changed on 12 Jan 2004 following Laurent's request to:
  while (i < ra_ncpn && ra_fixingdates[i][0] < today + eod_fix_flag) {
    if (ra_pay[i] >= today + eod_pay_flag) {
      j = 0;
      q = 0; // Index for quanto correl

      weight_fixing = coverage(ra_start[i], ra_pay[i], bas) /
                      ((ra_pay[i] - add_unit(ra_fixingdates[i][0], spotlag_f,
                                             SRT_BDAY, MODIFIED_SUCCEEDING)) *
                       1.0) *
                      ra_not * swp_f_df(today, ra_pay[i], dom_yc);

      while (j < ra_nfixings[i] &&
             ra_fixingdates[i][j] < today + eod_fix_flag) {
        // In case of floating coupons
        if (ra_cpn_type != 0) {
          if (ra_float_startl[i] <= ra_start[i]) {
            // Floating coupon fixed in advance (already fixed)
            ra_cpn_float = ra_float_past_fixings[i] * ra_float_gearings[i];
          } else {
            // Floating coupon fixed in arrears
            err = swp_f_get_ref_rate_details(ra_float_refrate, &basis_float_cpn,
                                             &freq_float_cpn);
            if (err)
              goto FREE_RETURN;

            tenor_float_cpn = 12 / freq_float_cpn;
            d_cpn = add_unit(ra_float_startl[i], tenor_float_cpn, SRT_MONTH,
                             MODIFIED_SUCCEEDING);

            if (ra_float_refrate_is_dom_for == 0)
              err = srt_f_get_spot_lag_from_refrate(dom_ref, &spotlag);
            else
              err = srt_f_get_spot_lag_from_refrate(for_ref, &spotlag);
            if (err)
              goto FREE_RETURN;

            fixing_cpn_float = add_unit(ra_float_startl[i], -spotlag, SRT_BDAY,
                                        MODIFIED_SUCCEEDING);

            // Get the FRA of the floating coupon
            if (ra_float_refrate_is_dom_for == 0) {
              yc_cpn_float = dom_yc;
              volc_cpn_float = dom_vc;
            } else {
              yc_cpn_float = for_yc;
              volc_cpn_float = for_vc;
            }

            ra_cpn_float = swp_f_fra(ra_float_startl[i], d_cpn, basis_float_cpn,
                                     yc_cpn_float, ra_float_refrate);
            err = swp_f_SABRvol(volc_cpn_float, ra_float_startl[i], d_cpn, 0.05,
                                &vol_fra_cpn, &power, SABR_ATMLOG);
            if (err)
              goto FREE_RETURN;

            // Quanto adjustment of the floating coupon
            if (ra_float_refrate_is_dom_for == 1) {
              // Get the FX ATM vol & FX correlations (using the same ones as
              // for the quanto index
              for (; q < n_fxvol_dates - 1 &&
                     fxvol_dates[q + 1] < ra_float_startl[i];
                   q++)
                ;
              if (fxvol_dates[q] > d_cpn || q == n_fxvol_dates - 1) {
                vol_fx_cpn = fxvol[q];
                rho_ffx_cpn = qtocorrel[q];
              } else {
                coef_cpn = ((double)(d_cpn - fxvol_dates[q])) /
                           (fxvol_dates[q + 1] - fxvol_dates[q]);
                vol_fx_cpn =
                    coef_cpn * fxvol[q + 1] + (1.0 - coef_cpn) * fxvol[q];
                rho_ffx_cpn = coef_cpn * qtocorrel[q + 1] +
                              (1.0 - coef_cpn) * qtocorrel[q];
              }

              quanto_corr_cpn = -rho_ffx_cpn * vol_fx_cpn * vol_fra_cpn;
            } else {
              quanto_corr_cpn = 0.0;
            }

            // DRS adjustment of the floating coupon
            if (ra_float_startl[i] == ra_start[i] && d_cpn == ra_start[i + 1])
              DRS_corr_cpn = 0;
            else {
              d1_cpn = (ra_start[i + 1] < d_cpn ? ra_start[i + 1] : d_cpn);
              d2_cpn = (ra_start[i + 1] > d_cpn ? ra_start[i + 1] : d_cpn);

              err = swp_f_get_ref_rate_details(dom_ref, &basis_d, &freq_d);
              if (err)
                goto FREE_RETURN;

              if (fabs(d2_cpn - d1_cpn) < 10) {
                d2_cpn = d1_cpn + 10;
              }
              fra_d_cpn = swp_f_fra(d1_cpn, d2_cpn, basis_d, dom_yc, dom_ref);
              cvg_d_cpn = coverage(d1_cpn, d2_cpn, basis_d);
              err = swp_f_SABRvol(dom_vc, d1_cpn, d2_cpn, fra_d_cpn, &vol_d_cpn,
                                  &power, SABR_LOGVOL);
              if (err)
                return err;

              DRS_corr_cpn = rho_df * vol_d_cpn * vol_fra_cpn * cvg_d_cpn *
                             fra_d_cpn / (1.0 + cvg_d_cpn * fra_d_cpn);
              if (ra_start[i + 1] > d_cpn)
                DRS_corr_cpn = -DRS_corr_cpn;
            }

            // DRS floating coupon is...
            ra_cpn_float *=
                ra_float_gearings[i] * exp((DRS_corr_cpn + quanto_corr_cpn) *
                                           (fixing_cpn_float - today) / 365.0);
          }
        } else {
          ra_cpn_float = 0.0;
        }

        if ((ra_fixings[i][j] > lower_barr[i]) &&
            (ra_fixings[i][j] < upper_barr[i]))
        // For testing: if( (ra_fixings[i][j] > -0.1) && (ra_fixings[i][j] >
        // lower_barr[i]) && (ra_fixings[i][j] < upper_barr[i]) )
        {
          if (j < ra_nfixings[i] - 1) {
            ra_leg_pv -= (add_unit(ra_fixingdates[i][j + 1], spotlag_f,
                                   SRT_BDAY, MODIFIED_SUCCEEDING) -
                          add_unit(ra_fixingdates[i][j], spotlag_f, SRT_BDAY,
                                   MODIFIED_SUCCEEDING)) *
                         weight_fixing * (ra_cpns[i] + ra_cpn_float);
          } else {
            ra_leg_pv -= (ra_pay[i] - add_unit(ra_fixingdates[i][j], spotlag_f,
                                               SRT_BDAY, MODIFIED_SUCCEEDING)) *
                         (ra_cpns[i] + ra_cpn_float) * weight_fixing *
                         (ra_cpns[i] + ra_cpn_float);
          }
        }

        j++;
      }

      if (j < ra_nfixings[i]) {
        /* remaining part must be valued */
        ra_leg_fixing =
            (RangeAccrualStruct *)calloc(1, sizeof(RangeAccrualStruct));
        if (!ra_leg_fixing) {
          err = "Memory Allocation Failed ra_leg";
          goto FREE_RETURN;
        }

        ra_leg_pv_fixing = 0.0;

        nb_fixing = ra_nfixings[i] - j;
        all_fixings = &(ra_fixingdates[i][j]);

        dates_fixing[0] = add_unit(ra_fixingdates[i][j], spotlag_f, SRT_BDAY,
                                   MODIFIED_SUCCEEDING);
        dates_fixing[1] = ra_pay[i];

        smessage("start ra_init_struct");

        err = ra_init_struct(today, dom_yc, for_yc, dom_vc, for_vc, dom_ref,
                             for_ref, ra_not, ra_cpn_type, n_fxvol_dates,
                             fxvol_dates, fxvol, qtocorrel, typeVol, 1,
                             &(dates_fixing[0]), &(ra_cpns[i]), ra_basis,
                             &nb_fixing, &all_fixings, &(ra_fixings[i]),

                             // RA float coupons
                             ra_float_refrate_is_dom_for, ra_float_refrate,
                             &(ra_float_startl[i]), &(ra_float_past_fixings[i]),
                             &(ra_float_gearings[i]),

                             &(upper_barr[i]), &(lower_barr[i]), ra_recpay,
                             ra_refrate, c_spread, obs_freq_iv, rho_df,

                             // Params for the floating coupon
                             correl_start, correl_end, float_adj_strike,

                             eod_fix_flag, ra_leg_fixing);

        if (err)
          goto FREE_RETURN;

        err = RA_FwdPV(dom_yc, for_yc, ra_leg_fixing, today, 1, FwdPVDates,
                       &ra_leg_pv_fixing);

        /* coverage adjustment */
        ra_leg_pv_fixing *= coverage(ra_start[i], ra_pay[i], bas) /
                            coverage(dates_fixing[0], ra_pay[i], bas);

        /* fixing adjustment */
        ra_leg_pv_fixing *=
            (dates_fixing[1] - dates_fixing[0]) /
            ((ra_pay[i] - add_unit(ra_fixingdates[i][0], spotlag_f, SRT_BDAY,
                                   MODIFIED_SUCCEEDING)) *
             1.0);

        if (err)
          goto FREE_RETURN;

        ra_leg_pv -= ra_leg_pv_fixing;
      }
    }

    i++;
  }

  /*	Initial and final exchange */
  if (for_fund) {
    //	Final
    if (ra_leg->n_periods > 0) {
      temp = swp_f_df(today, fund_leg->cpn[fund_leg->num_cpn - 1].pay_date,
                      dom_yc) *
             ra_not;
      ra_leg_pv += -fact * temp;

      if (calc_fwdiv) {
        for (i = 0; i < ctsqto->num_calls; i++) {
          und->market_fwdiv[i] += -fact * temp;
        }
      }
    } else {
      fin_not_date = fund_pay[fund_ncpn - 1];
      if (fin_not_date >= today + eod_pay_flag) {
        temp = swp_f_df(today, fin_not_date, dom_yc) * ra_not;
        ra_leg_pv += -fact * temp;

        if (calc_fwdiv) {
          for (i = 0; i < ctsqto->num_calls; i++) {
            und->market_fwdiv[i] += -fact * temp;
          }
        }
      }
    }

    //	Initial
    if (start_date >= today + eod_pay_flag) {
      ra_leg_pv -= -fact * swp_f_df(today, start_date, dom_yc) * ra_not;
    }

    if (calc_fwdiv) {
      for (i = 0; i < ctsqto->num_calls; i++) {
        callt = ctsqto->call + i;
        und->market_fwdiv[i] -=
            -fact *
            swp_f_df(today, ctsqto->ra_leg->dates[callt->ra_idx], dom_yc) *
            ra_not;
      }
    }
  }

  if (calc_fwdiv) {
    for (i = 0; i < ctsqto->num_calls; i++) {
      callt = ctsqto->call + i;
      und->market_fwdiv[i] += -fact * fwdPV[i + 1];
      // Was - (und-> replaced by +) to be consistent with convention for fee
      // (fee is PAID)
      und->extra_fees[i] = +(und->market_fwdiv[i] + und->model_fwdiv[i]) /
                           swp_f_df(today, callt->set_date, dom_yc);

      if (adjust_fee) {
        callt->fee += und->extra_fees[i];
      }
    }
  }

  /*	4)	If there is at least one call after today        , value call
   * feature
   */

  if (call_feat == 1) {
    smessage(
        "Launching adi        , time steps requested: %d        , actual: %d",
        req_stp, adi_arg->nstp);

    err = ctsqto_launch_adi(ctsqto, adi_arg, &call);

    if (err) {
      smessage("Error in ctsqto_launch_adi");
      goto FREE_RETURN;
    }
  } else {
    call = 0.0;
  }

  *fund_val = fund_leg_pv;
  *ra_val = -ra_leg_pv;
  *call_val = call;

  if (export_ts) {
    copy_lgmQto_und(und, und_exp);
  }

FREE_RETURN:

  smessage("Free CTSQTOCaller");
  ra_free(ra_leg);
  if (ra_leg)
    free(ra_leg);

  ra_free(ra_leg_fixing);
  if (ra_leg_fixing)
    free(ra_leg_fixing);

  smessage("Free CTSQTOCaller 1");
  if (rafreqchar)
    free(rafreqchar);

  smessage("Free CTSQTOCaller 2");
  if (dates)
    free_lvector(dates, 0, ra_ncpn);

  smessage("Free CTSQTOCaller 3");
  if (fwdPV)
    free_dvector(fwdPV, 0, ctsqto->num_calls);

  smessage("Free CTSQTOCaller 4");
  if (FwdPVDates)
    free(FwdPVDates);

  smessage("Free CTSQTOCaller 5");
  if (free_struct) {
    ctsqto_free_all_struct(ctsqto, und, call_feat, adi_arg);
  }

  smessage("Free CTSQTOCaller 6");
  if (ctsqto)
    free(ctsqto);

  smessage("Free CTSQTOCaller 6");
  if (und)
    free(und);

  smessage("Free CTSQTOCaller 6");
  if (adi_arg)
    free(adi_arg);

  smessage("Free CTSQTOCaller 6");

  if (fund_mrg2)
    free(fund_mrg2);
  if (fund_spr2)
    free(fund_spr2);
  if (fund_fix_cpn2)
    free(fund_fix_cpn2);

  return err;
}
