/* =========================================================================================
        FILE NAME				SRT_F_FWDVOL.C
        PURPOSE					COMPUTE FORWARD VOLATILITY  ,
   EITHER THROUGH AUTOCALIBRATION  , OR PARABOLAS
   =========================================================================================

   VISION Information Consulting

   Y2K Compliance

   Programmer  : Neil Glover

   Date        : 12/11/1998

   Fixes To    : srt_f_fwdvol_parabola

   Total Fixes : 1

   =========================================================================================
 */

/*	Srt Library */
#include "srt_h_all.h"
/*	For implied vol function */
#include "opfnctns.h"
/*	For underlying initialisation */
#include "srtaccess.h"
/* For CMS constant */
#include "math.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"

#define MAXNUMPERIODS 1200

/*	MEMORY FREE MACROS */

static void free_calib_instr(long *cal_start, long *cal_end, double *cal_str,
                             double *cal_bndstr, String *cal_type,
                             String *cal_recpay, String *cal_refrate,
                             String *cal_freq, String *cal_basis,
                             double *cal_price, double *cal_vega,
                             long cal_numinst, long num_period,
                             int fix_tau_flag) {
  if (cal_start)
    free_lngvector(cal_start, 0, cal_numinst - 1);
  cal_start = NULL;
  if (cal_end)
    free_lngvector(cal_end, 0, cal_numinst - 1);
  cal_end = NULL;
  if (cal_str)
    free_dvector(cal_str, 0, cal_numinst - 1);
  cal_str = NULL;
  if (cal_bndstr)
    free_dvector(cal_bndstr, 0, cal_numinst - 1);
  cal_bndstr = NULL;
  if (cal_type)
    free_svector_size(cal_type, 0, cal_numinst - 1, 32);
  cal_type = NULL;
  if (cal_recpay)
    free_svector_size(cal_recpay, 0, cal_numinst - 1, 32);
  cal_recpay = NULL;
  if (cal_refrate)
    free_svector_size(cal_refrate, 0, cal_numinst - 1, 32);
  cal_refrate = NULL;
  if (cal_price)
    free_dvector(cal_price, 0, cal_numinst - 1);
  cal_price = NULL;
  if (cal_vega)
    free_dvector(cal_vega, 0, cal_numinst - 1);
  cal_vega = NULL;
  if (cal_freq)
    free(cal_freq);
  if (cal_basis)
    free(cal_basis);
}

static void free_option_computation(Date *und_end_th, double *premium_array) {
  if (und_end_th)
    free(und_end_th);
  if (premium_array)
    free(premium_array);
}

#define FREE_GRFN_CALIB_PARAM                                                  \
  {                                                                            \
    free_svector_size(grfn_param, 0, num_grfn_param - 1, 32);                  \
    grfn_param = NULL;                                                         \
    free_svector_size(grfn_value, 0, num_grfn_param - 1, 32);                  \
    grfn_value = NULL;                                                         \
    if (calib_param)                                                           \
      free_svector_size(calib_param, 0, num_calib_param - 1, 32);              \
    calib_param = NULL;                                                        \
    if (calib_value)                                                           \
      free_svector_size(calib_value, 0, num_calib_param - 1, 32);              \
    calib_value = NULL;                                                        \
  }

/* ------------------------------------------------------------------------------------
 */
/*	Main function for autocal option */

Err srt_f_fwdvol_autocal(
    String model,   /* Model (LGM  , CHEY  , CHEY_BETA  , 1F or 2F) */
    long *tau_date, /* Tau Term Struct */
    double *tau_value, long num_tau, int fix_tau_flag,

    double alpha,                              /* 2F parameters */
    double beta, double rho, double chey_beta, /* Chey Beta parameter */
    long num_period, /* Number of required forward vols */
    long *vol_start, /* Vol Start Dates array */
    long *vol_end,   /* Vol End Dates array */

    int option_type, /* Option Type */

    long *und_start_act,                   /* Und start date array */
    char *und_tenor, SrtCompounding compd, /* Und compd */
    SrtBasisCode basis,                    /* Und basi */
    String ref_rate,                       /* Reference Rate Code */
    double *und_val,

    String yc_id,       /* YC id */
    String bs_vol_type, /* Lognormal or normal */
    Err(*GetVol)        /* GetVol function */
    (Ddate start, Ddate end, double strike, double dForward, double dSpread,
     double *bs_vol),

    SrtPriceType price_type, double *output_array, /* Output array */
    String fwd_vol_type)                           /* Lognormal or normal */
{

  /*	---------- Declarations ---------- */
  SrtGrfnParam grfnparam;
  SrtMdlType mdl_type;
  SrtMdlDim mdl_dim;
  Err err;
  SrtUndPtr undptr;
  SrtCurvePtr yldcrv;
  Date today;
  char *basisStr, *compdStr;
  long i;
  SrtDiffusionType srt_bs_vol_type, srt_fwd_vol_type;
  int sig_col, tau_col;
  double **sig_date_value, **tau_date_value;
  long cal_numinst;
  Date *cal_start = NULL;
  Date *cal_end = NULL;
  String *cal_freq = NULL;
  String *cal_basis = NULL;
  String *cal_type = NULL;
  String *cal_recpay = NULL;
  String *cal_refrate = NULL;
  double *cal_str = NULL, *cal_bndstr = NULL, *cal_price = NULL,
         *cal_vega = NULL;
  String *grfn_param = NULL, *grfn_value = NULL, *calib_param = NULL,
         *calib_value = NULL;
  long num_grfn_param, num_calib_param;
  double premium, lvl;
  char und_id[32];
  double chisq;
  long tmpl1 = 0, tmpl2 = 0, tmpl3 = 0, tmpl4 = 0;
  double **tmpsig = NULL, **tmptau = NULL;
  double cmsvol;
  Date *und_end_th = NULL;
  double *premium_array = NULL;
  double SwpSpread;
  char yc_name[32];
  Date clcn_date;
  double *pdTheoprices;
  int *indexUsedInstr;

  rem_tick_string(yc_id, yc_name);
  yldcrv = lookup_curve(yc_id);
  if (!yldcrv) {
    return serror("Error in reset_option  , YC not found: %s", yc_id);
  }

  clcn_date = get_clcndate_from_yldcrv(yldcrv);

  /*	---------- Put information in relevant format ---------- */

  /*	Number of periods */
  if (num_period > MAXNUMPERIODS) {
    return serror("Error in srt_f_fwdvol_autocal (num_period): too large");
  }

  /*	Model */
  err = srt_f_interp_model(model, &mdl_type, &mdl_dim);
  if (err) {
    return serror("Error in srt_f_fwdvol_autocal (srt_f_interp_model): %s",
                  err);
  }
  if ((mdl_dim == TWO_FAC) &&
      (fabs(alpha) < EPS || fabs(beta) < EPS || fabs(rho) < EPS)) {
    return serror("Error in srt_f_fwdvol_autocal: two factor model requires "
                  "alpha  , beta and rho");
  }
  if ((mdl_type == CHEY_BETA) && (fabs(chey_beta) < EPS)) {
    return serror("Error in srt_f_fwdvol_autocal: CHEY_BETA model requires "
                  "chey beta parameter");
  }

  /*	Today */
  yldcrv = lookup_curve(yc_id);
  if (!yldcrv) {
    return serror("Error in srt_f_fwdvol_autocal  , YC not found: %s", yc_id);
  }
  today = get_clcndate_from_yldcrv(yldcrv);

  /*	Lognormal / Normal */
  err = interp_diffusion_type(bs_vol_type, &srt_bs_vol_type);
  err = interp_diffusion_type(fwd_vol_type, &srt_fwd_vol_type);
  if (err) {
    return serror("Error in srt_f_fwdvol_autocal (interp_diffusion_type): %s",
                  err);
  }

  /*	Sigma curve */
  if (mdl_dim == ONE_FAC) {
    sig_col = 2;
    sig_date_value = dmatrix(0, 1, 0, 0);
    if (!sig_date_value) {
      return serror(
          "Error in srt_f_fwdvol_autocal (sig_date_value): memory allocation");
    }
    sig_date_value[0][0] = today + 365;
    if (mdl_type == LGM) {
      sig_date_value[1][0] = 0.01;
    } else if (mdl_type == CHEY) {
      sig_date_value[1][0] = 0.20;
    } else if (mdl_type == CHEY_BETA) {
      sig_date_value[1][0] = 0.01 / pow(0.05, chey_beta);
    }
  } else
  /* TWO FAC */
  {
    sig_col = 4;
    sig_date_value = dmatrix(0, 3, 0, 0);
    if (!sig_date_value) {
      return serror(
          "Error in srt_f_fwdvol_autocal (sig_date_value): memory allocation");
    }
    sig_date_value[0][0] = today + 365;
    sig_date_value[1][0] = 0.05;
    sig_date_value[2][0] = alpha * sig_date_value[1][0];
    sig_date_value[3][0] = rho;
  }

  /*	Tau curve */
  if (mdl_dim == ONE_FAC) {
    tau_col = 2;
    tau_date_value = dmatrix(0, 1, 0, num_tau - 1);
    if (!tau_date_value) {
      free_dmatrix(sig_date_value, 0, sig_col - 1, 0, 0);
      return serror(
          "Error in srt_f_fwdvol_autocal (tau_date_value): memory allocation");
    }
    for (i = 0; i < num_tau; i++) {
      tau_date_value[0][i] = tau_date[i];
      tau_date_value[1][i] = tau_value[i];
    }
  } else
  /* TWO FAC */
  {
    tau_col = 3;
    tau_date_value = dmatrix(0, 2, 0, num_tau - 1);
    if (!tau_date_value) {
      free_dmatrix(sig_date_value, 0, sig_col - 1, 0, 0);
      return serror(
          "Error in srt_f_fwdvol_autocal (tau_date_value): memory allocation");
    }
    for (i = 0; i < num_tau; i++) {
      tau_date_value[0][i] = tau_date[i];
      tau_date_value[1][i] = tau_value[i];
      tau_date_value[2][i] = 1.0 / (1.0 / tau_value[i] + beta);
    }
  }

  /*	---------- Initialise underlying ---------- */

  err = SrtInitIRUnd("TMPUND", yc_id, model, 1, sig_col, sig_date_value,
                     num_tau, tau_col, tau_date_value, chey_beta, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, /* vasicek parms */
                     0, 0, NULL);

  free_dmatrix(sig_date_value, 0, sig_col - 1, 0, 0);
  free_dmatrix(tau_date_value, 0, tau_col - 1, 0, num_tau - 1);

  if (err) {
    return serror("Error in srt_f_fwdvol_autocal (SrtInitIRUnd): %s", err);
  }

  /*	---------- Set calibration instruments ---------- */

  if (option_type == RESETCAPFLOOR) {
    /*	One caplet at each vol start date + one caplet at each end date */
    err = capfwd_vol_set_calibration_instruments(
        num_period, vol_start, vol_end, ref_rate, compd, basis, yc_id,
        bs_vol_type, GetVol, &cal_numinst, &cal_start, &cal_end, &cal_freq,
        &cal_basis, &cal_str, &cal_bndstr, &cal_type, &cal_recpay, &cal_refrate,
        &cal_price, &cal_vega);
  }

  if (option_type == RESETCMSOPTION) {

    err = cmsfwd_vol_set_calibration_instruments(
        num_period, fix_tau_flag, vol_start, vol_end, und_tenor, compd, basis,
        ref_rate, yc_id, bs_vol_type, GetVol, &cal_numinst, &cal_start,
        &cal_end, &cal_freq, &cal_basis, &cal_str, &cal_bndstr, &cal_type,
        &cal_recpay, &cal_refrate, &cal_price, &cal_vega);
  }

  if (err) {
    free_calib_instr(cal_start, cal_end, cal_str, cal_bndstr, cal_type,
                     cal_recpay, cal_refrate, cal_freq, cal_basis, cal_price,
                     cal_vega, cal_numinst, num_period, fix_tau_flag);

    srt_f_destroy_und("TMPUND");
    return serror("Error in srt_f_fwdvol_autocal "
                  "(fwd_vol_set_calibration_instruments): %s",
                  err);
  }

  /*	---------- Set grfn & calib options ---------- */

  /* ------- grfn -------- */
  fwd_vol_set_grfn_options(mdl_type, &grfn_param, &grfn_value, &num_grfn_param);

  /* ------- calib ------- */
  fwd_vol_set_calib_options(mdl_type, option_type, fix_tau_flag, &calib_param,
                            &calib_value, &num_calib_param);

  if (err = srt_f_set_GrfnParams(num_grfn_param, grfn_param, grfn_value,
                                 &grfnparam)) {
    free_calib_instr(cal_start, cal_end, cal_str, cal_bndstr, cal_type,
                     cal_recpay, cal_refrate, cal_freq, cal_basis, cal_price,
                     cal_vega, cal_numinst, num_period, fix_tau_flag);
    FREE_GRFN_CALIB_PARAM;
    srt_f_destroy_und("TMPUND");
    return ("Error in srt_f_fwdvol_autocal (srt_f_set_GrfnParams): %s", err);
  }

  /*	---------- Calibrate ---------- */
  pdTheoprices = dvector(1, cal_numinst);
  indexUsedInstr = ivector(0, cal_numinst - 1);

  err = SrtNewCalibrateAll(
      grfn_param, grfn_value, num_grfn_param, cal_start, cal_end, cal_freq,
      cal_basis, cal_str, cal_bndstr, cal_type, cal_recpay, cal_refrate,
      cal_price, cal_vega, GetVol, &bs_vol_type, cal_numinst, NULL, NULL, 0,
      "TMPUND", und_id, calib_param, calib_value, num_calib_param, &tmpsig,
      &tmpl1, &tmpl2, &tmptau, &tmpl3, &tmpl4, &chisq, &pdTheoprices,
      indexUsedInstr);

  free_dvector(pdTheoprices, 1, cal_numinst);
  pdTheoprices = NULL;
  free_ivector(indexUsedInstr, 0, cal_numinst - 1);
  indexUsedInstr = NULL;

  free_calib_instr(cal_start, cal_end, cal_str, cal_bndstr, cal_type,
                   cal_recpay, cal_refrate, cal_freq, cal_basis, cal_price,
                   cal_vega, cal_numinst, num_period, fix_tau_flag);

  FREE_GRFN_CALIB_PARAM;
  if (err) {
    return serror("Error in srt_f_fwdvol_autocal (srt_f_CalibrateAll): %s",
                  err);
  }
  if (chisq > 0.25) {
    srt_f_destroy_und(und_id);
    return serror("Error in srt_f_fwdvol_autocal (srt_f_CalibrateAll): "
                  "calibration failed");
  }

  premium_array = (double *)malloc(num_period * sizeof(double));

  if (option_type == RESETCAPFLOOR) {

    und_end_th = (Date *)malloc(num_period * sizeof(Date));

    /*	---------- Price forward reset caplets and invert BS to get fwd vol
     * ---------- */
    /*	Forward Reset Caplet is the following profile:
            At T4  , pay: MAX (F(T2  ,T3  ,T4) - F(T1  ,T3  ,T4)  , 0) * cvg (T3
       , T4) Where	T1 = vol start date T2 = vol end date T3 = fra start
       date T4 = fra end date */

    /* Lookup und */
    undptr = lookup_und(und_id);

    /*	Loop on periods */
    for (i = 0; i < num_period; i++) {
      und_end_th[i] = add_unit(DTOL(und_start_act[i]), 12 / (int)compd,
                               SRT_MONTH, NO_BUSDAY_CONVENTION);

      err = swp_f_get_ref_rate_details(ref_rate, &basis, &compd);

      err = translate_compounding(&compdStr, compd);

      if (err) {
        free_option_computation(und_end_th, premium_array);
        return err;
      }
      err = translate_basis(&basisStr, basis);
      if (err) {
        free_option_computation(und_end_th, premium_array);
        return err;
      }

      /* Price forward reset caplet in the calibrated model */
      err = srt_f_fwd_resetcaplet(undptr, &grfnparam, yc_id, today, ref_rate,
                                  mdl_type, SRT_PAYER, vol_start[i], vol_end[i],
                                  und_start_act[i], und_end_th[i], compdStr,
                                  basisStr, &output_array[i]);

      if (err) {
        free_option_computation(und_end_th, premium_array);
        srt_f_destroy_und(und_id);
        return serror(
            "Error in srt_f_fwdvol_autocal (srt_f_fwd_resetcaplet): %s", err);
      }

      err = swp_f_LevelPayment(DTOL(und_start_act[i]), DTOL(und_end_th[i]),
                               compdStr, basisStr, yc_id, ref_rate, &lvl);
      if (err) {
        free_option_computation(und_end_th, premium_array);
        srt_f_destroy_und(und_id);
        return serror("Error in srt_f_fwdvol_autocal (srt_f_ForwardRate): %s",
                      err);
      }

      if (price_type == SRT_VOLATILITY) {

        premium = output_array[i];

        err = srt_f_optimpvol(premium, und_val[i], und_val[i],
                              (vol_end[i] - vol_start[i]) * YEARS_IN_DAY, lvl,
                              SRT_CALL, srt_fwd_vol_type, &(output_array[i]));
      }

      if (err) {
        free_option_computation(und_end_th, premium_array);
        srt_f_destroy_und(und_id);
        return serror("Error in srt_f_fwdvol_autocal (srt_f_optimpvol): %s",
                      err);
      }
    }

  } else if (option_type == RESETCMSOPTION) {

    und_end_th = (Date *)malloc(num_period * sizeof(Date));

    undptr = lookup_und(und_id);
    for (i = 0; i < num_period; i++) {
      err = translate_compounding(&compdStr, compd);
      if (err) {
        free_option_computation(und_end_th, premium_array);
        return err;
      }
      err = translate_basis(&basisStr, basis);
      if (err) {
        free_option_computation(und_end_th, premium_array);
        return err;
      }

      add_tenor(DTOL(und_start_act[i]), und_tenor, NO_BUSDAY_CONVENTION,
                &und_end_th[i]);

      err = srt_f_fwd_resetcms(
          undptr, &grfnparam, (long)vol_start[i], (long)vol_end[i],
          (long)und_start_act[i], (long)und_end_th[i], basisStr, compdStr,
          ref_rate, yc_id, GetVol, bs_vol_type, &output_array[i]);

      if (err) {
        free_option_computation(und_end_th, premium_array);
        srt_f_destroy_und(und_id);
        return serror("Error in srt_f_fwdvol_autocal (srt_f_fwd_resetcms): %s",
                      err);
      }

      SwpSpread = swp_f_spread(und_start_act[i], und_end_th[i], ref_rate);
      err = GetVol(und_start_act[i], und_end_th[i], und_val[i], und_val[i],
                   SwpSpread, &cmsvol);

      if (err) {
        free_option_computation(und_end_th, premium_array);
        srt_f_destroy_und(und_id);
        return serror("Error in srt_f_fwdvol_autocal (srt_f_ForwardRate): %s",
                      err);
      }

      if (price_type == SRT_VOLATILITY) {

        premium = output_array[i];
        err = srt_f_optimpvol(
            premium, und_val[i], und_val[i],
            (vol_end[i] - vol_start[i]) * YEARS_IN_DAY,
            swp_f_df(clcn_date, (long)und_start_act[i], yc_name), SRT_CALL,
            srt_fwd_vol_type, &output_array[i]);
      }

      if (err) {
        free_option_computation(und_end_th, premium_array);
        srt_f_destroy_und(und_id);
        return serror("Error in srt_f_fwdvol_autocal (srt_f_optimpvol): %s",
                      err);
      }
    }
  }

  srt_f_destroy_und(und_id);
  free_option_computation(und_end_th, premium_array);

  return NULL;
}

/* ---------------------------------------------------------------------------------------
 */

/* Function to compute forward volatility through parabolas for normal model on
 * one period */

/*******************************************************************************
 *
 * FUNCTION     	: srt_f_spfwdnormvol()
 *
 * PURPOSE     	: compute fwd normal vol of CMS(T2) between T1 and T2
 *				  using parabolic profiles
 *				  parabolic profiles are priced using (CMS)
 *caplets caplets are approximated by calls on CMS rate caplets are priced with
 *normal or lognormal smile
 *
 * AUTHOR		: Antoine Savine
 *
 * PARAMETERS  	: today 			- today (in days)
 *				: vol_start		    - start date (in
 *days) : vol_end			- end date (in days)
 *				: rate_tenor		- maturity of the desired
 *rate : fwd_cms_start
 *				  fwd_cms_end 		- forward rates at for start
 *and end dates CMS adjusted : level_slope_correl
 *				: slope_vol			- correl level/slope and
 *vol of the slope used for adj : GetVolTenor	    - function that enables to
 *retrieve implied vols : vol_id     		- ivol curve identification :
 *vol_type			- "normal" or "lognormal"
 *
 * RETURNS      	: *ans				- forward normal vol
 *
 *******************************************************************************/

/*	Core function for parabolic option */
Err srt_f_spfwdnormvol(
    Date today, Date vol_start, Date vol_end, SrtCompounding compd,
    double fwd_cms_start, double fwd_cms_end, double level_slope_correl,
    double slope_vol,
    Err (*GetVol)(Ddate, Ddate, double, double, double, double *),
    char *vol_type, double *fwd_vol, long max_strike, long max_vol,
    double delta_strike, double tol, long nvol) {
  double pc1, pc2, pp1, pp2;
  double time;
  double adjustment;
  Err err;

  if (vol_start > today) {
    /* Parabola atm call maturity 1 */
    err =
        parabola(today, vol_start, compd, fwd_cms_start, SRT_CALL, GetVol,
                 vol_type, &pc1, max_strike, max_vol, delta_strike, tol, nvol);
    if (err) {
      return serror(err);
    }

    /* Parabola atm put maturity 1 */
    err =
        parabola(today, vol_start, compd, fwd_cms_start, SRT_PUT, GetVol,
                 vol_type, &pp1, max_strike, max_vol, delta_strike, tol, nvol);
    if (err) {
      return serror(err);
    }
  } else {
    pc1 = pp1 = 0.0;
  }

  if (vol_end > today) {
    /* Parabola atm call maturity 2 */
    err =
        parabola(today, vol_end, compd, fwd_cms_end, SRT_CALL, GetVol, vol_type,
                 &pc2, max_strike, max_vol, delta_strike, tol, nvol);
    if (err) {
      return serror(err);
    }

    /* Parabola atm put maturity 2 */
    err =
        parabola(today, vol_end, compd, fwd_cms_end, SRT_PUT, GetVol, vol_type,
                 &pp2, max_strike, max_vol, delta_strike, tol, nvol);
    if (err) {
      return serror(err);
    }
  } else {
    pc2 = pp2 = 0.0;
  }

  *fwd_vol = pc2 + pp2 - pc1 - pp1;

  adjustment = level_slope_correl * slope_vol *
               sqrt(YEARS_IN_DAY * (vol_start - today) * 4 * (pc1 + pp1));

  *fwd_vol -= adjustment;

  if ((*fwd_vol) < 0) {
    return serror("negative variance");
  }

  time = YEARS_IN_DAY * (vol_end - vol_start);
  if (time <= 0) {
    return serror("invalid maturities");
  }

  *fwd_vol = sqrt((*fwd_vol) / time);

  return NULL;
}

/* ---------------------------------------------------------------------------------
 */

/*	Main function to compute the forwad volatility for a string of period in
 * a normal model  */
Err srt_f_fwdvol_parabola(double level_slope_correl, /* Parameters */
                          double slope_vol,
                          long num_period, /* Number of required forward vols */
                          Ddate *vol_start,   /* Vol Start Dates array */
                          Ddate *vol_end,     /* Vol End Dates array */
                          String ref_rate,    /* Reference Rate Code */
                          String yc_id,       /* YC id */
                          String bs_vol_type, /* Lognormal or normal */
                          Err(*GetVol)        /* GetVol function */
                          (Ddate start, Ddate end, double strike,
                           double dForward, double dSpread, double *bs_vol),
                          double *fwd_vol, /* Output array */
                          long max_strike, long max_vol, double delta_strike,
                          double tol, long nvol) {
  Err err;
  SrtCurvePtr yldcrv;
  Date today;
  long spot_lag;
  SrtCompounding compd;
  BasisCode basis;
  SrtDiffusionType srt_bs_vol_type;
  SwapDP swapdp_st, swapdp_nd;
  double start_rate_st, end_rate_st, atm_vol_st, cms_st, start_rate_nd,
      end_rate_nd, atm_vol_nd, cms_nd;
  long i;

  /* Get Today and Spot Lag */
  yldcrv = lookup_curve(yc_id);
  if (!yldcrv) {
    return serror("Error in srt_f_fwdvol_parabola: YC not found");
  }
  today = get_clcndate_from_yldcrv(yldcrv);
  spot_lag = get_spotlag_from_curve(yldcrv);

  /* Get the details of the reference rate: compounding (in value and string)
   * and basis */
  err = swp_f_get_ref_rate_details(ref_rate, &basis, &compd);
  if (err) {
    return err;
  }

  /* Get the diffusion type */
  if (err = interp_diffusion_type(bs_vol_type, &srt_bs_vol_type)) {
    return serror("Error in srt_f_fwdvol_parabola (interp_vol_type): %s", err);
  }

  /*VISION Y2K Fix Number 1
  Check dates in volatility arrays

  err = srt_f_test_date_dblarr ((double *) vol_start  , (int) num_period);
  if (err)
  {
          return serror ("Volatility start dates should be between 01-Jan-1983
  and 01-Jan-2070");
  }

  err = srt_f_test_date_dblarr ((double *) vol_end  , (int) num_period);
  if (err)
  {
          return serror ("Volatility end dates should be between 01-Jan-1983 and
  01-Jan-2070");
  }

  End of Y2K Fix Block	*/

  /* Compute the fwd volatility for each period */
  for (i = 0; i < num_period; i++) {
    /* Get the CMS rate adjusted at vol start date */
    start_rate_st =
        add_unit((Date)vol_start[i], spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
    end_rate_st = add_unit((Date)start_rate_st, 12 / (int)compd, SRT_MONTH,
                           NO_BUSDAY_CONVENTION);

    err = GetVol(start_rate_st, end_rate_st, 0.0, 1.0, 0.0, &atm_vol_st);
    if (err)
      return err;

    err = swp_f_setSwapDP((Date)start_rate_st, 1, compd, basis, &swapdp_st);
    if (err)
      return err;

    /* Set the spotlag in the SwapDP */
    swapdp_st.spot_lag = spot_lag;

    err = swp_f_cmswp(&swapdp_st, ref_rate, 0.0, atm_vol_st, SRT_PAYER,
                      (Date)end_rate_st, DEFAULT_CMS_DELTA, MAX_CMS_SWAPS,
                      yc_id, srt_bs_vol_type, &cms_st);
    if (err)
      return err;

    /* Get the CMS rate adjusted at vol end date */
    start_rate_nd =
        add_unit((Date)vol_end[i], spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
    end_rate_nd = add_unit((Date)start_rate_nd, 12 / (int)compd, SRT_MONTH,
                           NO_BUSDAY_CONVENTION);

    err = GetVol(start_rate_nd, end_rate_nd, 0.0, 1.0, 0.0, &atm_vol_nd);
    if (err)
      return err;

    err = swp_f_setSwapDP((Date)start_rate_nd, 1, compd, basis, &swapdp_nd);
    if (err)
      return err;

    /* Set the spotlag in the SwapDP */
    swapdp_nd.spot_lag = spot_lag;

    err = swp_f_cmswp(&swapdp_nd, ref_rate, 0.0, atm_vol_nd, SRT_PAYER,
                      (Date)end_rate_nd, DEFAULT_CMS_DELTA, MAX_CMS_SWAPS,
                      yc_id, srt_bs_vol_type, &cms_nd);
    if (err)
      return err;

    /* Compute forward normal vol for this period */
    err = srt_f_spfwdnormvol(today, (Date)(vol_start[i]), (Date)(vol_end[i]),
                             compd, cms_st, cms_nd, level_slope_correl,
                             slope_vol, GetVol, bs_vol_type, &(fwd_vol[i]),
                             max_strike, max_vol, delta_strike, tol, nvol);
    if (err)
      return err;
  }

  return NULL;
}