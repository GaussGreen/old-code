/*	New SABR pricing functions	*/
/*	Dec99	*/

#include "BGMUtils.h"
#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/*	Interpret volatility type	*/
Err interp_sabr_vol_type(char *str, SABR_VOL_TYPE *val) {
  strupper(str);
  strip_white_space(str);

  if (!strcmp(str, "LOGNORMAL") || !strcmp(str, "LOG") || !strcmp(str, "L") ||
      !strcmp(str, "0") || !strcmp(str, "ATMLOG")) {
    *val = SABR_ATM_LOG;
    return NULL;
  }

  if (!strcmp(str, "NORMAL") || !strcmp(str, "NORM") || !strcmp(str, "N") ||
      !strcmp(str, "1") || !strcmp(str, "ATMNORM")) {
    *val = SABR_ATM_NORM;
    return NULL;
  }

  if (!strcmp(str, "BETAVOL") || !strcmp(str, "BETA") || !strcmp(str, "BV") ||
      !strcmp(str, "B") || !strcmp(str, "2")) {
    *val = SABR_ATM_BETA;
    return NULL;
  }

  if (!strcmp(str, "LOGNORMALATSTRIKE") || !strcmp(str, "BLACKSCHOLES") ||
      !strcmp(str, "BS") || !strcmp(str, "STRLOG")) {
    *val = SABR_STR_LOG;
    return NULL;
  }

  if (!strcmp(str, "NORMALATSTRIKE") || !strcmp(str, "BLACKSCHOLESNORMAL") ||
      !strcmp(str, "BSN") || !strcmp(str, "STRNORM")) {
    *val = SABR_STR_NORM;
    return NULL;
  }

  return serror("Vol type: enter L (LOGNORMAL)  , N (NORMAL) or B (Beta)");
}

/*	Setup parameters	*/
Err op_sabr_set_param(int num_param, char **param_str, char **value_str,
                      SABR_RISK_PARAM *param) {
  int i;
  double d;
  SABR_VOL_TYPE vt;
  Err err;

  /*	Set defaults	*/

  param->delta_shift = srt_f_getdeltashift();
  param->gamma_shift = srt_f_getdeltashift();
  param->vanna_shift = srt_f_getvannashift();
  param->vega_shift = srt_f_getvegashift();
  param->volga_shift = srt_f_getvegashift();
  param->vanna_shift = srt_f_getvannashift();
  param->theta_shift = srt_f_getthetashift();
  param->alpha_shift = srt_f_getalphashift();
  param->beta_shift = srt_f_getbetashift();
  param->rho_shift = srt_f_getrhoshift();

  param->vega_mult = 1;
  param->volga_mult = 0;

  param->delta_fac = 1.0; /* 0.0001; */
  param->gamma_fac = 1.0; /* 0.0001; */
  param->vega_fac = 10;
  param->volga_fac = 10;
  param->vanna_fac = 10;
  param->theta_fac = 1.0; /*YEARS_IN_DAY;*/
  param->alpha_fac = 1.0; /* 0.01; */
  param->beta_fac = 1.0;  /* 0.01; */
  param->rho_fac = 1.0;   /*0.01; */

  param->delta_freeze_vol = SABR_ATM_BETA;
  param->theta_freeze_vol = SABR_ATM_BETA;
  param->alpha_beta_rho_freeze_vol = SABR_ATM_LOG;
  param->vega_bump_vol = SABR_ATM_LOG;

  /*	Read parameters	*/
  for (i = 0; i < num_param; i++) {
    strupper(param_str[i]);
    strip_white_space(param_str[i]);
    strupper(value_str[i]);
    strip_white_space(value_str[i]);

    if (strlen(param_str[i]) == 0 || strlen(value_str[i]) == 0)
      continue;

    if (!strcmp(param_str[i], "DELTASHIFT")) {
      sscanf(value_str[i], "%lf", &d);
      param->delta_shift = d;
    } else if (!strcmp(param_str[i], "GAMMASHIFT")) {
      sscanf(value_str[i], "%lf", &d);
      param->gamma_shift = d;
    } else if (!strcmp(param_str[i], "VEGASHIFT")) {
      sscanf(value_str[i], "%lf", &d);
      param->vega_shift = d;
    } else if (!strcmp(param_str[i], "VOLGASHIFT")) {
      sscanf(value_str[i], "%lf", &d);
      param->volga_shift = d;
    } else if (!strcmp(param_str[i], "VANNASHIFT")) {
      sscanf(value_str[i], "%lf", &d);
      param->vanna_shift = d;
    } else if (!strcmp(param_str[i], "THETASHIFT")) {
      sscanf(value_str[i], "%lf", &d);
      param->theta_shift = d;
    } else if (!strcmp(param_str[i], "ALPHASHIFT")) {
      sscanf(value_str[i], "%lf", &d);
      param->alpha_shift = d;
    } else if (!strcmp(param_str[i], "BETASHIFT")) {
      sscanf(value_str[i], "%lf", &d);
      param->beta_shift = d;
    } else if (!strcmp(param_str[i], "RHOSHIFT")) {
      sscanf(value_str[i], "%lf", &d);
      param->rho_shift = d;
    } else if (!strcmp(param_str[i], "DELTAFAC")) {
      sscanf(value_str[i], "%lf", &d);
      param->delta_fac = d;
    } else if (!strcmp(param_str[i], "GAMMAFAC")) {
      sscanf(value_str[i], "%lf", &d);
      param->gamma_fac = d;
    } else if (!strcmp(param_str[i], "VEGAFAC")) {
      sscanf(value_str[i], "%lf", &d);
      param->vega_fac = d;
    } else if (!strcmp(param_str[i], "VANNAFAC")) {
      sscanf(value_str[i], "%lf", &d);
      param->vanna_fac = d;
    } else if (!strcmp(param_str[i], "VOLGAFAC")) {
      sscanf(value_str[i], "%lf", &d);
      param->volga_fac = d;
    } else if (!strcmp(param_str[i], "THETAFAC")) {
      sscanf(value_str[i], "%lf", &d);
      param->theta_fac = d;
    } else if (!strcmp(param_str[i], "ALPHAFAC")) {
      sscanf(value_str[i], "%lf", &d);
      param->alpha_fac = d;
    } else if (!strcmp(param_str[i], "BETAFAC")) {
      sscanf(value_str[i], "%lf", &d);
      param->beta_fac = d;
    } else if (!strcmp(param_str[i], "RHOFAC")) {
      sscanf(value_str[i], "%lf", &d);
      param->rho_fac = d;
    } else if (!strcmp(param_str[i], "VEGAMULT")) {
      if (!strcmp(value_str[i], "YES")) {
        param->vega_mult = 1;
      } else if (!strcmp(value_str[i], "NO")) {
        param->vega_mult = 0;
      } else {
        return serror("VEGAMULT must be YES or NO");
      }
    } else if (!strcmp(param_str[i], "VOLGAMULT")) {
      if (!strcmp(value_str[i], "YES")) {
        param->volga_mult = 1;
      } else if (!strcmp(value_str[i], "NO")) {
        param->volga_mult = 0;
      } else {
        return serror("VOLGAMULT must be YES or NO");
      }
    } else if (!strcmp(param_str[i], "DELTAFREEZEVOL")) {
      if (err = interp_sabr_vol_type(value_str[i], &vt)) {
        return err;
      } else {
        param->delta_freeze_vol = vt;
      }
    } else if (!strcmp(param_str[i], "THETAFREEZEVOL")) {
      if (err = interp_sabr_vol_type(value_str[i], &vt)) {
        return err;
      } else {
        param->theta_freeze_vol = vt;
      }
    } else if (!strcmp(param_str[i], "SMILEFREEZEVOL")) {
      if (err = interp_sabr_vol_type(value_str[i], &vt)) {
        return err;
      } else {
        param->alpha_beta_rho_freeze_vol = vt;
      }
    } else if (!strcmp(param_str[i], "VEGABUMPVOL")) {
      if (err = interp_sabr_vol_type(value_str[i], &vt)) {
        return err;
      } else {
        param->vega_bump_vol = vt;
      }
    } else {
      return serror("Unknown SABR parameter: %s", param_str[i]);
    }
  }

  return NULL;
}

/*	Conversion from vol type to diffusion type and ATM	*/
static Err vol_type_to_diff_type_atm(SABR_VOL_TYPE vol_type,
                                     SrtDiffusionType *diff_type,
                                     /*	0: Strike  , 1: ATM	*/
                                     int *atm) {
  switch (vol_type) {
  case SABR_ATM_BETA:
    *diff_type = SRT_BETAVOL;
    *atm = 1;
    break;
  case SABR_ATM_LOG:
    *diff_type = SRT_LOGNORMAL;
    *atm = 1;
    break;
  case SABR_ATM_NORM:
    *diff_type = SRT_NORMAL;
    *atm = 1;
    break;
  case SABR_STR_LOG:
    *diff_type = SRT_LOGNORMAL;
    *atm = 0;
    break;
  case SABR_STR_NORM:
    *diff_type = SRT_NORMAL;
    *atm = 0;
    break;
  default:
    return serror("Bad SABR vol type %d", vol_type);
    break;
  }

  return NULL;
}

/*	Conversion from one vol type to another	*/
Err vol_conv(double input_vol, SABR_VOL_TYPE input_vol_type, double *output_vol,
             SABR_VOL_TYPE output_vol_type,
             /*	Forward	*/
             double forward,
             /*	Option specs	*/
             double strike, double maturity,
             /*	Model specs	*/
             double alpha, double beta, double rho) {
  SrtDiffusionType input_diff_type, output_diff_type;
  int input_vol_atm, output_vol_atm;
  double tmp_vol;
  Err err;

  /*	Get input diffusion type and strike	*/
  if (err = vol_type_to_diff_type_atm(input_vol_type, &input_diff_type,
                                      &input_vol_atm)) {
    return err;
  }

  /*	Get output diffusion type and strike	*/
  if (err = vol_type_to_diff_type_atm(output_vol_type, &output_diff_type,
                                      &output_vol_atm)) {
    return err;
  }

  /*	Case 1: vol input is ATM  , we can use Jerome;s function
          (srt_f_optsarbvol)	*/
  if (input_vol_atm) {
    /*	If output vol is also ATM  , strike information is irrelevant	*/
    if (output_vol_atm) {
      strike = forward;
    }

    /*	Call to srt_f_optsarbvol	*/
    if (err = srt_f_optsarbvol(forward, strike, maturity, input_vol, alpha,
                               beta, rho, input_diff_type, output_diff_type,
                               output_vol)) {
      return err;
    }

    return NULL;
  }

  /*	Case 2: vol input is not ATM  , we cannot use Jerome;s function
          (srt_f_optsarbvol)	*/

  /*	Input is lognormal	*/
  if (input_diff_type == SRT_LOGNORMAL) {
    /*	Transform into beta vol	*/
    if (err = srt_f_optblkvolATMtobetavolStochVol(
            forward, strike, input_vol, maturity, alpha, beta, rho, &tmp_vol)) {
      return err;
    }

    /*	Transform beta vol into what is required	*/
    return vol_conv(tmp_vol, SABR_ATM_BETA, output_vol, output_vol_type,
                    forward, strike, maturity, alpha, beta, rho);
  }

  /*	Input is normal: transform into lognormal	*/
  if (err = srt_f_optbetavoltoblkvol(forward, strike, input_vol, maturity, 0.0,
                                     &tmp_vol)) {
    return err;
  }

  /*	Call function with lognormal input	*/
  return vol_conv(tmp_vol, SABR_STR_LOG, output_vol, output_vol_type, forward,
                  strike, maturity, alpha, beta, rho);
}

/*	Pricing	*/
/*	All parameters checks are supposed to have been done before call
 */
Err op_sabr_pricing(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha, double beta, double rho,
    /*	Vol input	*/
    double input_vol, SABR_VOL_TYPE input_vol_type,
    /*	Premium	*/
    double *prem) {
  double str_log_vol;
  Err err;

  /*	Transforms vol into the equivalent BS vol for the considered strike
   */
  if (err = vol_conv(input_vol, input_vol_type, &str_log_vol, SABR_STR_LOG,
                     forward, strike, maturity, alpha, beta, rho)) {
    return err;
  }

  /*	Compute premium using Black-Scholes	*/

  *prem = srt_f_optblksch(forward, strike, str_log_vol, maturity, disc,
                          call_put, PREMIUM);

  return NULL;
}

/*	Delta/Gamma	*/
/*	All parameters checks are supposed to have been done before call
 */
Err op_sabr_delta_gamma(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha, double beta, double rho,
    /*	Vol input	*/
    double input_vol, SABR_VOL_TYPE input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double delta_shift, double gamma_shift,
    /*	Which vol to freeze when bumping	*/
    SABR_VOL_TYPE freeze_vol_type,
    /*	Calculate delta or gamma	*/
    /*	0: Delta  , 1: Gamma	*/
    int delta_gamma,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answer	*/
    double *sens) {
  double frozen_vol;
  double prem1, prem2, delta1, delta2;
  Err err;

  /*	Calculate equivalent vol of the type to be frozen	*/
  if (err = vol_conv(input_vol, input_vol_type, &frozen_vol, freeze_vol_type,
                     forward, strike, maturity, alpha, beta, rho)) {
    return err;
  }

  /*	Calculate premium before shift	*/
  if (err = op_sabr_pricing(forward, strike, maturity, disc, call_put, alpha,
                            beta, rho, frozen_vol, freeze_vol_type, &prem1)) {
    return err;
  }

  /*	Calculate premium after shift	*/
  if (err = op_sabr_pricing(forward + delta_shift, strike, maturity, disc,
                            call_put, alpha, beta, rho, frozen_vol,
                            freeze_vol_type, &prem2)) {
    return err;
  }

  /*	Calculate delta before shift	*/
  delta1 = (prem2 - prem1) / delta_shift;
  *sens = fac * delta1;

  /*	If gamma is to be calculated  , calculate delta after shift	*/
  if (delta_gamma == 1) {
    if (err = op_sabr_delta_gamma(forward + gamma_shift, strike, maturity, disc,
                                  call_put, alpha, beta, rho, frozen_vol,
                                  freeze_vol_type, delta_shift, gamma_shift,
                                  freeze_vol_type, 0, 1.0, &delta2)) {
      return err;
    }

    /*	Calculate gamma	*/
    *sens = fac * (delta2 - delta1) / gamma_shift;
  }

  return NULL;
}

/*	Vega/Volga	*/
/*	All parameters checks are supposed to have been done before call
 */
Err op_sabr_vega_volga(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha, double beta, double rho,
    /*	Vol input	*/
    double input_vol, SABR_VOL_TYPE input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double vega_shift, double volga_shift,
    /*	Multiplicative or additive	*/
    /*	0: Additive  , 1: Multiplicative	*/
    int vega_mult, int volga_mult,
    /*	Which vol to bump	*/
    SABR_VOL_TYPE bump_vol_type,
    /*	Calculate vega or volga	*/
    /*	0: Vega  , 1: Volga	*/
    int vega_volga,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answer	*/
    double *sens) {
  double bump_vol;
  double prem1, prem2, vega1, vega2;
  double act_vega_shift, act_volga_shift;
  Err err;

  /*	Calculate equivalent vol of the type to be bumped	*/
  if (err = vol_conv(input_vol, input_vol_type, &bump_vol, bump_vol_type,
                     forward, strike, maturity, alpha, beta, rho)) {
    return err;
  }

  /*	Calculate premium before shift	*/
  if (err = op_sabr_pricing(forward, strike, maturity, disc, call_put, alpha,
                            beta, rho, bump_vol, bump_vol_type, &prem1)) {
    return err;
  }

  /*	Calculate actual shifts	*/

  if (vega_mult == 1) {
    act_vega_shift = vega_shift * bump_vol;
  } else {
    act_vega_shift = vega_shift;
  }

  if (volga_mult == 1) {
    act_volga_shift = volga_shift * bump_vol;
  } else {
    act_volga_shift = volga_shift;
  }

  /*	Calculate premium after shift	*/
  if (err = op_sabr_pricing(forward, strike, maturity, disc, call_put, alpha,
                            beta, rho, bump_vol + act_vega_shift, bump_vol_type,
                            &prem2)) {
    return err;
  }

  /*	Calculate vega before shift	*/
  vega1 = (prem2 - prem1) / vega_shift;
  *sens = fac * vega1;

  /*	If volga is to be calculated  , calculate delta after shift	*/
  if (vega_volga == 1) {
    if (err = op_sabr_vega_volga(
            forward, strike, maturity, disc, call_put, alpha, beta, rho,
            bump_vol + act_volga_shift, bump_vol_type, vega_shift, volga_shift,
            vega_mult, volga_mult, bump_vol_type, 0, 1.0, &vega2)) {
      return err;
    }

    /*	Calculate volga	*/
    *sens = fac * (vega2 - vega1) / volga_shift;
  }

  return NULL;
}

/*	Vanna	*/
/*	All parameters checks are supposed to have been done before call
 */
Err op_sabr_vanna(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha, double beta, double rho,
    /*	Vol input	*/
    double input_vol, SABR_VOL_TYPE input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double delta_shift, double vega_shift,
    /*	Multiplicative or additive vega shift */
    /*	0: Additive  , 1: Multiplicative	*/
    int vega_mult,
    /*	Which vol to bump	*/
    SABR_VOL_TYPE bump_vol_type,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answer	*/
    double *vanna) {
  double bump_vol;
  double delta1, delta2;
  double act_vega_shift;
  Err err;

  /*	Calculate equivalent vol of the type to be bumped	*/
  if (err = vol_conv(input_vol, input_vol_type, &bump_vol, bump_vol_type,
                     forward, strike, maturity, alpha, beta, rho)) {
    return err;
  }

  /*	Calculate delta before shift	*/
  if (err = op_sabr_delta_gamma(
          forward, strike, maturity, disc, call_put, alpha, beta, rho, bump_vol,
          bump_vol_type, delta_shift, 0.0, bump_vol_type, 0, 1.0, &delta1)) {
    return err;
  }

  /*	Calculate actual shifts	*/

  if (vega_mult == 1) {
    act_vega_shift = vega_shift * bump_vol;
  } else {
    act_vega_shift = vega_shift;
  }

  /*	Calculate delta after shift	*/
  if (err = op_sabr_delta_gamma(forward, strike, maturity, disc, call_put,
                                alpha, beta, rho, bump_vol + act_vega_shift,
                                bump_vol_type, delta_shift, 0.0, bump_vol_type,
                                0, 1.0, &delta2)) {
    return err;
  }

  /*	Calculate vanna	*/
  *vanna = fac * (delta2 - delta1) / vega_shift;

  return NULL;
}

/*	Theta	*/
/*	All parameters checks are supposed to have been done before call
 */
Err op_sabr_theta(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha, double beta, double rho,
    /*	Vol input	*/
    double input_vol, SABR_VOL_TYPE input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double time_shift,
    /*	Which vol to freeze when bumping	*/
    SABR_VOL_TYPE freeze_vol_type,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answers	*/
    double *theta, double *thetagamma, double *thetavolga, double *thetavanna,
    double *thetacarry) {
  double frozen_vol;
  double prem1, prem2;
  double r, sigma_beta, gamma, volga, vanna;
  Err err;

  /*	Calculate equivalent vol of the type to be frozen	*/
  if (err = vol_conv(input_vol, input_vol_type, &frozen_vol, freeze_vol_type,
                     forward, strike, maturity, alpha, beta, rho)) {
    return err;
  }

  /*	Calculate premium before shift	*/
  if (err = op_sabr_pricing(forward, strike, maturity, disc, call_put, alpha,
                            beta, rho, frozen_vol, freeze_vol_type, &prem1)) {
    return err;
  }

  /*	Calculate premium after shift	*/
  if (err = op_sabr_pricing(forward, strike, maturity - time_shift, disc,
                            call_put, alpha, beta, rho, frozen_vol,
                            freeze_vol_type, &prem2)) {
    return err;
  }

  /*	Calculate theta	*/
  *theta = (prem2 - prem1) / time_shift;

  /*	Calculate carry	*/
  r = -log(disc) / maturity;
  *thetacarry = r * prem1;
  *theta += *thetacarry;

  /*	Calculate sigma beta  , gamma  , volga and vanna	*/

  if (err = vol_conv(input_vol, input_vol_type, &sigma_beta, SABR_ATM_BETA,
                     forward, strike, maturity, alpha, beta, rho)) {
    return err;
  }

  if (err =
          op_sabr_delta_gamma(forward, strike, maturity, disc, call_put, alpha,
                              beta, rho, sigma_beta, SABR_ATM_BETA, 0.0001,
                              0.0001, SABR_ATM_BETA, 1, 1.0, &gamma)) {
    return err;
  }

  if (err = op_sabr_vega_volga(forward, strike, maturity, disc, call_put, alpha,
                               beta, rho, sigma_beta, SABR_ATM_BETA, 0.0001,
                               0.0001, 0, 0, SABR_ATM_BETA, 1, 1.0, &volga)) {
    return err;
  }

  if (err = op_sabr_vanna(forward, strike, maturity, disc, call_put, alpha,
                          beta, rho, sigma_beta, SABR_ATM_BETA, 0.0001, 0.0001,
                          0, SABR_ATM_BETA, 1.0, &vanna)) {
    return err;
  }

  *thetagamma = -0.5 * pow(sigma_beta, 2.0) * pow(forward, 2.0 * beta) * gamma;

  *thetavolga = -0.5 * pow(alpha, 2.0) * pow(sigma_beta, 2.0) * volga;

  *thetavanna =
      -rho * pow(sigma_beta, 2.0) * pow(forward, beta) * alpha * vanna;

  *theta *= fac;
  *thetagamma *= fac;
  *thetavolga *= fac;
  *thetavanna *= fac;
  *thetacarry *= fac;

  return NULL;
}

/*	Sensitivity to model parameters	*/
/*	All parameters checks are supposed to have been done before call
 */
Err op_sabr_model_sens(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha, double beta, double rho,
    /*	Vol input	*/
    double input_vol, SABR_VOL_TYPE input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double alpha_shift, double beta_shift, double rho_shift,
    /*	Which vol to freeze when bumping	*/
    SABR_VOL_TYPE freeze_vol_type,
    /*	Factor that multiplies results	*/
    double alpha_fac, double beta_fac, double rho_fac,
    /*	Answers	*/
    double *alpha_sens, double *beta_sens, double *rho_sens) {
  double frozen_vol;
  double prem1, prem2, prem3, prem4;
  Err err;

  /*	Calculate equivalent vol of the type to be frozen	*/
  if (err = vol_conv(input_vol, input_vol_type, &frozen_vol, freeze_vol_type,
                     forward, strike, maturity, alpha, beta, rho)) {
    return err;
  }

  /*	Calculate premium before shift	*/
  if (err = op_sabr_pricing(forward, strike, maturity, disc, call_put, alpha,
                            beta, rho, frozen_vol, freeze_vol_type, &prem1)) {
    return err;
  }

  /*	Calculate premium after shift of alpha	*/
  if (err = op_sabr_pricing(forward, strike, maturity, disc, call_put,
                            alpha + alpha_shift, beta, rho, frozen_vol,
                            freeze_vol_type, &prem2)) {
    return err;
  }

  /*	Calculate premium after shift of beta	*/
  if (err = op_sabr_pricing(forward, strike, maturity, disc, call_put, alpha,
                            beta + beta_shift, rho, frozen_vol, freeze_vol_type,
                            &prem3)) {
    return err;
  }

  /*	Calculate premium after shift of rho	*/
  if (err = op_sabr_pricing(forward, strike, maturity, disc, call_put, alpha,
                            beta, rho + rho_shift, frozen_vol, freeze_vol_type,
                            &prem4)) {
    return err;
  }

  /*	Calculate sensitivites	*/
  *alpha_sens = alpha_fac * (prem2 - prem1) / alpha_shift;
  *beta_sens = beta_fac * (prem3 - prem1) / beta_shift;
  *rho_sens = rho_fac * (prem4 - prem1) / rho_shift;

  return NULL;
}

#undef R
#undef C
#undef SHIFT2
#undef SHIFT3
