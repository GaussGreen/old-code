/*	New SABRAF pricing functions	*/
/*	Aug04	*/
/*	J.B / P.T*/

#include "BGMUtils.h"
#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"
#include <opsabraf.h"
#include <opsabrgenericinterp.h"

/*	Interpret volatility type	*/
Err interp_sabraf_vol_type(char *str, SrtDiffusionType *val) {
  Err err = NULL;
  if (strcmp(str, "") != 0)
    err = interp_diffusion_type(str, val);
  return err;
}

/*	Setup parameters	*/
Err op_sabraf_set_param(int num_param, char **param_str, char **value_str,
                        SABRAF_RISK_PARAM *param) {
  int i;
  double d;
  SrtDiffusionType vt;
  Err err;

  /*	Set defaults	*/

  param->delta_shift = srt_f_getdeltashift();
  param->gamma_shift = srt_f_getdeltashift();
  param->vega_shift = srt_f_getvegashift();
  param->volga_shift = srt_f_getvegashift();
  param->theta_shift = srt_f_getthetashift();
  param->alpha_shift = srt_f_getalphashift();
  param->beta_shift = srt_f_getbetashift();
  param->rho_shift = srt_f_getrhoshift();
  param->zeta_shift = srt_f_getzetashift();

  param->vega_mult = 1;
  param->volga_mult = 0;
  param->alpha_mult = 0;

  param->alpha_mult = 0;
  param->beta_mult = 0;
  param->rho_mult = 0;
  param->zeta_mult = 0;

  param->delta_fac = 1.0; /* 0.0001; */
  param->gamma_fac = 1.0; /* 0.0001; */
  param->vega_fac = 10;
  param->volga_fac = 10;
  param->vanna_fac = 10;
  param->theta_fac = 1.0; /*YEARS_IN_DAY;*/
  param->alpha_fac = 1.0; /* 0.01; */
  param->beta_fac = 1.0;  /* 0.01; */
  param->rho_fac = 1.0;   /*0.01; */
  param->zeta_fac = 1.0;

  param->delta_freeze_vol = SABR_ATM_BETA;
  param->theta_freeze_vol = SABR_ATM_BETA;
  param->alpha_beta_rho_zeta_freeze_vol = SABR_ATM_LOG;
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
    } else if (!strcmp(param_str[i], "ZETASHIFT")) {
      sscanf(value_str[i], "%lf", &d);
      param->zeta_shift = d;
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
    } else if (!strcmp(param_str[i], "ZETAFAC")) {
      sscanf(value_str[i], "%lf", &d);
      param->zeta_fac = d;
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
    } else if (!strcmp(param_str[i], "ALPHAMULT")) {
      if (!strcmp(value_str[i], "YES")) {
        param->alpha_mult = 1;
      } else if (!strcmp(value_str[i], "NO")) {
        param->alpha_mult = 0;
      } else {
        return serror("ALPHAMULT must be YES or NO");
      }
    } else if (!strcmp(param_str[i], "BETAMULT")) {
      if (!strcmp(value_str[i], "YES")) {
        param->alpha_mult = 1;
      } else if (!strcmp(value_str[i], "NO")) {
        param->alpha_mult = 0;
      } else {
        return serror("BETAMULT must be YES or NO");
      }
    } else if (!strcmp(param_str[i], "RHOMULT")) {
      if (!strcmp(value_str[i], "YES")) {
        param->alpha_mult = 1;
      } else if (!strcmp(value_str[i], "NO")) {
        param->alpha_mult = 0;
      } else {
        return serror("RHOMULT must be YES or NO");
      }
    } else if (!strcmp(param_str[i], "ZETAMULT")) {
      if (!strcmp(value_str[i], "YES")) {
        param->alpha_mult = 1;
      } else if (!strcmp(value_str[i], "NO")) {
        param->alpha_mult = 0;
      } else {
        return serror("ZETAMULT must be YES or NO");
      }
    } else if (!strcmp(param_str[i], "DELTAFREEZEVOL")) {
      if (err = interp_sabraf_vol_type(value_str[i], &vt)) {
        return err;
      } else {
        param->delta_freeze_vol = vt;
      }
    } else if (!strcmp(param_str[i], "THETAFREEZEVOL")) {
      if (err = interp_sabraf_vol_type(value_str[i], &vt)) {
        return err;
      } else {
        param->theta_freeze_vol = vt;
      }
    } else if (!strcmp(param_str[i], "SMILEFREEZEVOL")) {
      if (err = interp_sabraf_vol_type(value_str[i], &vt)) {
        return err;
      } else {
        param->alpha_beta_rho_zeta_freeze_vol = vt;
      }
    } else if (!strcmp(param_str[i], "VEGABUMPVOL")) {
      if (err = interp_sabraf_vol_type(value_str[i], &vt)) {
        return err;
      } else {
        param->vega_bump_vol = vt;
      }
    } else {
      return serror("Unknown SABRAF parameter: %s", param_str[i]);
    }
  }

  return NULL;
}

/*	Pricing	*/
/*	All parameters checks are supposed to have been done before call
 */

/////////////////////////////////////////////////////////////////////////////////////////////
//
//	srt_f_op_sabraf_price()  ,
//
//	Returns Price(K) Given one ATM Vol and the smiles parameters
//	input SigmaBeta ATMLog  , ATmNorm
//	output Log or Norm
//
/////////////////////////////////////////////////////////////////////////////////////////////

Err srt_f_op_sabraf_price(double Forward, double Strike, double Maturity,
                          double Disc, SrtCallPutType call_put, double Sigma,
                          double Alpha, double Beta, double Rho, double Zeta,
                          SrtDiffusionType input, double *price) {
  double SigmaBeta;
  Err err = NULL;

  // Process According to imput
  err = srt_f_optbmm2vol(Forward, Strike, Maturity, Sigma, Alpha, Beta, Rho,
                         Zeta, input, SRT_BETAVOL, &SigmaBeta);
  if (err)
    goto FREE_RETURN;

  err = BMM_Option_Price2(Forward, Strike, Maturity, SigmaBeta, Alpha, Beta,
                          Rho, Zeta, call_put, price);
  if (err)
    goto FREE_RETURN;
  *price *= Disc;

FREE_RETURN:
  return err;
}

/*	Delta/Gamma	*/
/*	All parameters checks are supposed to have been done before call
 */
Err op_sabraf_delta_gamma(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha, double beta, double rho, double zeta,
    /*	Vol input	*/
    double input_vol, SrtDiffusionType input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double delta_shift, double gamma_shift,
    /*	Which vol to freeze when bumping	*/
    SrtDiffusionType freeze_vol_type,
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

  /*	Calculate equivalent vol of the type to be frozen	only ATM types or
   * BETA*/
  err =
      srt_f_optbmm2vol(forward, forward, maturity, input_vol, alpha, beta, rho,
                       zeta, input_vol_type, freeze_vol_type, &frozen_vol);
  if (err)
    return err;

  /*	Calculate premium before shift	*/
  if (err = srt_f_op_sabraf_price(forward, strike, maturity, disc, call_put,
                                  frozen_vol, alpha, beta, rho, zeta,
                                  freeze_vol_type, &prem1)) {
    return err;
  }

  /*	Calculate premium after shift	*/
  if (err = srt_f_op_sabraf_price(forward + delta_shift, strike, maturity, disc,
                                  call_put, frozen_vol, alpha, beta, rho, zeta,
                                  freeze_vol_type, &prem2)) {
    return err;
  }

  /*	Calculate delta before shift	*/
  delta1 = (prem2 - prem1) / delta_shift;
  *sens = fac * delta1;

  /*	If gamma is to be calculated  , calculate delta after shift	*/
  if (delta_gamma == 1) {
    if (err = op_sabraf_delta_gamma(
            forward + gamma_shift, strike, maturity, disc, call_put, alpha,
            beta, rho, zeta, frozen_vol, freeze_vol_type, delta_shift,
            gamma_shift, freeze_vol_type, 0, 1.0, &delta2)) {
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
Err op_sabraf_vega_volga(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha, double beta, double rho, double zeta,
    /*	Vol input	*/
    double input_vol, SrtDiffusionType input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double vega_shift, double volga_shift,
    /*	Multiplicative or additive	*/
    /*	0: Additive  , 1: Multiplicative	*/
    int vega_mult, int volga_mult,
    /*	Which vol to bump	*/
    SrtDiffusionType bump_vol_type,
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
  err = srt_f_optbmm2vol(forward, forward, maturity, input_vol, alpha, beta,
                         rho, zeta, input_vol_type, bump_vol_type, &bump_vol);
  if (err)
    return err;

  /*	Calculate premium before shift	*/
  if (err = srt_f_op_sabraf_price(forward, strike, maturity, disc, call_put,
                                  bump_vol, alpha, beta, rho, zeta,
                                  bump_vol_type, &prem1)) {
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
  if (err = srt_f_op_sabraf_price(forward, strike, maturity, disc, call_put,
                                  bump_vol + act_vega_shift, alpha, beta, rho,
                                  zeta, bump_vol_type, &prem2)) {
    return err;
  }

  /*	Calculate vega before shift	*/
  vega1 = (prem2 - prem1) / vega_shift;
  *sens = fac * vega1;

  /*	If volga is to be calculated  , calculate delta after shift	*/
  if (vega_volga == 1) {
    if (err = op_sabraf_vega_volga(
            forward, strike, maturity, disc, call_put, alpha, beta, rho, zeta,
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
Err op_SABRAF_vanna(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha, double beta, double rho, double zeta,
    /*	Vol input	*/
    double input_vol, SrtDiffusionType input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double delta_shift, double vega_shift,
    /*	Multiplicative or additive vega shift */
    /*	0: Additive  , 1: Multiplicative	*/
    int vega_mult,
    /*	Which vol to bump	*/
    SrtDiffusionType bump_vol_type,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answer	*/
    double *vanna) {
  double bump_vol;
  double delta1, delta2;
  double act_vega_shift;
  Err err;

  /*	Calculate equivalent vol of the type to be bumped	*/
  err = srt_f_optbmm2vol(forward, forward, maturity, input_vol, alpha, beta,
                         rho, zeta, input_vol_type, bump_vol_type, &bump_vol);
  if (err)
    return err;

  /*	Calculate delta before shift	*/
  if (err = op_sabraf_delta_gamma(forward, strike, maturity, disc, call_put,
                                  alpha, beta, rho, zeta, bump_vol,
                                  bump_vol_type, delta_shift, 0.0,
                                  bump_vol_type, 0, 1.0, &delta1)) {
    return err;
  }

  /*	Calculate actual shifts	*/

  if (vega_mult == 1) {
    act_vega_shift = vega_shift * bump_vol;
  } else {
    act_vega_shift = vega_shift;
  }

  /*	Calculate delta after shift	*/
  if (err = op_sabraf_delta_gamma(
          forward, strike, maturity, disc, call_put, alpha, beta, rho, zeta,
          bump_vol + act_vega_shift, bump_vol_type, delta_shift, 0.0,
          bump_vol_type, 0, 1.0, &delta2)) {
    return err;
  }

  /*	Calculate vanna	*/
  *vanna = fac * (delta2 - delta1) / vega_shift;

  return NULL;
}

/*	Theta	*/
/*	All parameters checks are supposed to have been done before call
 */
Err op_SABRAF_theta(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha, double beta, double rho, double zeta,
    /*	Vol input	*/
    double input_vol, SrtDiffusionType input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double time_shift,
    /*	Which vol to freeze when bumping	*/
    SrtDiffusionType freeze_vol_type,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answers	*/
    double *theta, double *thetacarry) {
  double frozen_vol;
  double prem1, prem2;
  double r;
  Err err;

  /*	Calculate equivalent vol of the type to be frozen	*/
  err =
      srt_f_optbmm2vol(forward, forward, maturity, input_vol, alpha, beta, rho,
                       zeta, input_vol_type, freeze_vol_type, &frozen_vol);
  if (err)
    return err;

  /*	Calculate premium before shift	*/
  if (err = srt_f_op_sabraf_price(forward, strike, maturity, disc, call_put,
                                  frozen_vol, alpha, beta, rho, zeta,
                                  freeze_vol_type, &prem1)) {
    return err;
  }

  /*	Calculate premium after shift	*/
  if (err = srt_f_op_sabraf_price(forward, strike, maturity - time_shift, disc,
                                  call_put, frozen_vol, alpha, beta, rho, zeta,
                                  freeze_vol_type, &prem2)) {
    return err;
  }

  /*	Calculate theta	*/
  *theta = (prem2 - prem1) / time_shift;

  /*	Calculate carry	*/
  r = -log(disc) / maturity;
  *thetacarry = r * prem1;
  *theta += *thetacarry;

  *theta *= fac;

  return NULL;
}

/*	Sensitivity to model parameters	*/
/*	All parameters checks are supposed to have been done before call
 */
Err op_SABRAF_model_sens(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double alpha, double beta, double rho, double zeta,
    /*	Vol input	*/
    double input_vol, SrtDiffusionType input_vol_type,
    /*	Options	*/
    /*	Size of bumps	*/
    double alpha_shift, double beta_shift, double rho_shift, double zeta_shift,
    /* mult or add */
    double alpha_mult, double beta_mult, double rho_mult, double zeta_mult,
    /*	Which vol to freeze when bumping	*/
    SrtDiffusionType freeze_vol_type,
    /*	Factor that multiplies results	*/
    double alpha_fac, double beta_fac, double rho_fac, double zeta_fac,
    /*	Answers	*/
    double *alpha_sens, double *beta_sens, double *rho_sens,
    double *zeta_sens) {
  double frozen_vol;
  double prem1, prem2, prem3, prem4, prem5;
  double act_shift;
  Err err;

  /*	Calculate equivalent vol of the type to be frozen	*/
  err =
      srt_f_optbmm2vol(forward, forward, maturity, input_vol, alpha, beta, rho,
                       zeta, input_vol_type, freeze_vol_type, &frozen_vol);
  if (err)
    return err;

  /*	Calculate premium before shift	*/
  if (err = srt_f_op_sabraf_price(forward, strike, maturity, disc, call_put,
                                  frozen_vol, alpha, beta, rho, zeta,
                                  freeze_vol_type, &prem1)) {
    return err;
  }

  /*	Calculate premium after shift of alpha	*/
  act_shift = alpha_mult == 1 ? alpha * alpha_shift : alpha_shift;
  if (err = srt_f_op_sabraf_price(forward, strike, maturity, disc, call_put,
                                  frozen_vol, alpha + act_shift, beta, rho,
                                  zeta, freeze_vol_type, &prem2)) {
    return err;
  }

  /*	Calculate premium after shift of beta	*/
  act_shift = beta_mult == 1 ? beta * beta_shift : beta_shift;
  if (err = srt_f_op_sabraf_price(forward, strike, maturity, disc, call_put,
                                  frozen_vol, alpha, beta + act_shift, rho,
                                  zeta, freeze_vol_type, &prem3)) {
    return err;
  }

  /*	Calculate premium after shift of rho	*/
  act_shift = rho_mult == 1 ? rho * rho_shift : rho_shift;
  if (err = srt_f_op_sabraf_price(forward, strike, maturity, disc, call_put,
                                  frozen_vol, alpha, beta, rho + act_shift,
                                  zeta, freeze_vol_type, &prem4)) {
    return err;
  }

  /*	Calculate premium after shift of zeta	*/
  act_shift = zeta_mult == 1 ? zeta * zeta_shift : zeta_shift;
  if (err = srt_f_op_sabraf_price(forward, strike, maturity, disc, call_put,
                                  frozen_vol, alpha, beta, rho,
                                  zeta + act_shift, freeze_vol_type, &prem5)) {
    return err;
  }

  /*	Calculate sensitivites	*/
  *alpha_sens = alpha_fac * (prem2 - prem1) / alpha_shift;
  *beta_sens = beta_fac * (prem3 - prem1) / beta_shift;
  *rho_sens = rho_fac * (prem4 - prem1) / rho_shift;
  *zeta_sens = zeta_fac * (prem5 - prem1) / zeta_shift;

  return NULL;
}

#undef R
#undef C
#undef SHIFT2
#undef SHIFT3
