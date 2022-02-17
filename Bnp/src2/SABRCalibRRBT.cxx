
#include "SABRCalibRRBT.h"
#include "math.h"
#include "opfnctns.h"

#define SABR_RRBT_DEFAULT_ALPHA 0.35

/* All the Functions for calibration of RHO on Risk Reversal */
/* ********************************************************* */

Err SABR_RRBT_GetTarget_RR(void *Inst_, void *Params_, void *Model_,
                           CALIBGEN_PARAMS CalibConsts, double *target) {
  SABR_RRBT_INST Inst = (SABR_RRBT_INST)Inst_;
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;
  SABR_RRBT_PARAMS Params = (SABR_RRBT_PARAMS)Params_;

  if (Params->solve_on_vol) {
    *target = Inst->target_vol_up - Inst->target_vol_down;
  } else {
    *target = Inst->target_price_up - Inst->target_price_down;
  }

  return NULL;
}

Err SABR_RRBT_GetFirstGuess_RR(void *Model_, void *Params_, int index_param,
                               double target, double *param1) {
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;

  *param1 = Model->rhoalpha;

  return NULL;
}

Err SABR_RRBT_GetSecondGuess_RR(void *Model_, void *Params_, int index_param,
                                double param1, double price1, double target,
                                double *param2) {
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;
  SABR_RRBT_PARAMS Params = (SABR_RRBT_PARAMS)Params_;

  if (price1 > target) {
    *param2 = param1 - (param1 + 0.99 * Model->alpha) / 4.0;

    if (fabs(*param2 + 0.99 * Model->alpha) < 1.0E-08) {
      *param2 = 0.0;
    }
  } else {
    *param2 = param1 + (0.99 * Model->alpha - param1) / 4.0;

    if (fabs(*param2 - 0.99 * Model->alpha) < 1.0E-08) {
      *param2 = 0.0;
    }
  }

  return NULL;
}

Err SABR_RRBT_GetLimitAndLastParam_RR(void *Model_, CALIBGEN_PARAMS CalibConsts,
                                      void *Params_, int index_param,
                                      double *last_param, double *limit_down,
                                      double *limit_up) {
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;
  SABR_RRBT_PARAMS Params = (SABR_RRBT_PARAMS)Params_;

  *limit_down = -0.99 * 2.0;
  *limit_up = 0.99 * 2.0;

  return NULL;
}

Err SABR_RRBT_BumpParam_RR(void *Model_, void *Params_, int index_param,
                           double param1) {
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;

  Model->rhoalpha = param1;
  Model->rho = param1 / Model->alpha;

  if (Model->rho > 0.99) {
    Model->rho = 0.99;
    Model->alpha = Model->rhoalpha / 0.99;
  } else if (Model->rho < -0.99) {
    Model->rho = -0.99;
    Model->alpha = -Model->rhoalpha / 0.99;
  }

  return NULL;
}

Err SABR_RRBT_SetParam_RR(void *Model_, CALIBGEN_PARAMS CalibConsts,
                          void *Params_, int index_param, double param1) {
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;

  Model->rhoalpha = param1;
  Model->rho = param1 / Model->alpha;

  if (Model->rho > 0.99) {
    Model->rho = 0.99;
    Model->alpha = Model->rhoalpha / 0.99;
  } else if (Model->rho < -0.99) {
    Model->rho = -0.99;
    Model->alpha = -Model->rhoalpha / 0.99;
  }

  return NULL;
}

Err SABR_RRBT_ExtrapolParam_RR(void *Model_, void *Params_, int index_param) {
  return NULL;
}

Err SABR_RRBT_UpdateConstsAfterParam_RR(void *Inst_, void *InstConst_,
                                        void *Params_, void *Model_,
                                        CALIBGEN_PARAMS CalibConsts) {
  SABR_RRBT_PARAMS Params = (SABR_RRBT_PARAMS)Params_;
  Err err = NULL;

  err = CalibrateParamTS(0, 0, &Inst_, &Params_, Params_, Model_,
                         Params->CalibParams_ForAlpha,
                         Params->AllFunctions_ForAlpha);

  return err;
}

Err SABR_RRBT_PriceInst_RR(void *Inst_, void *InstConst_, void *Params_,
                           void *Model_, double *InstPrice) {
  SABR_RRBT_PARAMS Params = (SABR_RRBT_PARAMS)Params_;
  SABR_RRBT_INST Inst = (SABR_RRBT_INST)Inst_;

  if (Params->solve_on_vol) {
    *InstPrice = Inst->vol_up - Inst->vol_down;
  } else {
    *InstPrice = Inst->price_up - Inst->price_down;
  }

  return NULL;
}

/* All the Functions for calibration of ALPHA on Butterfly */
/* ********************************************************* */

Err SABR_RRBT_GetTarget_BT(void *Inst_, void *Params_, void *Model_,
                           CALIBGEN_PARAMS CalibConsts, double *target) {
  SABR_RRBT_INST Inst = (SABR_RRBT_INST)Inst_;
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;
  SABR_RRBT_PARAMS Params = (SABR_RRBT_PARAMS)Params_;

  if (Params->solve_on_vol) {
    *target = Inst->target_vol_up + Inst->target_vol_down - 2.0 * Inst->vol_atm;
  } else {
    *target =
        Inst->target_price_up - Inst->target_price_down - 2.0 * Inst->price_atm;
  }

  return NULL;
}

Err SABR_RRBT_GetFirstGuess_BT(void *Model_, void *Params_, int index_param,
                               double target, double *param1) {
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;

  *param1 = Model->alpha * Model->alpha;

  return NULL;
}

Err SABR_RRBT_GetSecondGuess_BT(void *Model_, void *Params_, int index_param,
                                double param1, double price1, double target,
                                double *param2) {
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;
  SABR_RRBT_PARAMS Params = (SABR_RRBT_PARAMS)Params_;

  if (price1 > target) {
    *param2 = param1 - 0.2 * 0.2;
  } else {
    *param2 = param1 + 0.2 * 0.2;
  }

  return NULL;
}

Err SABR_RRBT_GetLimitAndLastParam_BT(void *Model_, CALIBGEN_PARAMS CalibConsts,
                                      void *Params_, int index_param,
                                      double *last_param, double *limit_down,
                                      double *limit_up) {
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;
  SABR_RRBT_PARAMS Params = (SABR_RRBT_PARAMS)Params_;

  *limit_down = max(pow(fabs(Model->rhoalpha) / 0.99, 2.0), 0.01 * 0.01);
  *limit_up = 2.00 * 2.00;

  return NULL;
}

Err SABR_RRBT_BumpParam_BT(void *Model_, void *Params_, int index_param,
                           double param1) {
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;

  Model->alpha = sqrt(param1);
  Model->rho = Model->rhoalpha / Model->alpha;

  return NULL;
}

Err SABR_RRBT_SetParam_BT(void *Model_, CALIBGEN_PARAMS CalibConsts,
                          void *Params_, int index_param, double param1) {
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;

  Model->alpha = sqrt(param1);
  Model->rho = Model->rhoalpha / Model->alpha;

  return NULL;
}

Err SABR_RRBT_ExtrapolParam_BT(void *Model_, void *Params_, int index_param) {
  return NULL;
}

Err SABR_RRBT_UpdateConstsAfterParam_BT(void *Inst_, void *InstConst_,
                                        void *Params_, void *Model_,
                                        CALIBGEN_PARAMS CalibConsts) {
  SABR_RRBT_PARAMS Params = (SABR_RRBT_PARAMS)Params_;
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;
  Err err = NULL;

  err = srt_f_optsarbvol(Model->fwd, Model->fwd, Model->mat, Model->vol_atm,
                         Model->alpha, Model->beta, Model->rho, Model->vol_type,
                         SRT_BETAVOL, &(Model->sigma_beta));

  if (err)
    return err;

  return err;
}

Err SABR_RRBT_PriceInst_BT(void *Inst_, void *InstConst_, void *Params_,
                           void *Model_, double *InstPrice) {
  SABR_RRBT_PARAMS Params = (SABR_RRBT_PARAMS)Params_;
  SABR_RRBT_INST Inst = (SABR_RRBT_INST)Inst_;
  SABR_RRBT_MODEL Model = (SABR_RRBT_MODEL)Model_;

  Err err = NULL;

  /* First Get the Vols */
  err =
      srt_f_optsarbvol(Model->fwd, Inst->strike_down, Model->mat,
                       Model->sigma_beta, Model->alpha, Model->beta, Model->rho,
                       SRT_BETAVOL, Model->vol_type, &(Inst->vol_down));

  if (err)
    return err;

  err = srt_f_optsarbvol(
      Model->fwd, Inst->strike_up, Model->mat, Model->sigma_beta, Model->alpha,
      Model->beta, Model->rho, SRT_BETAVOL, Model->vol_type, &(Inst->vol_up));

  if (err)
    return err;

  if (Params->solve_on_vol) {
    *InstPrice = Inst->vol_up + Inst->vol_down - 2.0 * Inst->vol_atm;
  } else {
    /* Get the price ... */
    *InstPrice = Inst->price_up - Inst->price_down - 2.0 * Inst->price_atm;
  }

  return err;
}

Err calib_sabr_rr_bt_given_beta(
    double fwd, double mat, double vol_atm, double strike_down, double vol_down,
    double strike_up, double vol_up, double *sigma_beta, double *alpha,
    double fixed_beta, double *rho, SrtDiffusionType vol_type, int solve_on_vol,
    double precision, int nb_iter_max, double *calib_err) {
  CALIBGEN_PARAMS CalibParamsRho = NULL;
  CALIBGEN_PARAMS CalibParamsAlpha = NULL;
  CALIBFUNCTIONS AllFunctionsRho = NULL;
  CALIBFUNCTIONS AllFunctionsAlpha = NULL;
  SABR_RRBT_INST Inst = NULL;
  SABR_RRBT_MODEL Model = NULL;
  SABR_RRBT_PARAMS Params = NULL;

  Err err = NULL;

  CalibParamsAlpha = calloc(1, sizeof(CALIBGEN_Params));
  CalibParamsRho = calloc(1, sizeof(CALIBGEN_Params));
  AllFunctionsAlpha = calloc(1, sizeof(CalibFunctions));
  AllFunctionsRho = calloc(1, sizeof(CalibFunctions));
  Inst = calloc(1, sizeof(sabr_rrbt_inst));
  Model = calloc(1, sizeof(sabr_rrbt_model));
  Params = calloc(1, sizeof(sabr_rrbt_params));

  if (!CalibParamsRho || !CalibParamsAlpha || !AllFunctionsAlpha ||
      !AllFunctionsRho || !Inst || !Model || !Params) {
    err = "Memory allocation faillure in calib_sabr_rr_bt_given_beta";
    goto FREE_RETURN;
  }

  /* Initialise Newton params */
  err = Initialise_CalibParams(1, precision, nb_iter_max, 1, 1, 1, 0, 0.0,
                               CalibParamsRho);

  if (err)
    goto FREE_RETURN;

  err = Initialise_CalibParams(1, precision, nb_iter_max, 1, 1, 1, 0, 0.0,
                               CalibParamsAlpha);

  if (err)
    goto FREE_RETURN;

  /* Initialise Model */
  Model->fwd = fwd;
  Model->mat = mat;
  Model->vol_atm = vol_atm;

  if (fabs(*alpha) > 1.0E-04) {
    Model->alpha = *alpha;
  } else {
    Model->alpha = SABR_RRBT_DEFAULT_ALPHA;
  }

  Model->beta = fixed_beta;
  Model->rho = *rho;
  Model->rhoalpha = Model->rho * Model->alpha;
  Model->vol_type = vol_type;

  /* Initialise Instrument */
  Inst->vol_atm = vol_atm;
  Inst->strike_down = strike_down;
  Inst->target_vol_down = vol_down;
  Inst->strike_up = strike_up;
  Inst->target_vol_up = vol_up;

  if (!solve_on_vol) {
    /* Transform into price */
  }

  /* Initialise Params */
  Params->solve_on_vol = solve_on_vol;

  /* Initialise Functions */
  AllFunctionsRho->GetTarget = SABR_RRBT_GetTarget_RR;
  AllFunctionsRho->BumpParam = SABR_RRBT_BumpParam_RR;
  AllFunctionsRho->ExtrapolParam = SABR_RRBT_ExtrapolParam_RR;
  AllFunctionsRho->GetFirstGuess = SABR_RRBT_GetFirstGuess_RR;
  AllFunctionsRho->GetLimitAndLastParam = SABR_RRBT_GetLimitAndLastParam_RR;
  AllFunctionsRho->GetSecondGuess = SABR_RRBT_GetSecondGuess_RR;
  AllFunctionsRho->PriceInst = SABR_RRBT_PriceInst_RR;
  AllFunctionsRho->SetParam = SABR_RRBT_SetParam_RR;
  AllFunctionsRho->UpdateConstsAfterParam = SABR_RRBT_UpdateConstsAfterParam_RR;

  AllFunctionsAlpha->GetTarget = SABR_RRBT_GetTarget_BT;
  AllFunctionsAlpha->BumpParam = SABR_RRBT_BumpParam_BT;
  AllFunctionsAlpha->ExtrapolParam = SABR_RRBT_ExtrapolParam_BT;
  AllFunctionsAlpha->GetFirstGuess = SABR_RRBT_GetFirstGuess_BT;
  AllFunctionsAlpha->GetLimitAndLastParam = SABR_RRBT_GetLimitAndLastParam_BT;
  AllFunctionsAlpha->GetSecondGuess = SABR_RRBT_GetSecondGuess_BT;
  AllFunctionsAlpha->PriceInst = SABR_RRBT_PriceInst_BT;
  AllFunctionsAlpha->SetParam = SABR_RRBT_SetParam_BT;
  AllFunctionsAlpha->UpdateConstsAfterParam =
      SABR_RRBT_UpdateConstsAfterParam_BT;

  Params->AllFunctions_ForAlpha = AllFunctionsAlpha;
  Params->CalibParams_ForAlpha = CalibParamsAlpha;
  Params->AllFunctions_ForRho = AllFunctionsRho;
  Params->CalibParams_ForRho = CalibParamsRho;

  /* Solve !!! */
  err = CalibrateParamTS(0, 0, &Inst, &Params, Params, Model, CalibParamsRho,
                         AllFunctionsRho);

  if (err)
    goto FREE_RETURN;

  /* Outputs */
  *sigma_beta = Model->sigma_beta;
  *alpha = Model->alpha;
  *rho = Model->rho;

  if (solve_on_vol) {
    *calib_err = pow(Inst->vol_down - Inst->target_vol_down, 2.0) +
                 pow(Inst->vol_up - Inst->target_vol_up, 2.0);
  } else {
    *calib_err = pow(Inst->price_down - Inst->target_price_down, 2.0) +
                 pow(Inst->price_up - Inst->target_price_up, 2.0);
  }

  *calib_err = sqrt(*calib_err);

FREE_RETURN:

  if (CalibParamsRho) {
    Free_CalibParams(CalibParamsRho);
    free(CalibParamsRho);
  }

  if (CalibParamsAlpha) {
    Free_CalibParams(CalibParamsAlpha);
    free(CalibParamsAlpha);
  }

  if (AllFunctionsRho)
    free(AllFunctionsRho);
  if (AllFunctionsAlpha)
    free(AllFunctionsAlpha);
  if (Inst)
    free(Inst);
  if (Model)
    free(Model);
  if (Params)
    free(Params);

  return err;
}

Err transform_sabr_beta(double mat, double fwd, double init_sigma,
                        double init_alpha, double init_beta, double init_rho,
                        SrtDiffusionType vol_type, double *new_sigmabeta,
                        double *new_alpha, double new_beta, double *new_rho,
                        double nb_std, int solve_on_vol, double precision,
                        int nb_iter_max, double *calib_err) {
  double sigma_beta, vol_atm, vol_down, vol_up, strike_down, strike_up;
  Err err = NULL;

  err = srt_f_optsarbvol(fwd, fwd, mat, init_sigma, init_alpha, init_beta,
                         init_rho, vol_type, SRT_BETAVOL, &sigma_beta);

  err = srt_f_optsarbvol(fwd, fwd, mat, sigma_beta, init_alpha, init_beta,
                         init_rho, SRT_BETAVOL, SRT_LOGNORMAL, &vol_atm);

  if (err)
    return err;

  strike_down = fwd * exp(-nb_std * vol_atm * sqrt(mat));
  strike_up = fwd * exp(nb_std * vol_atm * sqrt(mat));

  err =
      srt_f_optsarbvol(fwd, strike_down, mat, sigma_beta, init_alpha, init_beta,
                       init_rho, SRT_BETAVOL, SRT_LOGNORMAL, &vol_down);

  if (err)
    return err;

  err = srt_f_optsarbvol(fwd, strike_up, mat, sigma_beta, init_alpha, init_beta,
                         init_rho, SRT_BETAVOL, SRT_LOGNORMAL, &vol_up);

  if (err)
    return err;

  /* First Guess */
  *new_alpha = init_alpha;
  *new_rho = init_rho + (init_beta - new_beta) * vol_atm / init_alpha;

  err = srt_f_optsarbvol(fwd, fwd, mat, vol_atm, *new_alpha, new_beta, *new_rho,
                         SRT_LOGNORMAL, SRT_BETAVOL, new_sigmabeta);

  if (err)
    return err;

  err = calib_sabr_rr_bt_given_beta(
      fwd, mat, vol_atm, strike_down, vol_down, strike_up, vol_up,
      new_sigmabeta, new_alpha, new_beta, new_rho, SRT_LOGNORMAL, solve_on_vol,
      precision, nb_iter_max, calib_err);

  return err;
}