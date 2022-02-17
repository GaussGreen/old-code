/*	==========================================================================
        FILE_NAME:	Fx3FBetaDLMCalculations.c

        PURPOSE:	All the calculations in the Fx3DBetaDLM model

        DATE:		10/15/03

        AUTHOR:		L.C.
        ==========================================================================
 */

#include "Fx3FBetaDLMCalibration.h"
#include "DiagCalibGen.h"
#include "Fx3FBetaDLMCalculations.h"
#include "Fx3FBetaDLMUtil.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "Math.h"
#include "opfnctns.h"

void FxBetaDLM_Free_CalibConst(FxBetaDLM_CalibConst *CalibConst) {
  int i;

  if (CalibConst) {
    if (CalibConst->cPrecalc) {
      FxBetaDLM_Free_Precalculations(CalibConst->cPrecalc);
      free(CalibConst->cPrecalc);
      CalibConst->cPrecalc = NULL;
    }

    if (CalibConst->lSigIndex) {
      free(CalibConst->lSigIndex);
      CalibConst->lSigIndex = NULL;
    }

    if (CalibConst->cAllInstPrecalc) {
      for (i = 0; i < CalibConst->iNbOpt; i++) {
        if (CalibConst->cAllInstPrecalc[i]) {
          free(CalibConst->cAllInstPrecalc[i]);
        }
      }

      free(CalibConst->cAllInstPrecalc);
      CalibConst->cAllInstPrecalc = NULL;
    }

    if (CalibConst->CalibParams) {
      Free_CalibParams(CalibConst->CalibParams);
      free(CalibConst->CalibParams);
      CalibConst->CalibParams = NULL;
    }

    if (CalibConst->CalibParamsForBeta) {
      Free_CalibParams(CalibConst->CalibParamsForBeta);
      free(CalibConst->CalibParamsForBeta);
      CalibConst->CalibParamsForBeta = NULL;
    }

    if (CalibConst->CalibFunctionsForVol)
      free(CalibConst->CalibFunctionsForVol);
    CalibConst->CalibFunctionsForVol = NULL;

    if (CalibConst->CalibFunctionsForBeta)
      free(CalibConst->CalibFunctionsForBeta);
    CalibConst->CalibFunctionsForBeta = NULL;
  }
}

Err FxBetaDLM_Fill_CalibConst(int iNbOpt, void **cAllInst, int iCalibSmile,
                              int iSmileIndex, FxBetaDLM_model *cModel,
                              FxBetaDLM_OptNumerParams *cNumParams,
                              FxBetaDLM_Hermite *cHermite,
                              FxBetaDLM_CalibConst *cCalibConst) {
  int i;
  FxBetaDLM_FxOptInst *Inst;
  double precision;
  Err err = NULL;

  cCalibConst->cPrecalc = calloc(1, sizeof(FxBetaDLM_ModelPrecalc));
  cCalibConst->lSigIndex = calloc(iNbOpt, sizeof(long));
  cCalibConst->cAllInstPrecalc =
      calloc(iNbOpt, sizeof(FxBetaDLM_InstPrecalc *));

  if (!cCalibConst->cPrecalc || !cCalibConst->lSigIndex ||
      !cCalibConst->cAllInstPrecalc) {
    err = "Memory allocation faillure in FxBetaDLM_Fill_CalibConst";
    goto FREE_RETURN;
  }

  cCalibConst->iCalibSmile = iCalibSmile;
  cCalibConst->iNbOpt = iNbOpt;
  cCalibConst->iSmileIndex = iSmileIndex;
  cCalibConst->cHermite = cHermite;
  cCalibConst->cNumParams = cNumParams;

  err = FxBetaDLM_Allocation_Precalculations(iNbOpt, cCalibConst->cPrecalc);

  if (err)
    goto FREE_RETURN;

  for (i = 0; i < iNbOpt; i++) {
    Inst = (FxBetaDLM_FxOptInst *)(cAllInst[i]);
    cCalibConst->lSigIndex[i] =
        Get_Index(Inst->dExeTime, cModel->dPWTime, cModel->iNbPWTime);

    cCalibConst->cAllInstPrecalc[i] = calloc(1, sizeof(FxBetaDLM_InstPrecalc));

    if (!cCalibConst->cAllInstPrecalc[i]) {
      err = "Memory allocation faillure in FxBetaDLM_Fill_CalibConst";
      goto FREE_RETURN;
    }

    err = FxBetaDLM_Fill_InstPrecalc(Inst, cModel,
                                     cCalibConst->cAllInstPrecalc[i]);

    if (err)
      return err;
  }

  err = FxBetaDLM_Fill_ModelPrecalc(cModel, cCalibConst);

  if (err)
    goto FREE_RETURN;

  cCalibConst->CalibParams = calloc(1, sizeof(CALIBGEN_Params));
  cCalibConst->CalibParamsForBeta = calloc(1, sizeof(CALIBGEN_Params));

  if (!cCalibConst->CalibParams || !cCalibConst->CalibParamsForBeta) {
    err = "Memory allocation faillure in FxBetaDLM_Fill_CalibConst";
    goto FREE_RETURN;
  }

  precision = cModel->dCashFx / 100000.0;

  cCalibConst->dMaxFact = 0.001;
  cCalibConst->dMinFact = 0.001;

  err = Initialise_CalibParams(1, precision, 10, 0, 1, 0, 0, 0.0,
                               cCalibConst->CalibParams);

  if (err)
    goto FREE_RETURN;

  err = Initialise_CalibParams(1, precision, 10, 1, iCalibSmile, 1, 0, 0.0,
                               cCalibConst->CalibParamsForBeta);

  if (err)
    goto FREE_RETURN;

  cCalibConst->CalibFunctionsForVol = calloc(1, sizeof(CalibFunctions));
  cCalibConst->CalibFunctionsForBeta = calloc(1, sizeof(CalibFunctions));

  if (!cCalibConst->CalibFunctionsForVol ||
      !cCalibConst->CalibFunctionsForBeta) {
    err = "Memory allocation faillure in FxBetaDLM_CalibrationModel";
    goto FREE_RETURN;
  }

  cCalibConst->cAllInst = cAllInst;

FREE_RETURN:

  if (err) {
    FxBetaDLM_Free_CalibConst(cCalibConst);
  }

  return err;
}

Err FxBetaDLM_Calibration(/* Model informations */
                          char *undname_fx, char *undname_dom,
                          char *undname_for, double tstar, double B0,
                          double *C0, double alpha, double lambda,
                          double spot_fx,

                          /*	Correlations	*/
                          int nb_3F_corr, long *date_3F_corr,
                          double *dom_for_3F_corr, double *dom_fx_3F_corr,
                          double *for_fx_3F_corr,

                          /* Volatility informations */
                          int nb_opt, long *opt_settlmt_date,
                          long *opt_fixing_date, double *opt_strikes,
                          double *opt_bs_vols,

                          /* Smile informations */
                          int calib_smile, long smile_settlmt_date,
                          double smile_strikes[2], double smile_bs_vols[2],

                          /* Numerical parameters */
                          FxBetaDLM_OptNumerParams *NumParams,

                          /* Output */
                          double **fx_vols, int *nb_DLM_corr,
                          long **date_DLM_corr, double **dom_for_DLM_corr,
                          double **dom_fx_DLM_corr, double **for_fx_DLM_corr,
                          double **dom_fx_DLM_corr_down,
                          double **for_fx_DLM_corr_down)

{
  FXBETADLM_STR str;
  Err err = NULL;
  long nb_merge_dates;
  double *null_sigma_fx = NULL, *exercise_opt = NULL, *maturity_opt = NULL,
         *merge_dates = NULL, *sig_dom = NULL, *sig_for = NULL,
         *corr_dom_for_3F = NULL, *corr_dom_fx_3F = NULL,
         *corr_for_fx_3F = NULL, max_time_correl = 0.0,
         min_time_space_correl = 0.0;
  FxBetaDLM_model *model = NULL;
  FxBetaDLM_Hermite *hermite = NULL;

  int i, index;

  model = calloc(1, sizeof(FxBetaDLM_model));
  hermite = calloc(1, sizeof(FxBetaDLM_Hermite));

  str = calloc(1, sizeof(FxBetaDLM_str));
  null_sigma_fx = calloc(nb_opt, sizeof(double));
  exercise_opt = calloc(nb_opt, sizeof(double));
  maturity_opt = calloc(nb_opt, sizeof(double));

  if (!str || !null_sigma_fx || !model || !hermite || !maturity_opt ||
      !exercise_opt) {
    err = "Memory allocation faillure (1) in SrtInitFXBETADLMUnd";
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_opt; i++) {
    null_sigma_fx[i] = 0.1;
  }

  /* Initialisation of the model */

  err = FxBetaDLM_fill_str(
      undname_fx, undname_dom, undname_for,

      /*	FX Underlying	*/
      tstar, B0, *C0, alpha, lambda, spot_fx, nb_opt, opt_fixing_date,
      null_sigma_fx,

      /*	Correlations	*/
      /*	3 Factor correlations */
      0, /* compute 3F */
      nb_3F_corr, date_3F_corr, dom_for_3F_corr, dom_fx_3F_corr, for_fx_3F_corr,

      /*	Model correlations */
      0,                                         /* compute_DLM_corr*/
      max_time_correl, min_time_space_correl, 0, /* nb_DLM_corr */
      NULL, NULL, NULL, NULL,

      /* For Alpha case */
      NULL, NULL,

      str, NumParams);

  if (err)
    goto FREE_RETURN;

  /* Merge the rates */
  err = merge_rates_fx_corr_ts(
      str->time_dom, str->sig_dom, str->nb_dom, str->time_for, str->sig_for,
      str->nb_for, NULL, NULL, 0, str->time_3F_corr, str->dom_for_3F_corr,
      str->dom_fx_3F_corr, str->for_fx_3F_corr, str->nb_3F_corr, &merge_dates,
      &sig_dom, &sig_for, NULL, &corr_dom_for_3F, &corr_dom_fx_3F,
      &corr_for_fx_3F, &nb_merge_dates);

  if (err)
    goto FREE_RETURN;

  /* First guess */

  for (i = 0; i < nb_opt; i++) {
    exercise_opt[i] = (opt_fixing_date[i] - str->today) * YEARS_IN_DAY;
    maturity_opt[i] = (opt_settlmt_date[i] - str->today) * YEARS_IN_DAY;
  }

  err = Fx3DtsCalibration_corr(
      exercise_opt, maturity_opt, opt_bs_vols, nb_opt, merge_dates,
      nb_merge_dates, sig_dom, str->lam_dom, sig_for, str->lam_for,
      str->time_3F_corr, str->dom_for_3F_corr, str->dom_fx_3F_corr,
      str->for_fx_3F_corr, str->nb_3F_corr, fx_vols);

  if (err)
    goto FREE_RETURN;

  /*
          if(fabs(str->c0) < 10e-6)
          {
  */
  if (fabs(str->b0 - 1.0) < 10e-6) {
    memcpy(str->sig_fx, *fx_vols, nb_opt * sizeof(double));
  } else {
    err = FxBetaDLM_GetFirstGuessFromB0(
        str->b0, str->tstar, exercise_opt, maturity_opt, nb_opt, merge_dates,
        nb_merge_dates, sig_dom, str->lam_dom, sig_for, str->lam_for,
        corr_dom_for_3F, corr_dom_fx_3F, corr_for_fx_3F, *fx_vols, str->sig_fx);

    if (err)
      goto FREE_RETURN;
  }
  /*
          }
          else
          {
                  err = FxBetaDLM_GetFirstGuessFromForwardCalibration();
          }

  */

  /* Initialisation of the model */

  err = init_FxBetaDLM_model(
      NumParams->dMinTimeSpaceCorrel, NumParams->dMaxTimeCorrel, str->today,
      str->yc_dom, str->tau_dom, str->nb_dom, str->time_dom, str->sig_dom,
      str->yc_for, str->tau_for, str->nb_for, str->time_for, str->sig_for,
      str->spot_fx, str->tstar, str->b0, str->c0, str->alpha, str->lambda,
      str->nb_fx, str->time_fx, *fx_vols, str->sig_fx, str->nb_3F_corr,
      str->time_3F_corr, str->dom_for_3F_corr, str->dom_fx_3F_corr,
      str->for_fx_3F_corr, model);

  if (err)
    goto FREE_RETURN;

  err = initialise_FxBetaDLM_Hermite(NumParams, hermite);

  if (err)
    goto FREE_RETURN;

  index = 0;
  for (i = 0; i < model->iNbPWTime; i++) {
    index = Get_Index(model->dPWTime[i], exercise_opt, nb_opt);
    model->dSigmaFx[i] = str->sig_fx[index];
    if (fabs(model->dLambda) < TINY) {
      model->dSigmaFxUp[i] =
          model->dSigmaFx[i] * exp(model->dAlpha * sqrt(exercise_opt[index]));
      model->dSigmaFxDown[i] =
          model->dSigmaFx[i] * exp(-model->dAlpha * sqrt(exercise_opt[index]));
    } else {
      model->dSigmaFxUp[i] =
          model->dSigmaFx[i] *
          exp(model->dAlpha *
              sqrt((1.0 - exp(-2.0 * model->dLambda * exercise_opt[index])) /
                   (2.0 * model->dLambda)));
      model->dSigmaFxDown[i] =
          model->dSigmaFx[i] *
          exp(-model->dAlpha *
              sqrt((1.0 - exp(-2.0 * model->dLambda * exercise_opt[index])) /
                   (2.0 * model->dLambda)));
    }
  }

  err = FxBetaDLM_CalibrationModel(nb_opt, opt_settlmt_date, opt_fixing_date,
                                   opt_strikes, opt_bs_vols, calib_smile,
                                   smile_settlmt_date, smile_strikes,
                                   smile_bs_vols, model, NumParams, hermite);

  if (err)
    goto FREE_RETURN;

  for (i = 0; i < nb_opt; i++) {
    (*fx_vols)[i] = model->dSigmaFx[Get_Index(exercise_opt[i], model->dPWTime,
                                              model->iNbPWTime)];
  }

  /* allocate */
  (*date_DLM_corr) = (long *)calloc(model->iNbPWTime, sizeof(long));
  (*dom_for_DLM_corr) = (double *)calloc(model->iNbPWTime, sizeof(double));
  (*dom_fx_DLM_corr) = (double *)calloc(model->iNbPWTime, sizeof(double));
  (*for_fx_DLM_corr) = (double *)calloc(model->iNbPWTime, sizeof(double));
  (*dom_fx_DLM_corr_down) = (double *)calloc(model->iNbPWTime, sizeof(double));
  (*for_fx_DLM_corr_down) = (double *)calloc(model->iNbPWTime, sizeof(double));

  if (!dom_fx_DLM_corr || !for_fx_DLM_corr || !dom_fx_DLM_corr_down ||
      !for_fx_DLM_corr_down || !date_DLM_corr) {
    err = "Memory allocation faillure in FxBetaDLM_Calibration";
    goto FREE_RETURN;
  }

  *nb_DLM_corr = model->iNbPWTime;
  for (i = 0; i < *nb_DLM_corr; i++) {
    (*date_DLM_corr)[i] =
        (long)(model->lToday + model->dPWTime[i] * DAYS_IN_YEAR);
    (*dom_for_DLM_corr)[i] = model->dCorrDomFor[i];
    (*dom_fx_DLM_corr)[i] = model->dCorrDomFx[i];
    (*for_fx_DLM_corr)[i] = model->dCorrForFx[i];
    (*dom_fx_DLM_corr_down)[i] = model->dCorrDomFxDown[i];
    (*for_fx_DLM_corr_down)[i] = model->dCorrForFxDown[i];
  }

  *C0 = model->dC0;

  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (str) {
    FxBetaDLM_free_str(str);
    free(str);
  }

  if (null_sigma_fx) {
    free(null_sigma_fx);
  }

  if (model) {
    free_FxBetaDLM_model(model);
    free(model);
  }

  if (hermite) {
    free_FxBetaDLM_Hermite(hermite);
    free(hermite);
  }

  if (exercise_opt)
    free(exercise_opt);
  if (maturity_opt)
    free(maturity_opt);
  if (sig_dom)
    free(sig_dom);
  if (sig_for)
    free(sig_for);
  if (corr_dom_for_3F)
    free(corr_dom_for_3F);
  if (corr_dom_fx_3F)
    free(corr_dom_fx_3F);
  if (corr_for_fx_3F)
    free(corr_for_fx_3F);

  if (merge_dates)
    free(merge_dates);

  return err;
}

Err FxBetaDLM_CalibrationModel(/* Volatility informations */
                               int nb_opt, long *opt_settlmt_date,
                               long *opt_fixing_date, double *opt_strikes,
                               double *opt_bs_vols,

                               /* Smile informations */
                               int calib_smile, long smile_settlmt_date,
                               double smile_strikes[2], double smile_bs_vols[2],

                               /* Model informations */
                               FxBetaDLM_model *model,
                               FxBetaDLM_OptNumerParams *NumParams,
                               FxBetaDLM_Hermite *hermite) {
  Err err = NULL;
  int i, index_smile, nb_strike;
  void **AllInst = NULL;
  double AllStrikes[3], AllVols[3];
  FxBetaDLM_CalibConst *CalibConst = NULL;

  /* Find the index of the smile information */
  if (calib_smile) {
    index_smile = 0;
    while (index_smile < nb_opt &&
           opt_settlmt_date[index_smile] < smile_settlmt_date) {
      index_smile++;
    }

    if (index_smile == nb_opt ||
        fabs(opt_settlmt_date[index_smile] - smile_settlmt_date) > 0.5) {
      err = "Cannot calibrate smile: maturity is not included in volatility TS";
      goto FREE_RETURN;
    }
  }

  /* Setup the instruments */
  AllInst = calloc(nb_opt, sizeof(void *));

  if (!AllInst) {
    err = "Memory allocation faillure in FxBetaDLM_CalibrationModel";
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_opt; i++) {
    AllInst[i] = calloc(1, sizeof(FxBetaDLM_FxOptInst));

    if (!AllInst[i]) {
      err = "Cannot calibrate smile: maturity is not included in volatility TS";
      goto FREE_RETURN;
    }

    nb_strike = 1;

    if (calib_smile && i == index_smile) {
      nb_strike = 3;
    }

    err = FxBetaDLM_Allocate_FxOptInst(i, nb_strike, AllInst[i]);

    if (err)
      goto FREE_RETURN;

    if (calib_smile && i == index_smile) {
      AllStrikes[0] = opt_strikes[index_smile];
      AllStrikes[1] = smile_strikes[0];
      AllStrikes[2] = smile_strikes[1];

      AllVols[0] = opt_bs_vols[index_smile];
      AllVols[1] = smile_bs_vols[0];
      AllVols[2] = smile_bs_vols[1];

      err = FxBetaDLM_Setup_FxOptInst(opt_settlmt_date[i], opt_settlmt_date[i],
                                      model->lToday, opt_fixing_date[i], 3,
                                      &AllStrikes[0], &AllVols[0], SRT_PUT,
                                      model, AllInst[i]);
    } else {
      err = FxBetaDLM_Setup_FxOptInst(opt_settlmt_date[i], opt_settlmt_date[i],
                                      model->lToday, opt_fixing_date[i], 1,
                                      &opt_strikes[i], &opt_bs_vols[i], SRT_PUT,
                                      model, AllInst[i]);
    }

    if (err)
      goto FREE_RETURN;
  }

  /* Setup the calibration constants */

  CalibConst = calloc(1, sizeof(FxBetaDLM_CalibConst));

  if (!CalibConst) {
    err = "Memory allocation faillure in FxBetaDLM_CalibrationModel";
    goto FREE_RETURN;
  }

  err = FxBetaDLM_Fill_CalibConst(nb_opt, AllInst, calib_smile, index_smile,
                                  model, NumParams, hermite, CalibConst);

  if (err)
    goto FREE_RETURN;

  CalibConst->CalibFunctionsForVol->GetTarget = FxBetaDLM_GetTargetVol;
  CalibConst->CalibFunctionsForVol->BumpParam = FxBetaDLM_BumpVol;
  CalibConst->CalibFunctionsForVol->ExtrapolParam = FxBetaDLM_ExtrapolVol;
  CalibConst->CalibFunctionsForVol->GetFirstGuess = FxBetaDLM_GetFirstGuessVol;
  CalibConst->CalibFunctionsForVol->GetLimitAndLastParam =
      FxBetaDLM_GetLimitAndLastVol;
  CalibConst->CalibFunctionsForVol->GetSecondGuess =
      FxBetaDLM_GetSecondGuessVol;
  CalibConst->CalibFunctionsForVol->PriceInst = FxBetaDLM_PriceInstVol;
  CalibConst->CalibFunctionsForVol->SetParam = FxBetaDLM_SetVol;
  CalibConst->CalibFunctionsForVol->UpdateConstsAfterParam =
      FxBetaDLM_UpdateConstsAfterVol;

  CalibConst->CalibFunctionsForBeta->GetTarget = FxBetaDLM_GetTargetBeta;
  CalibConst->CalibFunctionsForBeta->BumpParam = FxBetaDLM_BumpBeta;
  CalibConst->CalibFunctionsForBeta->ExtrapolParam = FxBetaDLM_ExtrapolBeta;
  CalibConst->CalibFunctionsForBeta->GetFirstGuess =
      FxBetaDLM_GetFirstGuessBeta;
  CalibConst->CalibFunctionsForBeta->GetLimitAndLastParam =
      FxBetaDLM_GetLimitAndLastBeta;
  CalibConst->CalibFunctionsForBeta->GetSecondGuess =
      FxBetaDLM_GetSecondGuessBeta;
  CalibConst->CalibFunctionsForBeta->PriceInst = FxBetaDLM_PriceInstBeta;
  CalibConst->CalibFunctionsForBeta->SetParam = FxBetaDLM_SetBeta;
  CalibConst->CalibFunctionsForBeta->UpdateConstsAfterParam =
      FxBetaDLM_UpdateConstsAfterBeta;

  if (calib_smile) {
    CalibConst->first_guess_for_smile =
        model->dSigmaFx[CalibConst->lSigIndex[index_smile]];

    CalibConst->ExtrapolAfterSetBeta = FxBetaDLM_ExtrapolVol;
    CalibConst->CalibFunctionsForVol->ExtrapolParam =
        FxBetaDLM_ExtrapolVol_WhenBetaCalib;

    err = CalibrateParamTS(0, 0, &(AllInst[index_smile]),
                           &(CalibConst->cAllInstPrecalc[index_smile]),
                           CalibConst, model, CalibConst->CalibParamsForBeta,
                           CalibConst->CalibFunctionsForBeta);
  } else {
    err = CalibrateParamTS(0, nb_opt - 1, AllInst, CalibConst->cAllInstPrecalc,
                           CalibConst, model, CalibConst->CalibParams,
                           CalibConst->CalibFunctionsForVol);
  }

  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (AllInst) {
    for (i = 0; i < nb_opt; i++) {
      if (AllInst[i]) {
        FxBetaDLM_Free_FxOptInst(AllInst[i]);
        free(AllInst[i]);
        AllInst[i] = NULL;
      }
    }

    free(AllInst);
  }

  if (CalibConst) {
    FxBetaDLM_Free_CalibConst(CalibConst);
    free(CalibConst);
  }

  return err;
}

Err FxBetaDLM_Fill_ModelPrecalc(FxBetaDLM_model *model,
                                FxBetaDLM_CalibConst *CalibConst) {
  Err err = NULL;
  int i, start_index, end_index, index;
  double last_time;
  double var_rates_dom, var_rates_for, covar_rates;
  double sig_dom, sig_for, sig_dom2, sig_for2, sig_domfor;
  double t1, t2;

  start_index = 0;
  end_index = CalibConst->lSigIndex[CalibConst->iNbOpt - 1];
  last_time = model->dPWTime[end_index];

  var_rates_dom = 0.0;
  var_rates_for = 0.0;
  covar_rates = 0.0;

  index = 0;

  for (i = start_index; i < end_index + 1; i++) {
    if (i > start_index) {
      t1 = model->dPWTime[i - 1];
    } else {
      /* First part */
      t1 = 0.0;
    }

    if (i == end_index || start_index == end_index) {
      /* Last part */
      t2 = last_time;
    } else {
      t2 = model->dPWTime[i];
    }

    sig_dom = model->dSigmaDom[i];
    sig_for = model->dSigmaFor[i];
    sig_dom2 = sig_dom * sig_dom;
    sig_for2 = sig_for * sig_for;
    sig_domfor = model->dCorrDomFor[i] * sig_dom * sig_for;

    var_rates_dom +=
        sig_dom2 * Phi_Func(2.0 * model->dLambdaDom, model->dTStar, t1, t2);
    var_rates_for +=
        sig_for2 * Phi_Func(2.0 * model->dLambdaFor, model->dTStar, t1, t2);
    covar_rates += sig_domfor * Phi_Func(model->dLambdaDom + model->dLambdaFor,
                                         model->dTStar, t1, t2);

    if (index < CalibConst->iNbOpt &&
        fabs(t2 - model->dPWTime[CalibConst->lSigIndex[index]]) < 1.0E-10) {
      CalibConst->cPrecalc->dVarRatesDom[index] = var_rates_dom;
      CalibConst->cPrecalc->dVarRatesFor3D[index] = var_rates_for;
      CalibConst->cPrecalc->dCovarRates3D[index] = covar_rates;

      index++;
    }
  }

  return err;
}

Err FxBetaDLM_Fill_InstPrecalc(FxBetaDLM_FxOptInst *Inst,
                               FxBetaDLM_model *model,
                               FxBetaDLM_InstPrecalc *CalibConst) {
  CalibConst->dBetaDom = (1.0 - exp(-model->dLambdaDom *
                                    (Inst->dSettlementTime - model->dTStar))) /
                         model->dLambdaDom;
  CalibConst->dBetaFor = (1.0 - exp(-model->dLambdaFor *
                                    (Inst->dSettlementTime - model->dTStar))) /
                         model->dLambdaFor;

  return NULL;
}

Err FxBetaDLM_Calculate_AllConst(FxBetaDLM_FxOptInst *Inst,
                                 FxBetaDLM_model *model,
                                 FxBetaDLM_OptNumerParams *NumParams,
                                 FxBetaDLM_InstPrecalc *InstConst) {
  Err err = NULL;
  int j, k;
  double t1, t2, dt;
  double for_fwd;
  double var_rates_dom, var_rates_for, var_ffx;
  double covar_dom_for, covar_dom_ffx, covar_for_ffx;
  double covar_dom_for_adjust;
  double sig_dom, sig_for, sig_fx, sig_dom2, sig_for2, sig_fx2, sig_domfor,
      sig_domfx, sig_forfx;
  double dom_lam, for_lam, tstar;

  int start_index, end_index;

  int nb_points;
  double min_dt, old_t, new_t;
  double temp_dom, temp_for, temp_bracket;
  double expect_for, total_var_ffx, var_0_fwdstart;
  double old_integ, new_integ, old_integ2, new_integ2;

  /* Parameters for the alpha fudge */
  double temp_bracket_down;
  double for_fwd_down;
  double var_ffx_down;
  double covar_dom_ffx_down, covar_for_ffx_down;
  double covar_dom_for_adjust_down;
  double sig_fx_down, sig_fx2_down, sig_domfx_down, sig_forfx_down;
  double expect_for_down, total_var_ffx_down, var_0_fwdstart_down;
  double old_integ_down, new_integ_down, old_integ2_down, new_integ2_down;

  dom_lam = model->dLambdaDom;
  for_lam = model->dLambdaFor;
  tstar = model->dTStar;

  start_index =
      Get_Index(Inst->dFwdStartTime, model->dPWTime, model->iNbPWTime);
  end_index = Get_Index(Inst->dExeTime, model->dPWTime, model->iNbPWTime);

  var_rates_dom = 0.0;
  var_rates_for = 0.0;
  covar_dom_for = 0.0;
  covar_dom_for_adjust = 0.0;
  var_ffx = 0.0;
  covar_dom_ffx = 0.0;
  covar_for_ffx = 0.0;
  expect_for = 0.0;
  old_integ = 0.0;
  old_integ2 = 0.0;
  old_t = 0.0;
  total_var_ffx = 0.0;

  covar_dom_for_adjust_down = 0.0;
  var_ffx_down = 0.0;
  covar_dom_ffx_down = 0.0;
  covar_for_ffx_down = 0.0;
  expect_for_down = 0.0;
  old_integ_down = 0.0;
  old_integ2_down = 0.0;
  total_var_ffx_down = 0.0;

  /* First update the total_var_ffx for C(t) */
  if (fabs(model->dAlpha) < TINY) {
    for (j = 0; j < start_index + 1; j++) {
      if (j > 0) {
        t1 = model->dPWTime[j - 1];
      } else {
        /* First part */
        t1 = 0.0;
      }

      if (j == start_index || 0 == start_index) {
        /* Last part */
        t2 = Inst->dFwdStartTime;
      } else {
        t2 = model->dPWTime[j];
      }

      dt = t2 - t1;

      if (dt > 0.0) {
        sig_dom = model->dSigmaDom[j];
        sig_dom2 = sig_dom * sig_dom;
        sig_for = model->dSigmaFor[j];
        sig_for2 = sig_for * sig_for;
        sig_fx = model->dSigmaFx[j];
        sig_fx2 = sig_fx * sig_fx;
        sig_domfor = model->dCorrDomFor[j] * sig_dom * sig_for;
        sig_domfx = model->dCorrDomFx[j] * sig_dom * sig_fx;
        sig_forfx = model->dCorrForFx[j] * sig_for * sig_fx;
        total_var_ffx += sig_fx2 * dt;
        total_var_ffx += -2.0 * sig_forfx * Etha_Func(for_lam, tstar, t1, t2) +
                         2.0 * sig_domfx * Etha_Func(dom_lam, tstar, t1, t2);
        total_var_ffx +=
            sig_for2 * Psi_Func(for_lam, for_lam, tstar, t1, t2) +
            sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, t1, t2) -
            2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, tstar, t1, t2);
      }
    }

    var_0_fwdstart = total_var_ffx;

    /* Calculate the starting points */
    old_t = Inst->dFwdStartTime;

    temp_dom = exp(-dom_lam * (model->dTStar - old_t));
    temp_for = exp(-for_lam * (model->dTStar - old_t));

    temp_bracket =
        temp_for * (sig_forfx - sig_for2 * (1.0 - temp_for) / for_lam +
                    sig_domfor * (1.0 - temp_dom) / dom_lam);

    old_integ = temp_bracket / (1.0 + 2.0 * model->dC0 * total_var_ffx);
    old_integ2 = covar_dom_ffx * old_integ;
    old_integ *= total_var_ffx;

    /* then caluclate the variances between fwd_start and maturity */
    for (j = start_index; j < end_index + 1; j++) {
      if (j > start_index) {
        t1 = model->dPWTime[j - 1];
      } else {
        /* First part */
        t1 = Inst->dFwdStartTime;
      }

      if (j == end_index || start_index == end_index) {
        /* Last part */
        t2 = Inst->dExeTime;
      } else {
        t2 = model->dPWTime[j];
      }

      dt = t2 - t1;

      if (dt > 0.0) {
        sig_dom = model->dSigmaDom[j];
        sig_for = model->dSigmaFor[j];
        sig_fx = model->dSigmaFx[j];
        sig_dom2 = sig_dom * sig_dom;
        sig_for2 = sig_for * sig_for;
        sig_fx2 = sig_fx * sig_fx;
        sig_domfor = model->dCorrDomFor[j] * sig_dom * sig_for;
        sig_domfx = model->dCorrDomFx[j] * sig_dom * sig_fx;
        sig_forfx = model->dCorrForFx[j] * sig_for * sig_fx;

        /* Expectation and Variance approximation for foreign */
        nb_points = max((int)(dt / NumParams->dMinTime + 0.5), 1);
        min_dt = dt / nb_points;

        for (k = 0; k < nb_points; k++) {
          new_t = old_t + min_dt;

          /* update var_ffx */
          total_var_ffx += sig_fx2 * min_dt;
          total_var_ffx +=
              -2.0 * sig_forfx * Etha_Func(for_lam, tstar, old_t, new_t) +
              2.0 * sig_domfx * Etha_Func(dom_lam, tstar, old_t, new_t);
          total_var_ffx +=
              sig_for2 * Psi_Func(for_lam, for_lam, tstar, old_t, new_t) +
              sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, old_t, new_t) -
              2.0 * sig_domfor *
                  Psi_Func(dom_lam, for_lam, tstar, old_t, new_t);

          covar_dom_ffx +=
              sig_domfx * Phi_Func(dom_lam, tstar, old_t, new_t) -
              sig_domfor * Gamma_Func(for_lam, dom_lam, tstar, old_t, new_t) +
              sig_dom2 * Gamma_Func(dom_lam, dom_lam, tstar, old_t, new_t);

          temp_dom = exp(-dom_lam * (model->dTStar - new_t));
          temp_for = exp(-for_lam * (model->dTStar - new_t));

          temp_bracket =
              temp_for * (sig_forfx - sig_for2 * (1.0 - temp_for) / for_lam +
                          sig_domfor * (1.0 - temp_dom) / dom_lam);

          new_integ = temp_bracket / (1.0 + 2.0 * model->dC0 * total_var_ffx);
          new_integ2 = covar_dom_ffx * new_integ;
          new_integ *= total_var_ffx;

          expect_for += (old_integ + new_integ) * min_dt / 2.0;
          covar_dom_for_adjust += (old_integ2 + new_integ2) * min_dt / 2.0;

          old_t = new_t;
          old_integ = new_integ;
          old_integ2 = new_integ2;
        }

        /* Update */
        var_rates_dom += sig_dom2 * Phi_Func(2.0 * dom_lam, tstar, t1, t2);
        var_rates_for += sig_for2 * Phi_Func(2.0 * for_lam, tstar, t1, t2);
        covar_dom_for +=
            sig_domfor * Phi_Func(dom_lam + for_lam, tstar, t1, t2);

        /* Fx component */
        var_ffx += sig_fx2 * dt;

        /* Fx / Rates component */
        var_ffx += -2.0 * sig_forfx * Etha_Func(for_lam, tstar, t1, t2) +
                   2.0 * sig_domfx * Etha_Func(dom_lam, tstar, t1, t2);

        /* Rates / Rates component */
        var_ffx += sig_for2 * Psi_Func(for_lam, for_lam, tstar, t1, t2) +
                   sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, t1, t2) -
                   2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, tstar, t1, t2);

        covar_for_ffx +=
            sig_forfx * Phi_Func(for_lam, tstar, t1, t2) -
            sig_for2 * Gamma_Func(for_lam, for_lam, tstar, t1, t2) +
            sig_domfor * Gamma_Func(dom_lam, for_lam, tstar, t1, t2);
      }
    }

    for_fwd = -covar_for_ffx + 2.0 * model->dC0 * expect_for;

    /* Save the coefficients */

    err = FxBetaDLM_Fill_InstPrecalc(Inst, model, InstConst);

    if (err)
      return err;

    InstConst->dVarRatesDom = var_rates_dom;
    InstConst->dVarRatesFor =
        var_rates_for - 2.0 * model->dC0 *
                            (1.0 + 2.0 * model->dC0 * var_0_fwdstart) *
                            for_fwd * for_fwd;
    InstConst->dCovarRates =
        covar_dom_for - 2.0 * model->dC0 * covar_dom_for_adjust;

    InstConst->dTotalVarFFx = var_ffx + var_0_fwdstart;
    InstConst->dVarFFx = var_ffx;
    InstConst->dStdFFX = sqrt(var_ffx);
    InstConst->dCovarFFxDom = covar_dom_ffx;
    InstConst->dCovarFFxFor =
        -(1.0 + 2.0 * model->dC0 * var_0_fwdstart) * for_fwd;
  } else /* Fill the constants with alpha fudge */
  {
    for (j = 0; j < start_index + 1; j++) {
      if (j > 0) {
        t1 = model->dPWTime[j - 1];
      } else {
        /* First part */
        t1 = 0.0;
      }

      if (j == start_index || 0 == start_index) {
        /* Last part */
        t2 = Inst->dFwdStartTime;
      } else {
        t2 = model->dPWTime[j];
      }

      dt = t2 - t1;

      if (dt > 0.0) {
        sig_dom = model->dSigmaDom[j];
        sig_dom2 = sig_dom * sig_dom;
        sig_for = model->dSigmaFor[j];
        sig_for2 = sig_for * sig_for;
        sig_fx = model->dSigmaFxUp[j];
        sig_fx2 = sig_fx * sig_fx;
        sig_domfor = model->dCorrDomFor[j] * sig_dom * sig_for;
        sig_domfx = model->dCorrDomFx[j] * sig_dom * sig_fx;
        sig_forfx = model->dCorrForFx[j] * sig_for * sig_fx;
        total_var_ffx += sig_fx2 * dt;
        total_var_ffx += -2.0 * sig_forfx * Etha_Func(for_lam, tstar, t1, t2) +
                         2.0 * sig_domfx * Etha_Func(dom_lam, tstar, t1, t2);
        total_var_ffx +=
            sig_for2 * Psi_Func(for_lam, for_lam, tstar, t1, t2) +
            sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, t1, t2) -
            2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, tstar, t1, t2);

        sig_fx_down = model->dSigmaFxDown[j];
        sig_fx2_down = sig_fx_down * sig_fx_down;
        sig_domfx_down = model->dCorrDomFxDown[j] * sig_dom * sig_fx_down;
        sig_forfx_down = model->dCorrForFxDown[j] * sig_for * sig_fx_down;
        total_var_ffx_down += sig_fx2_down * dt;
        total_var_ffx_down +=
            -2.0 * sig_forfx_down * Etha_Func(for_lam, tstar, t1, t2) +
            2.0 * sig_domfx_down * Etha_Func(dom_lam, tstar, t1, t2);
        total_var_ffx_down +=
            sig_for2 * Psi_Func(for_lam, for_lam, tstar, t1, t2) +
            sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, t1, t2) -
            2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, tstar, t1, t2);
      }
    }

    var_0_fwdstart = total_var_ffx;
    var_0_fwdstart_down = total_var_ffx_down;

    /* Calculate the starting points */
    old_t = Inst->dFwdStartTime;

    temp_dom = exp(-dom_lam * (model->dTStar - old_t));
    temp_for = exp(-for_lam * (model->dTStar - old_t));

    temp_bracket =
        temp_for * (sig_forfx - sig_for2 * (1.0 - temp_for) / for_lam +
                    sig_domfor * (1.0 - temp_dom) / dom_lam);
    temp_bracket_down =
        temp_for * (sig_forfx_down - sig_for2 * (1.0 - temp_for) / for_lam +
                    sig_domfor * (1.0 - temp_dom) / dom_lam);

    old_integ = temp_bracket / (1.0 + 2.0 * model->dC0 * total_var_ffx);
    old_integ2 = covar_dom_ffx * old_integ;
    old_integ *= total_var_ffx;

    old_integ_down =
        temp_bracket_down / (1.0 + 2.0 * model->dC0 * total_var_ffx_down);
    old_integ2_down = covar_dom_ffx_down * old_integ_down;
    old_integ_down *= total_var_ffx_down;

    /* then caluclate the variances between fwd_start and maturity */
    for (j = start_index; j < end_index + 1; j++) {
      if (j > start_index) {
        t1 = model->dPWTime[j - 1];
      } else {
        /* First part */
        t1 = Inst->dFwdStartTime;
      }

      if (j == end_index || start_index == end_index) {
        /* Last part */
        t2 = Inst->dExeTime;
      } else {
        t2 = model->dPWTime[j];
      }

      dt = t2 - t1;

      if (dt > 0.0) {
        sig_dom = model->dSigmaDom[j];
        sig_for = model->dSigmaFor[j];
        sig_dom2 = sig_dom * sig_dom;
        sig_for2 = sig_for * sig_for;
        sig_domfor = model->dCorrDomFor[j] * sig_dom * sig_for;

        sig_fx = model->dSigmaFxUp[j];
        sig_fx2 = sig_fx * sig_fx;
        sig_domfx = model->dCorrDomFx[j] * sig_dom * sig_fx;
        sig_forfx = model->dCorrForFx[j] * sig_for * sig_fx;

        sig_fx_down = model->dSigmaFxDown[j];
        sig_fx2_down = sig_fx_down * sig_fx_down;
        sig_domfx_down = model->dCorrDomFxDown[j] * sig_dom * sig_fx_down;
        sig_forfx_down = model->dCorrForFxDown[j] * sig_for * sig_fx_down;

        /* Expectation and Variance approximation for foreign */
        nb_points = max((int)(dt / NumParams->dMinTime + 0.5), 1);
        min_dt = dt / nb_points;

        for (k = 0; k < nb_points; k++) {
          new_t = old_t + min_dt;
          temp_dom = exp(-dom_lam * (model->dTStar - new_t));
          temp_for = exp(-for_lam * (model->dTStar - new_t));

          /* update var_ffx up*/
          total_var_ffx += sig_fx2 * min_dt;
          total_var_ffx +=
              -2.0 * sig_forfx * Etha_Func(for_lam, tstar, old_t, new_t) +
              2.0 * sig_domfx * Etha_Func(dom_lam, tstar, old_t, new_t);
          total_var_ffx +=
              sig_for2 * Psi_Func(for_lam, for_lam, tstar, old_t, new_t) +
              sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, old_t, new_t) -
              2.0 * sig_domfor *
                  Psi_Func(dom_lam, for_lam, tstar, old_t, new_t);

          covar_dom_ffx +=
              sig_domfx * Phi_Func(dom_lam, tstar, old_t, new_t) -
              sig_domfor * Gamma_Func(for_lam, dom_lam, tstar, old_t, new_t) +
              sig_dom2 * Gamma_Func(dom_lam, dom_lam, tstar, old_t, new_t);

          temp_bracket =
              temp_for * (sig_forfx - sig_for2 * (1.0 - temp_for) / for_lam +
                          sig_domfor * (1.0 - temp_dom) / dom_lam);

          new_integ = temp_bracket / (1.0 + 2.0 * model->dC0 * total_var_ffx);
          new_integ2 = covar_dom_ffx * new_integ;
          new_integ *= total_var_ffx;

          expect_for += (old_integ + new_integ) * min_dt / 2.0;
          covar_dom_for_adjust += (old_integ2 + new_integ2) * min_dt / 2.0;

          old_integ = new_integ;
          old_integ2 = new_integ2;

          /* update var_ffx down*/
          total_var_ffx_down += sig_fx2_down * min_dt;
          total_var_ffx_down +=
              -2.0 * sig_forfx_down * Etha_Func(for_lam, tstar, old_t, new_t) +
              2.0 * sig_domfx_down * Etha_Func(dom_lam, tstar, old_t, new_t);
          total_var_ffx_down +=
              sig_for2 * Psi_Func(for_lam, for_lam, tstar, old_t, new_t) +
              sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, old_t, new_t) -
              2.0 * sig_domfor *
                  Psi_Func(dom_lam, for_lam, tstar, old_t, new_t);

          covar_dom_ffx_down +=
              sig_domfx_down * Phi_Func(dom_lam, tstar, old_t, new_t) -
              sig_domfor * Gamma_Func(for_lam, dom_lam, tstar, old_t, new_t) +
              sig_dom2 * Gamma_Func(dom_lam, dom_lam, tstar, old_t, new_t);

          temp_bracket_down =
              temp_for *
              (sig_forfx_down - sig_for2 * (1.0 - temp_for) / for_lam +
               sig_domfor * (1.0 - temp_dom) / dom_lam);

          new_integ_down =
              temp_bracket_down / (1.0 + 2.0 * model->dC0 * total_var_ffx_down);
          new_integ2_down = covar_dom_ffx_down * new_integ_down;
          new_integ_down *= total_var_ffx_down;

          expect_for_down += (old_integ_down + new_integ_down) * min_dt / 2.0;
          covar_dom_for_adjust_down +=
              (old_integ2_down + new_integ2_down) * min_dt / 2.0;

          old_t = new_t;
          old_integ_down = new_integ_down;
          old_integ2_down = new_integ2_down;
        }

        /* Update */
        var_rates_dom += sig_dom2 * Phi_Func(2.0 * dom_lam, tstar, t1, t2);
        var_rates_for += sig_for2 * Phi_Func(2.0 * for_lam, tstar, t1, t2);
        covar_dom_for +=
            sig_domfor * Phi_Func(dom_lam + for_lam, tstar, t1, t2);

        /* Fx Up component */
        var_ffx += sig_fx2 * dt;

        /* Fx Up / Rates component */
        var_ffx += -2.0 * sig_forfx * Etha_Func(for_lam, tstar, t1, t2) +
                   2.0 * sig_domfx * Etha_Func(dom_lam, tstar, t1, t2);

        /* Rates Up / Rates component */
        var_ffx += sig_for2 * Psi_Func(for_lam, for_lam, tstar, t1, t2) +
                   sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, t1, t2) -
                   2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, tstar, t1, t2);

        covar_for_ffx +=
            sig_forfx * Phi_Func(for_lam, tstar, t1, t2) -
            sig_for2 * Gamma_Func(for_lam, for_lam, tstar, t1, t2) +
            sig_domfor * Gamma_Func(dom_lam, for_lam, tstar, t1, t2);

        /* Fx Down component */
        var_ffx_down += sig_fx2_down * dt;

        /* Fx Down / Rates component */
        var_ffx_down +=
            -2.0 * sig_forfx_down * Etha_Func(for_lam, tstar, t1, t2) +
            2.0 * sig_domfx_down * Etha_Func(dom_lam, tstar, t1, t2);

        /* Rates Down / Rates component */
        var_ffx_down +=
            sig_for2 * Psi_Func(for_lam, for_lam, tstar, t1, t2) +
            sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, t1, t2) -
            2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, tstar, t1, t2);

        covar_for_ffx_down +=
            sig_forfx_down * Phi_Func(for_lam, tstar, t1, t2) -
            sig_for2 * Gamma_Func(for_lam, for_lam, tstar, t1, t2) +
            sig_domfor * Gamma_Func(dom_lam, for_lam, tstar, t1, t2);
      }
    }

    for_fwd = -covar_for_ffx + 2.0 * model->dC0 * expect_for;
    for_fwd_down = -covar_for_ffx_down + 2.0 * model->dC0 * expect_for_down;

    /* Save the coefficients */

    err = FxBetaDLM_Fill_InstPrecalc(Inst, model, InstConst);

    if (err)
      return err;

    InstConst->dVarRatesDom = var_rates_dom;

    /* Save Up */
    InstConst->dVarRatesFor =
        var_rates_for - 2.0 * model->dC0 *
                            (1.0 + 2.0 * model->dC0 * var_0_fwdstart) *
                            for_fwd * for_fwd;
    InstConst->dCovarRates =
        covar_dom_for - 2.0 * model->dC0 * covar_dom_for_adjust;

    InstConst->dTotalVarFFx = var_ffx + var_0_fwdstart;
    InstConst->dVarFFx = var_ffx;
    InstConst->dStdFFX = sqrt(var_ffx);
    InstConst->dCovarFFxDom = covar_dom_ffx;
    InstConst->dCovarFFxFor =
        -(1.0 + 2.0 * model->dC0 * var_0_fwdstart) * for_fwd;

    /* Save Down */
    InstConst->dVarRatesFor_down =
        var_rates_for - 2.0 * model->dC0 *
                            (1.0 + 2.0 * model->dC0 * var_0_fwdstart_down) *
                            for_fwd_down * for_fwd_down;
    InstConst->dCovarRates_down =
        covar_dom_for - 2.0 * model->dC0 * covar_dom_for_adjust_down;

    InstConst->dTotalVarFFx_down = var_ffx_down + var_0_fwdstart_down;
    InstConst->dVarFFx_down = var_ffx_down;
    InstConst->dStdFFX_down = sqrt(var_ffx_down);
    InstConst->dCovarFFxDom_down = covar_dom_ffx_down;
    InstConst->dCovarFFxFor_down =
        -(1.0 + 2.0 * model->dC0 * var_0_fwdstart_down) * for_fwd_down;
  }
  return err;
}

Err FxBetaDLM_Update_ModelPrecalc(int index, FxBetaDLM_model *model,
                                  FxBetaDLM_CalibConst *CalibConst) {
  Err err = NULL;
  int j, k;
  double t1, t2, dt, new_time, last_time;
  double for_fwd;
  double var_ffx;
  double covar_dom_ffx, covar_for_ffx;
  double sig_dom, sig_for, sig_fx, sig_dom2, sig_for2, sig_fx2, sig_domfor,
      sig_domfx, sig_forfx;
  double dom_lam, for_lam, tstar;

  int start_index, end_index;

  int nb_points;
  double min_dt, old_t, new_t;
  double temp_dom, temp_for, temp_bracket;
  double expect_for;
  double old_integ, new_integ, old_integ2, new_integ2;
  double covar_dom_for_adjust;

  double sig_fx_down, sig_fx2_down, sig_domfx_down, sig_forfx_down,
      var_ffx_down, covar_dom_ffx_down, covar_for_ffx_down, expect_for_down,
      old_integ_down, new_integ_down, old_integ2_down, new_integ2_down,
      covar_dom_for_adjust_down, temp_bracket_down, for_fwd_down;

  double sig_fx_mid, sig_fx2_mid, sig_domfx_mid, sig_forfx_mid, var_ffx_mid,
      covar_dom_ffx_mid, covar_for_ffx_mid, expect_for_mid, old_integ_mid,
      new_integ_mid, old_integ2_mid, new_integ2_mid, covar_dom_for_adjust_mid,
      temp_bracket_mid, for_fwd_mid;

  dom_lam = model->dLambdaDom;
  for_lam = model->dLambdaFor;
  tstar = model->dTStar;

  if (fabs(model->dAlpha) < TINY) {
    if (index == 0) {
      start_index = 0;
      last_time = 0.0;
      var_ffx = 0.0;
      covar_dom_ffx = 0.0;
      covar_for_ffx = 0.0;
      expect_for = 0.0;
      old_t = 0.0;
      old_integ = 0.0;
      old_integ2 = 0.0;
      covar_dom_for_adjust = 0.0;
    } else {
      start_index = CalibConst->lSigIndex[index - 1];
      last_time = model->dPWTime[start_index];
      var_ffx = CalibConst->cPrecalc->dVarFFx[index - 1];
      covar_dom_ffx = CalibConst->cPrecalc->dCovarFFxDom[index - 1];
      covar_for_ffx = CalibConst->cPrecalc->dCovarFFxFor3D[index - 1];
      expect_for = CalibConst->cPrecalc->dExpectFor[index - 1];
      old_t = last_time;
      old_integ = CalibConst->cPrecalc->dIntegral[index - 1];
      old_integ2 = CalibConst->cPrecalc->dIntegral2[index - 1];
      covar_dom_for_adjust = CalibConst->cPrecalc->dCovarRatesAdjust[index - 1];
    }

    end_index = CalibConst->lSigIndex[index];
    new_time = model->dPWTime[end_index];

    for (j = start_index; j < end_index + 1; j++) {
      if (j > start_index) {
        t1 = model->dPWTime[j - 1];
      } else {
        /* First part */
        t1 = last_time;
      }

      if (j == end_index || start_index == end_index) {
        /* Last part */
        t2 = new_time;
      } else {
        t2 = model->dPWTime[j];
      }

      dt = t2 - t1;

      if (dt > 0.0) {
        sig_dom = model->dSigmaDom[j];
        sig_for = model->dSigmaFor[j];
        sig_fx = model->dSigmaFx[j];
        sig_dom2 = sig_dom * sig_dom;
        sig_for2 = sig_for * sig_for;
        sig_fx2 = sig_fx * sig_fx;
        sig_domfor = model->dCorrDomFor[j] * sig_dom * sig_for;
        sig_domfx = model->dCorrDomFx[j] * sig_dom * sig_fx;
        sig_forfx = model->dCorrForFx[j] * sig_for * sig_fx;

        /* Expectation and Variance approximation for foreign */
        nb_points = max((int)(dt / CalibConst->cNumParams->dMinTime + 0.5), 1);
        min_dt = dt / nb_points;

        for (k = 0; k < nb_points; k++) {
          new_t = old_t + min_dt;

          /* update var_ffx */
          var_ffx += sig_fx2 * min_dt;
          var_ffx +=
              -2.0 * sig_forfx * Etha_Func(for_lam, tstar, old_t, new_t) +
              2.0 * sig_domfx * Etha_Func(dom_lam, tstar, old_t, new_t);
          var_ffx +=
              sig_for2 * Psi_Func(for_lam, for_lam, tstar, old_t, new_t) +
              sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, old_t, new_t) -
              2.0 * sig_domfor *
                  Psi_Func(dom_lam, for_lam, tstar, old_t, new_t);

          covar_dom_ffx +=
              sig_domfx * Phi_Func(dom_lam, tstar, old_t, new_t) -
              sig_domfor * Gamma_Func(for_lam, dom_lam, tstar, old_t, new_t) +
              sig_dom2 * Gamma_Func(dom_lam, dom_lam, tstar, old_t, new_t);

          temp_dom = exp(-dom_lam * (model->dTStar - new_t));
          temp_for = exp(-for_lam * (model->dTStar - new_t));

          temp_bracket =
              temp_for * (sig_forfx - sig_for2 * (1.0 - temp_for) / for_lam +
                          sig_domfor * (1.0 - temp_dom) / dom_lam);

          new_integ = temp_bracket / (1.0 + 2.0 * model->dC0 * var_ffx);
          new_integ2 = covar_dom_ffx * new_integ;
          new_integ *= var_ffx;

          expect_for += (old_integ + new_integ) * min_dt / 2.0;
          covar_dom_for_adjust += (old_integ2 + new_integ2) * min_dt / 2.0;

          old_t = new_t;
          old_integ = new_integ;
          old_integ2 = new_integ2;
        }

        covar_for_ffx +=
            sig_forfx * Phi_Func(for_lam, tstar, t1, t2) -
            sig_for2 * Gamma_Func(for_lam, for_lam, tstar, t1, t2) +
            sig_domfor * Gamma_Func(dom_lam, for_lam, tstar, t1, t2);
      }
    }

    for_fwd = -covar_for_ffx + 2.0 * model->dC0 * expect_for;

    /* Update the saved coefficients */
    if (CalibConst->cPrecalc) {
      CalibConst->cPrecalc->dVarFFx[CalibConst->cPrecalc->iNbDone] = var_ffx;
      CalibConst->cPrecalc->dCovarFFxDom[CalibConst->cPrecalc->iNbDone] =
          covar_dom_ffx;
      CalibConst->cPrecalc->dCovarFFxFor[CalibConst->cPrecalc->iNbDone] =
          -for_fwd;
      CalibConst->cPrecalc->dCovarFFxFor3D[CalibConst->cPrecalc->iNbDone] =
          covar_for_ffx;
      CalibConst->cPrecalc->dExpectFor[CalibConst->cPrecalc->iNbDone] =
          expect_for;
      CalibConst->cPrecalc->dIntegral[CalibConst->cPrecalc->iNbDone] =
          new_integ;
      CalibConst->cPrecalc->dIntegral2[CalibConst->cPrecalc->iNbDone] =
          new_integ2;
      CalibConst->cPrecalc->dCovarRatesAdjust[CalibConst->cPrecalc->iNbDone] =
          covar_dom_for_adjust;

      /* Add the adjustments to the FORVAR and RATESCOV*/
      CalibConst->cPrecalc->dVarRatesFor[CalibConst->cPrecalc->iNbDone] =
          CalibConst->cPrecalc->dVarRatesFor3D[CalibConst->cPrecalc->iNbDone] -
          2.0 * model->dC0 * for_fwd * for_fwd;
      CalibConst->cPrecalc->dCovarRates[CalibConst->cPrecalc->iNbDone] =
          CalibConst->cPrecalc->dCovarRates3D[CalibConst->cPrecalc->iNbDone] -
          2.0 * model->dC0 * covar_dom_for_adjust;
    }
  } else {
    if (index == 0) {
      start_index = 0;
      last_time = 0.0;
      old_t = 0.0;

      var_ffx = 0.0;
      covar_dom_ffx = 0.0;
      covar_for_ffx = 0.0;
      expect_for = 0.0;
      old_integ = 0.0;
      old_integ2 = 0.0;
      covar_dom_for_adjust = 0.0;

      var_ffx_down = 0.0;
      covar_dom_ffx_down = 0.0;
      covar_for_ffx_down = 0.0;
      expect_for_down = 0.0;
      old_integ_down = 0.0;
      old_integ2_down = 0.0;
      covar_dom_for_adjust_down = 0.0;

      var_ffx_mid = 0.0;
      covar_dom_ffx_mid = 0.0;
      covar_for_ffx_mid = 0.0;
      expect_for_mid = 0.0;
      old_integ_mid = 0.0;
      old_integ2_mid = 0.0;
      covar_dom_for_adjust_mid = 0.0;
    } else {
      start_index = CalibConst->lSigIndex[index - 1];
      last_time = model->dPWTime[start_index];
      old_t = last_time;

      var_ffx = CalibConst->cPrecalc->dVarFFx[index - 1];
      covar_dom_ffx = CalibConst->cPrecalc->dCovarFFxDom[index - 1];
      covar_for_ffx = CalibConst->cPrecalc->dCovarFFxFor3D[index - 1];
      expect_for = CalibConst->cPrecalc->dExpectFor[index - 1];
      old_integ = CalibConst->cPrecalc->dIntegral[index - 1];
      old_integ2 = CalibConst->cPrecalc->dIntegral2[index - 1];
      covar_dom_for_adjust = CalibConst->cPrecalc->dCovarRatesAdjust[index - 1];

      var_ffx_down = CalibConst->cPrecalc->dVarFFx_down[index - 1];
      covar_dom_ffx_down = CalibConst->cPrecalc->dCovarFFxDom_down[index - 1];
      covar_for_ffx_down = CalibConst->cPrecalc->dCovarFFxFor3D_down[index - 1];
      expect_for_down = CalibConst->cPrecalc->dExpectFor_down[index - 1];
      old_integ_down = CalibConst->cPrecalc->dIntegral_down[index - 1];
      old_integ2_down = CalibConst->cPrecalc->dIntegral2_down[index - 1];
      covar_dom_for_adjust_down =
          CalibConst->cPrecalc->dCovarRatesAdjust_down[index - 1];

      var_ffx_mid = CalibConst->cPrecalc->dVarFFx_mid[index - 1];
      covar_dom_ffx_mid = CalibConst->cPrecalc->dCovarFFxDom_mid[index - 1];
      covar_for_ffx_mid = CalibConst->cPrecalc->dCovarFFxFor3D_mid[index - 1];
      expect_for_mid = CalibConst->cPrecalc->dExpectFor_mid[index - 1];
      old_integ_mid = CalibConst->cPrecalc->dIntegral_mid[index - 1];
      old_integ2_mid = CalibConst->cPrecalc->dIntegral2_mid[index - 1];
      covar_dom_for_adjust_mid =
          CalibConst->cPrecalc->dCovarRatesAdjust_mid[index - 1];
    }

    end_index = CalibConst->lSigIndex[index];
    new_time = model->dPWTime[end_index];

    for (j = start_index; j < end_index + 1; j++) {
      if (j > start_index) {
        t1 = model->dPWTime[j - 1];
      } else {
        /* First part */
        t1 = last_time;
      }

      if (j == end_index || start_index == end_index) {
        /* Last part */
        t2 = new_time;
      } else {
        t2 = model->dPWTime[j];
      }

      dt = t2 - t1;

      if (dt > 0.0) {
        sig_dom = model->dSigmaDom[j];
        sig_for = model->dSigmaFor[j];
        sig_dom2 = sig_dom * sig_dom;
        sig_for2 = sig_for * sig_for;
        sig_domfor = model->dCorrDomFor[j] * sig_dom * sig_for;

        sig_fx = model->dSigmaFxUp[j];
        sig_fx2 = sig_fx * sig_fx;
        sig_domfx = model->dCorrDomFx[j] * sig_dom * sig_fx;
        sig_forfx = model->dCorrForFx[j] * sig_for * sig_fx;

        sig_fx_down = model->dSigmaFxDown[j];
        sig_fx2_down = sig_fx_down * sig_fx_down;
        sig_domfx_down = model->dCorrDomFxDown[j] * sig_dom * sig_fx_down;
        sig_forfx_down = model->dCorrForFxDown[j] * sig_for * sig_fx_down;

        sig_fx_mid = model->dSigmaFx[j];
        sig_fx2_mid = sig_fx_mid * sig_fx_mid;
        sig_domfx_mid = model->dCorrDomFxDown[j] * sig_dom * sig_fx_mid;
        sig_forfx_mid = model->dCorrForFxDown[j] * sig_for * sig_fx_mid;

        /* Expectation and Variance approximation for foreign */
        nb_points = max((int)(dt / CalibConst->cNumParams->dMinTime + 0.5), 1);
        min_dt = dt / nb_points;

        for (k = 0; k < nb_points; k++) {
          new_t = old_t + min_dt;
          temp_dom = exp(-dom_lam * (model->dTStar - new_t));
          temp_for = exp(-for_lam * (model->dTStar - new_t));

          /* update var_ffx up */
          var_ffx += sig_fx2 * min_dt;
          var_ffx +=
              -2.0 * sig_forfx * Etha_Func(for_lam, tstar, old_t, new_t) +
              2.0 * sig_domfx * Etha_Func(dom_lam, tstar, old_t, new_t);
          var_ffx +=
              sig_for2 * Psi_Func(for_lam, for_lam, tstar, old_t, new_t) +
              sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, old_t, new_t) -
              2.0 * sig_domfor *
                  Psi_Func(dom_lam, for_lam, tstar, old_t, new_t);

          covar_dom_ffx +=
              sig_domfx * Phi_Func(dom_lam, tstar, old_t, new_t) -
              sig_domfor * Gamma_Func(for_lam, dom_lam, tstar, old_t, new_t) +
              sig_dom2 * Gamma_Func(dom_lam, dom_lam, tstar, old_t, new_t);

          temp_bracket =
              temp_for * (sig_forfx - sig_for2 * (1.0 - temp_for) / for_lam +
                          sig_domfor * (1.0 - temp_dom) / dom_lam);

          new_integ = temp_bracket / (1.0 + 2.0 * model->dC0 * var_ffx);
          new_integ2 = covar_dom_ffx * new_integ;
          new_integ *= var_ffx;

          expect_for += (old_integ + new_integ) * min_dt / 2.0;
          covar_dom_for_adjust += (old_integ2 + new_integ2) * min_dt / 2.0;

          old_integ = new_integ;
          old_integ2 = new_integ2;

          /* update var_ffx down */
          var_ffx_down += sig_fx2_down * min_dt;
          var_ffx_down +=
              -2.0 * sig_forfx_down * Etha_Func(for_lam, tstar, old_t, new_t) +
              2.0 * sig_domfx_down * Etha_Func(dom_lam, tstar, old_t, new_t);
          var_ffx_down +=
              sig_for2 * Psi_Func(for_lam, for_lam, tstar, old_t, new_t) +
              sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, old_t, new_t) -
              2.0 * sig_domfor *
                  Psi_Func(dom_lam, for_lam, tstar, old_t, new_t);

          covar_dom_ffx_down +=
              sig_domfx_down * Phi_Func(dom_lam, tstar, old_t, new_t) -
              sig_domfor * Gamma_Func(for_lam, dom_lam, tstar, old_t, new_t) +
              sig_dom2 * Gamma_Func(dom_lam, dom_lam, tstar, old_t, new_t);

          temp_bracket_down =
              temp_for *
              (sig_forfx_down - sig_for2 * (1.0 - temp_for) / for_lam +
               sig_domfor * (1.0 - temp_dom) / dom_lam);

          new_integ_down =
              temp_bracket_down / (1.0 + 2.0 * model->dC0 * var_ffx_down);
          new_integ2_down = covar_dom_ffx_down * new_integ_down;
          new_integ_down *= var_ffx_down;

          expect_for_down += (old_integ_down + new_integ_down) * min_dt / 2.0;
          covar_dom_for_adjust_down +=
              (old_integ2_down + new_integ2_down) * min_dt / 2.0;

          old_integ_down = new_integ_down;
          old_integ2_down = new_integ2_down;

          /* update var_ffx mid */
          var_ffx_mid += sig_fx2_mid * min_dt;
          var_ffx_mid +=
              -2.0 * sig_forfx_mid * Etha_Func(for_lam, tstar, old_t, new_t) +
              2.0 * sig_domfx_mid * Etha_Func(dom_lam, tstar, old_t, new_t);
          var_ffx_mid +=
              sig_for2 * Psi_Func(for_lam, for_lam, tstar, old_t, new_t) +
              sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, old_t, new_t) -
              2.0 * sig_domfor *
                  Psi_Func(dom_lam, for_lam, tstar, old_t, new_t);

          covar_dom_ffx_mid +=
              sig_domfx_mid * Phi_Func(dom_lam, tstar, old_t, new_t) -
              sig_domfor * Gamma_Func(for_lam, dom_lam, tstar, old_t, new_t) +
              sig_dom2 * Gamma_Func(dom_lam, dom_lam, tstar, old_t, new_t);

          temp_bracket_mid =
              temp_for *
              (sig_forfx_mid - sig_for2 * (1.0 - temp_for) / for_lam +
               sig_domfor * (1.0 - temp_dom) / dom_lam);

          new_integ_mid =
              temp_bracket_mid / (1.0 + 2.0 * model->dC0 * var_ffx_mid);
          new_integ2_mid = covar_dom_ffx_mid * new_integ_mid;
          new_integ_mid *= var_ffx_mid;

          expect_for_mid += (old_integ_mid + new_integ_mid) * min_dt / 2.0;
          covar_dom_for_adjust_mid +=
              (old_integ2_mid + new_integ2_mid) * min_dt / 2.0;

          old_integ_mid = new_integ_mid;
          old_integ2_mid = new_integ2_mid;

          old_t = new_t;
        }

        covar_for_ffx +=
            sig_forfx * Phi_Func(for_lam, tstar, t1, t2) -
            sig_for2 * Gamma_Func(for_lam, for_lam, tstar, t1, t2) +
            sig_domfor * Gamma_Func(dom_lam, for_lam, tstar, t1, t2);
        covar_for_ffx_down +=
            sig_forfx_down * Phi_Func(for_lam, tstar, t1, t2) -
            sig_for2 * Gamma_Func(for_lam, for_lam, tstar, t1, t2) +
            sig_domfor * Gamma_Func(dom_lam, for_lam, tstar, t1, t2);
        covar_for_ffx_mid +=
            sig_forfx_mid * Phi_Func(for_lam, tstar, t1, t2) -
            sig_for2 * Gamma_Func(for_lam, for_lam, tstar, t1, t2) +
            sig_domfor * Gamma_Func(dom_lam, for_lam, tstar, t1, t2);
      }
    }

    for_fwd = -covar_for_ffx + 2.0 * model->dC0 * expect_for;
    for_fwd_down = -covar_for_ffx_down + 2.0 * model->dC0 * expect_for_down;
    for_fwd_mid = -covar_for_ffx_mid + 2.0 * model->dC0 * expect_for_mid;

    /* Update the saved coefficients */
    if (CalibConst->cPrecalc) {

      CalibConst->cPrecalc->dVarFFx[CalibConst->cPrecalc->iNbDone] = var_ffx;
      CalibConst->cPrecalc->dCovarFFxDom[CalibConst->cPrecalc->iNbDone] =
          covar_dom_ffx;
      CalibConst->cPrecalc->dCovarFFxFor[CalibConst->cPrecalc->iNbDone] =
          -for_fwd;
      CalibConst->cPrecalc->dCovarFFxFor3D[CalibConst->cPrecalc->iNbDone] =
          covar_for_ffx;
      CalibConst->cPrecalc->dExpectFor[CalibConst->cPrecalc->iNbDone] =
          expect_for;
      CalibConst->cPrecalc->dIntegral[CalibConst->cPrecalc->iNbDone] =
          new_integ;
      CalibConst->cPrecalc->dIntegral2[CalibConst->cPrecalc->iNbDone] =
          new_integ2;
      CalibConst->cPrecalc->dCovarRatesAdjust[CalibConst->cPrecalc->iNbDone] =
          covar_dom_for_adjust;

      CalibConst->cPrecalc->dVarFFx_down[CalibConst->cPrecalc->iNbDone] =
          var_ffx_down;
      CalibConst->cPrecalc->dCovarFFxDom_down[CalibConst->cPrecalc->iNbDone] =
          covar_dom_ffx_down;
      CalibConst->cPrecalc->dCovarFFxFor_down[CalibConst->cPrecalc->iNbDone] =
          -for_fwd_down;
      CalibConst->cPrecalc->dCovarFFxFor3D_down[CalibConst->cPrecalc->iNbDone] =
          covar_for_ffx_down;
      CalibConst->cPrecalc->dExpectFor_down[CalibConst->cPrecalc->iNbDone] =
          expect_for_down;
      CalibConst->cPrecalc->dIntegral_down[CalibConst->cPrecalc->iNbDone] =
          new_integ_down;
      CalibConst->cPrecalc->dIntegral2_down[CalibConst->cPrecalc->iNbDone] =
          new_integ2_down;
      CalibConst->cPrecalc
          ->dCovarRatesAdjust_down[CalibConst->cPrecalc->iNbDone] =
          covar_dom_for_adjust_down;

      CalibConst->cPrecalc->dVarFFx_mid[CalibConst->cPrecalc->iNbDone] =
          var_ffx_mid;
      CalibConst->cPrecalc->dCovarFFxDom_mid[CalibConst->cPrecalc->iNbDone] =
          covar_dom_ffx_mid;
      CalibConst->cPrecalc->dCovarFFxFor_mid[CalibConst->cPrecalc->iNbDone] =
          -for_fwd_mid;
      CalibConst->cPrecalc->dCovarFFxFor3D_mid[CalibConst->cPrecalc->iNbDone] =
          covar_for_ffx_mid;
      CalibConst->cPrecalc->dExpectFor_mid[CalibConst->cPrecalc->iNbDone] =
          expect_for_mid;
      CalibConst->cPrecalc->dIntegral_mid[CalibConst->cPrecalc->iNbDone] =
          new_integ_mid;
      CalibConst->cPrecalc->dIntegral2_mid[CalibConst->cPrecalc->iNbDone] =
          new_integ2_mid;
      CalibConst->cPrecalc
          ->dCovarRatesAdjust_mid[CalibConst->cPrecalc->iNbDone] =
          covar_dom_for_adjust_mid;

      /* Add the adjustments to the FORVAR and RATESCOV*/
      CalibConst->cPrecalc->dVarRatesFor[CalibConst->cPrecalc->iNbDone] =
          CalibConst->cPrecalc->dVarRatesFor3D[CalibConst->cPrecalc->iNbDone] -
          2.0 * model->dC0 * for_fwd * for_fwd;
      CalibConst->cPrecalc->dCovarRates[CalibConst->cPrecalc->iNbDone] =
          CalibConst->cPrecalc->dCovarRates3D[CalibConst->cPrecalc->iNbDone] -
          2.0 * model->dC0 * covar_dom_for_adjust;

      CalibConst->cPrecalc->dVarRatesFor_down[CalibConst->cPrecalc->iNbDone] =
          CalibConst->cPrecalc->dVarRatesFor3D[CalibConst->cPrecalc->iNbDone] -
          2.0 * model->dC0 * for_fwd_down * for_fwd_down;
      CalibConst->cPrecalc->dCovarRates_down[CalibConst->cPrecalc->iNbDone] =
          CalibConst->cPrecalc->dCovarRates3D[CalibConst->cPrecalc->iNbDone] -
          2.0 * model->dC0 * covar_dom_for_adjust_down;

      CalibConst->cPrecalc->dVarRatesFor_mid[CalibConst->cPrecalc->iNbDone] =
          CalibConst->cPrecalc->dVarRatesFor3D[CalibConst->cPrecalc->iNbDone] -
          2.0 * model->dC0 * for_fwd_mid * for_fwd_mid;
      CalibConst->cPrecalc->dCovarRates_mid[CalibConst->cPrecalc->iNbDone] =
          CalibConst->cPrecalc->dCovarRates3D[CalibConst->cPrecalc->iNbDone] -
          2.0 * model->dC0 * covar_dom_for_adjust_mid;
    }
  }
  return err;
}

Err FxBetaDLM_Update_InstPrecalc_FromMoment(FxBetaDLM_FxOptInst *Inst,
                                            FxBetaDLM_model *model,
                                            FxBetaDLM_OptNumerParams *NumParams,
                                            FxBetaDLM_InstPrecalc *InstConst) {
  Err err = NULL;
  double k1, k2, k3;
  double var_rates;
  double fixed_adjustment, lin_adjustment;
  double temp, temp_fwd, coef_b, coef_c;
  double expectX;
  double alpha;
  double varY, varZ, covFFxY, covFFxZ;

  /* Pricing With a first conditioning on the quadratic variable*/
  if (NumParams->iMethod == 1) {
    expectX = NumParams->dX0 - InstConst->dCovarFFxDom * InstConst->dBetaDom;

    k1 = (InstConst->dVarFFx * InstConst->dVarRatesDom -
          InstConst->dCovarFFxDom * InstConst->dCovarFFxDom) /
         InstConst->dVarFFx;
    k2 = (InstConst->dVarFFx * InstConst->dVarRatesFor -
          InstConst->dCovarFFxFor * InstConst->dCovarFFxFor) /
         InstConst->dVarFFx;
    k3 = (InstConst->dVarFFx * InstConst->dCovarRates -
          InstConst->dCovarFFxDom * InstConst->dCovarFFxFor) /
         InstConst->dVarFFx;

    varY =
        InstConst->dBetaDom * InstConst->dBetaDom * InstConst->dVarRatesDom +
        InstConst->dBetaFor * InstConst->dBetaFor * InstConst->dVarRatesFor -
        2 * InstConst->dBetaDom * InstConst->dBetaFor * InstConst->dCovarRates;

    covFFxY = InstConst->dBetaDom * InstConst->dCovarFFxDom -
              InstConst->dBetaFor * InstConst->dCovarFFxFor;

    var_rates = k1 * InstConst->dBetaDom * InstConst->dBetaDom +
                k2 * InstConst->dBetaFor * InstConst->dBetaFor -
                2.0 * k3 * InstConst->dBetaDom * InstConst->dBetaFor;

    temp = 1.0 + 2.0 * model->dC0 * InstConst->dTotalVarFFx;

    if (temp < 1.0E-10) {
      err = "invalid model parameters in FxBetaDLM_FxOption";
      return err;
    }

    coef_c = model->dC0 / temp;
    coef_b = model->dB0 / temp + 2.0 * expectX * coef_c;
    temp_fwd = 1.0 - 2.0 * coef_c * InstConst->dVarFFx;

    lin_adjustment = (InstConst->dBetaDom * InstConst->dCovarFFxDom -
                      InstConst->dBetaFor * InstConst->dCovarFFxFor) /
                     InstConst->dVarFFx;

    fixed_adjustment =
        Inst->dLnFwdFx + 0.5 * log(temp_fwd) -
        0.5 * InstConst->dVarFFx / temp_fwd * pow(coef_b + lin_adjustment, 2);

    lin_adjustment = (lin_adjustment + coef_b) * InstConst->dStdFFX;

    /* Save the information needed for the Forward into the structure*/
    if (InstConst->dVarFFx == 0) {
      InstConst->dFwdQuad2 = 1.0;
      InstConst->dFwdAlpha = 0.0;
      InstConst->dFwdBeta = exp(0.5 * varY * varY);
    } else {
      InstConst->dFwdQuad2 = 1 - 2.0 * InstConst->dVarFFx * coef_c;
      if (InstConst->dFwdQuad2 <= 0.0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Infinite variance "
              "process";
        return err;
      }
      InstConst->dFwdQuad2 = InstConst->dVarFFx / InstConst->dFwdQuad2;
      InstConst->dFwdAlpha = covFFxY / InstConst->dVarFFx;
      InstConst->dFwdBeta =
          exp(0.5 * (varY * InstConst->dVarFFx - covFFxY * covFFxY) /
              InstConst->dVarFFx) /
          sqrt(InstConst->dFwdQuad2);
    }

    /* Save results in the structure */
    InstConst->dStdRates = sqrt(var_rates);
    InstConst->dConstCoef = fixed_adjustment;
    InstConst->dLinCoef = lin_adjustment;
    InstConst->dQuadCoef = coef_c * InstConst->dVarFFx;
  } else /* Pricing With a first conditioning on the Linear varible*/
      if (NumParams->iMethod == 2) {
    expectX = NumParams->dX0 - InstConst->dCovarFFxDom * InstConst->dBetaDom;
    temp = 1.0 + 2.0 * model->dC0 * InstConst->dTotalVarFFx;
    coef_c = model->dC0 / temp;
    coef_b = model->dB0 / temp + 2.0 * coef_c * expectX;
    temp_fwd = 1.0 - 2.0 * coef_c * InstConst->dVarFFx;

    varY =
        InstConst->dBetaDom * InstConst->dBetaDom * InstConst->dVarRatesDom +
        InstConst->dBetaFor * InstConst->dBetaFor * InstConst->dVarRatesFor -
        2 * InstConst->dBetaDom * InstConst->dBetaFor * InstConst->dCovarRates;

    covFFxY = InstConst->dBetaDom * InstConst->dCovarFFxDom -
              InstConst->dBetaFor * InstConst->dCovarFFxFor;

    if (fabs(InstConst->dVarFFx) > 1e-10) {
      fixed_adjustment = Inst->dLnFwdFx + 0.5 * log(temp_fwd) -
                         0.5 * (InstConst->dVarFFx * varY - covFFxY * covFFxY) /
                             InstConst->dVarFFx -
                         0.5 * pow(coef_b + covFFxY / InstConst->dVarFFx, 2) *
                             InstConst->dVarFFx / temp_fwd;
    } else {
      fixed_adjustment = Inst->dLnFwdFx + 0.5 * log(temp_fwd);
    }

    /* Save the information needed for the Forward into the structure*/
    if (InstConst->dVarFFx == 0) {
      InstConst->dFwdQuad2 = 1.0;
      InstConst->dFwdAlpha = 0.0;
      InstConst->dFwdBeta = exp(0.5 * varY * varY);
    } else {
      InstConst->dFwdQuad2 = 1 - 2.0 * InstConst->dVarFFx * coef_c;
      if (InstConst->dFwdQuad2 <= 0.0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Infinite variance "
              "process";
        return err;
      }
      InstConst->dFwdQuad2 = InstConst->dVarFFx / InstConst->dFwdQuad2;
      InstConst->dFwdAlpha = covFFxY / InstConst->dVarFFx;
      InstConst->dFwdBeta =
          1 / sqrt(InstConst->dFwdQuad2) *
          exp(0.5 * (varY * InstConst->dVarFFx - covFFxY * covFFxY) /
              InstConst->dVarFFx);
    }

    /* Save results in the structure */
    InstConst->dQuadConst = fixed_adjustment;
    InstConst->dQuadLin = coef_b;
    InstConst->dQuadQuad = coef_c;
    if (varY == 0) {
      InstConst->dQuadMean = 0;
      InstConst->dQuadVar = InstConst->dVarFFx;
      if (InstConst->dQuadVar < 0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Negative Variance";
      }
    } else {
      InstConst->dQuadMean = covFFxY / varY;
      InstConst->dQuadVar =
          (varY * InstConst->dVarFFx - covFFxY * covFFxY) / varY;
      if (InstConst->dQuadVar < 0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Negative Variance";
      }
    }
    InstConst->dQuadStd = sqrt(InstConst->dQuadVar);
    InstConst->dQuadStdY = sqrt(varY);

    if (InstConst->iIsCPD) {
      InstConst->coef_b_const =
          model->dB0 / temp -
          2.0 * coef_c * InstConst->dCovarFFxDom * InstConst->dBetaDom;
      InstConst->coef_b_lin = 2.0 * coef_c;

      InstConst->fixed_adj_const =
          0.5 * log(temp_fwd) -
          0.5 * (InstConst->dVarFFx * varY - covFFxY * covFFxY) /
              InstConst->dVarFFx;

      InstConst->fixed_adj_lin = -0.5 * InstConst->dVarFFx / temp_fwd;
      InstConst->fixed_adj_lin_const = covFFxY / InstConst->dVarFFx;
    }

    if (fabs(model->dAlpha) > TINY) {
      expectX =
          NumParams->dX0 - InstConst->dCovarFFxDom_down * InstConst->dBetaDom;
      temp = 1.0 + 2.0 * model->dC0 * InstConst->dTotalVarFFx_down;
      coef_c = model->dC0 / temp;
      coef_b = model->dB0 / temp + 2.0 * coef_c * expectX;
      temp_fwd = 1.0 - 2.0 * coef_c * InstConst->dVarFFx_down;

      varY =
          InstConst->dBetaDom * InstConst->dBetaDom * InstConst->dVarRatesDom +
          InstConst->dBetaFor * InstConst->dBetaFor *
              InstConst->dVarRatesFor_down -
          2 * InstConst->dBetaDom * InstConst->dBetaFor *
              InstConst->dCovarRates_down;

      covFFxY = InstConst->dBetaDom * InstConst->dCovarFFxDom_down -
                InstConst->dBetaFor * InstConst->dCovarFFxFor_down;

      if (fabs(InstConst->dVarFFx_down) > 1e-10) {
        fixed_adjustment =
            Inst->dLnFwdFx + 0.5 * log(temp_fwd) -
            0.5 * (InstConst->dVarFFx_down * varY - covFFxY * covFFxY) /
                InstConst->dVarFFx_down -
            0.5 * pow(coef_b + covFFxY / InstConst->dVarFFx_down, 2) *
                InstConst->dVarFFx_down / temp_fwd;
      } else {
        fixed_adjustment = Inst->dLnFwdFx + 0.5 * log(temp_fwd);
      }

      /* Save the information needed for the Forward into the structure*/
      if (InstConst->dVarFFx_down == 0) {
        InstConst->dFwdQuad2_down = 1.0;
        InstConst->dFwdAlpha_down = 0.0;
        InstConst->dFwdBeta_down = exp(0.5 * varY * varY);
      } else {
        InstConst->dFwdQuad2_down = 1 - 2.0 * InstConst->dVarFFx_down * coef_c;
        if (InstConst->dFwdQuad2_down <= 0.0) {
          err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Infinite variance "
                "process";
          return err;
        }
        InstConst->dFwdQuad2_down =
            InstConst->dVarFFx_down / InstConst->dFwdQuad2_down;
        InstConst->dFwdAlpha_down = covFFxY / InstConst->dVarFFx_down;
        InstConst->dFwdBeta_down =
            1 / sqrt(InstConst->dFwdQuad2_down) *
            exp(0.5 * (varY * InstConst->dVarFFx_down - covFFxY * covFFxY) /
                InstConst->dVarFFx_down);
      }

      /* Save results in the structure */
      InstConst->dQuadConst_down = fixed_adjustment;
      InstConst->dQuadLin_down = coef_b;
      InstConst->dQuadQuad_down = coef_c;
      if (varY == 0) {
        InstConst->dQuadMean_down = 0;
        InstConst->dQuadVar_down = InstConst->dVarFFx_down;
      } else {
        InstConst->dQuadMean_down = covFFxY / varY;
        InstConst->dQuadVar_down =
            (varY * InstConst->dVarFFx_down - covFFxY * covFFxY) / varY;
      }
      InstConst->dQuadStd_down = sqrt(InstConst->dQuadVar_down);
      InstConst->dQuadStdY_down = sqrt(varY);

      if (InstConst->iIsCPD) {
        InstConst->coef_b_const_down =
            model->dB0 / temp -
            2.0 * coef_c * InstConst->dCovarFFxDom_down * InstConst->dBetaDom;
        InstConst->coef_b_lin_down = 2.0 * coef_c;

        InstConst->fixed_adj_const_down =
            0.5 * log(temp_fwd) -
            0.5 * (InstConst->dVarFFx * varY - covFFxY * covFFxY) /
                InstConst->dVarFFx_down;

        InstConst->fixed_adj_lin_down =
            -0.5 * InstConst->dVarFFx_down / temp_fwd;
        InstConst->fixed_adj_lin_const_down = covFFxY / InstConst->dVarFFx_down;
      }
    }
  } else /* Pricing With a first conditioning on a minimum variance Linear
            varible*/
      if (NumParams->iMethod == 3) {
    varY =
        InstConst->dBetaDom * InstConst->dBetaDom * InstConst->dVarRatesDom +
        InstConst->dBetaFor * InstConst->dBetaFor * InstConst->dVarRatesFor -
        2 * InstConst->dBetaDom * InstConst->dBetaFor * InstConst->dCovarRates;
    covFFxY = InstConst->dBetaDom * InstConst->dCovarFFxDom -
              InstConst->dBetaFor * InstConst->dCovarFFxFor;

    alpha = -covFFxY / InstConst->dVarFFx;
    varZ = varY - alpha * alpha * InstConst->dVarFFx;
    covFFxZ = alpha * InstConst->dVarFFx + covFFxY;

    expectX = -InstConst->dCovarFFxDom * InstConst->dBetaDom;
    temp = 1.0 + 2.0 * model->dC0 * InstConst->dTotalVarFFx;
    coef_c = model->dC0 / temp;
    coef_b = model->dB0 / temp + 2.0 * coef_c * expectX - alpha;
    temp_fwd = 1.0 - 2.0 * coef_c * InstConst->dVarFFx;

    if (fabs(InstConst->dVarFFx) > 1e-10) {
      fixed_adjustment = Inst->dLnFwdFx + 0.5 * log(temp_fwd) -
                         0.5 * (InstConst->dVarFFx * varZ - covFFxZ * covFFxZ) /
                             InstConst->dVarFFx -
                         0.5 * pow(coef_b + covFFxZ / InstConst->dVarFFx, 2) *
                             InstConst->dVarFFx / temp_fwd;
    } else {
      fixed_adjustment = Inst->dLnFwdFx + 0.5 * log(temp_fwd);
    }

    /* Save results in the structure */
    InstConst->dQuadConst = fixed_adjustment;
    InstConst->dQuadLin = coef_b;
    InstConst->dQuadQuad = coef_c;
    if (varY == 0) {
      InstConst->dQuadMean = 0;
      InstConst->dQuadVar = InstConst->dVarFFx;
      if (InstConst->dQuadVar < 0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Negative Variance";
      }
    } else {
      InstConst->dQuadMean = covFFxZ / varZ;
      InstConst->dQuadVar =
          (varZ * InstConst->dVarFFx - covFFxZ * covFFxZ) / varZ;
      if (InstConst->dQuadVar < 0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Negative Variance";
      }
    }
    InstConst->dQuadStd = sqrt(InstConst->dQuadVar);
    InstConst->dQuadStdY = sqrt(varZ);
  } else /* Pricing With a first conditioning on a minimum variance Linear
            varible*/
      if (NumParams->iMethod == 4) {
    varY =
        InstConst->dBetaDom * InstConst->dBetaDom * InstConst->dVarRatesDom +
        InstConst->dBetaFor * InstConst->dBetaFor * InstConst->dVarRatesFor -
        2 * InstConst->dBetaDom * InstConst->dBetaFor * InstConst->dCovarRates;
    covFFxY = InstConst->dBetaDom * InstConst->dCovarFFxDom -
              InstConst->dBetaFor * InstConst->dCovarFFxFor;

    alpha = NumParams->dAlpha;
    varZ = varY + alpha * alpha * InstConst->dVarFFx + 2.0 * alpha * covFFxY;
    covFFxZ = alpha * InstConst->dVarFFx + covFFxY;

    expectX = -InstConst->dCovarFFxDom * InstConst->dBetaDom;
    temp = 1.0 + 2.0 * model->dC0 * InstConst->dTotalVarFFx;
    coef_c = model->dC0 / temp;
    coef_b = model->dB0 / temp + 2.0 * coef_c * expectX - alpha;
    temp_fwd = 1.0 - 2.0 * coef_c * InstConst->dVarFFx;

    if (fabs(InstConst->dVarFFx) > 1e-10) {
      fixed_adjustment = Inst->dLnFwdFx + 0.5 * log(temp_fwd) -
                         0.5 * (InstConst->dVarFFx * varZ - covFFxZ * covFFxZ) /
                             InstConst->dVarFFx -
                         0.5 * pow(coef_b + covFFxZ / InstConst->dVarFFx, 2) *
                             InstConst->dVarFFx / temp_fwd;
    } else {
      fixed_adjustment = Inst->dLnFwdFx + 0.5 * log(temp_fwd);
    }

    /* Save results in the structure */
    InstConst->dQuadConst = fixed_adjustment;
    InstConst->dQuadLin = coef_b;
    InstConst->dQuadQuad = coef_c;
    if (varZ == 0) {
      InstConst->dQuadMean = 0;
      InstConst->dQuadVar = InstConst->dVarFFx;
      if (InstConst->dQuadVar < 0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Negative Variance";
      }
    } else {
      InstConst->dQuadMean = covFFxZ / varZ;
      InstConst->dQuadVar =
          (varZ * InstConst->dVarFFx - covFFxZ * covFFxZ) / varZ;
      if (InstConst->dQuadVar < 0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Negative Variance";
      }
    }
    InstConst->dQuadStd = sqrt(InstConst->dQuadVar);
    InstConst->dQuadStdY = sqrt(varZ);
  } else {
    err =
        "FxBetaDLM_Update_InstPrecalc_FromMoment: Parameter Method is invalid";
    return err;
  }
  return err;
}

Err FxBetaDLM_Update_InstPrecalc_FromX0(double X0, FxBetaDLM_FxOptInst *Inst,
                                        FxBetaDLM_model *model,
                                        FxBetaDLM_OptNumerParams *NumParams,
                                        FxBetaDLM_InstPrecalc *InstConst) {
  Err err = NULL;
  double k1, k2, k3;
  double var_rates;
  double fixed_adjustment, lin_adjustment;
  double temp, temp_fwd, coef_b, coef_c;
  double expectX;
  double alpha;
  double varY, varZ, covFFxY, covFFxZ;

  /* Pricing With a first conditioning on the quadratic variable*/
  if (NumParams->iMethod == 1) {
    expectX = X0 - InstConst->dCovarFFxDom * InstConst->dBetaDom;

    k1 = (InstConst->dVarFFx * InstConst->dVarRatesDom -
          InstConst->dCovarFFxDom * InstConst->dCovarFFxDom) /
         InstConst->dVarFFx;
    k2 = (InstConst->dVarFFx * InstConst->dVarRatesFor -
          InstConst->dCovarFFxFor * InstConst->dCovarFFxFor) /
         InstConst->dVarFFx;
    k3 = (InstConst->dVarFFx * InstConst->dCovarRates -
          InstConst->dCovarFFxDom * InstConst->dCovarFFxFor) /
         InstConst->dVarFFx;

    var_rates = k1 * InstConst->dBetaDom * InstConst->dBetaDom +
                k2 * InstConst->dBetaFor * InstConst->dBetaFor -
                2.0 * k3 * InstConst->dBetaDom * InstConst->dBetaFor;

    temp = 1.0 + 2.0 * model->dC0 * InstConst->dTotalVarFFx;

    if (temp < 1.0E-10) {
      err = "invalid model parameters in FxBetaDLM_FxOption";
      return err;
    }

    coef_c = model->dC0 / temp;
    coef_b = model->dB0 / temp + 2.0 * expectX * coef_c;
    temp_fwd = 1.0 - 2.0 * coef_c * InstConst->dVarFFx;

    lin_adjustment = (InstConst->dBetaDom * InstConst->dCovarFFxDom -
                      InstConst->dBetaFor * InstConst->dCovarFFxFor) /
                     InstConst->dVarFFx;

    fixed_adjustment =
        Inst->dLnFwdFx + 0.5 * log(temp_fwd) -
        0.5 * InstConst->dVarFFx / temp_fwd * pow(coef_b + lin_adjustment, 2);

    lin_adjustment = (lin_adjustment + coef_b) * InstConst->dStdFFX;

    /* Save results in the structure */
    InstConst->dStdRates = sqrt(var_rates);
    InstConst->dConstCoef = fixed_adjustment;
    InstConst->dLinCoef = lin_adjustment;
    InstConst->dQuadCoef = coef_c * InstConst->dVarFFx;
  } else /* Pricing With a first conditioning on the Linear varible*/
      if (NumParams->iMethod == 2) {
    InstConst->dQuadLin = InstConst->coef_b_const + InstConst->coef_b_lin * X0;
    InstConst->dQuadConst =
        Inst->dLnFwdFx + InstConst->fixed_adj_const +
        InstConst->fixed_adj_lin *
            pow(InstConst->dQuadLin + InstConst->fixed_adj_lin_const, 2);
  } else /* Pricing With a first conditioning on a minimum variance Linear
            varible*/
      if (NumParams->iMethod == 3) {
    varY =
        InstConst->dBetaDom * InstConst->dBetaDom * InstConst->dVarRatesDom +
        InstConst->dBetaFor * InstConst->dBetaFor * InstConst->dVarRatesFor -
        2 * InstConst->dBetaDom * InstConst->dBetaFor * InstConst->dCovarRates;
    covFFxY = InstConst->dBetaDom * InstConst->dCovarFFxDom -
              InstConst->dBetaFor * InstConst->dCovarFFxFor;

    alpha = -covFFxY / InstConst->dVarFFx;
    varZ = varY - alpha * alpha * InstConst->dVarFFx;
    covFFxZ = alpha * InstConst->dVarFFx + covFFxY;

    expectX = -InstConst->dCovarFFxDom * InstConst->dBetaDom;
    temp = 1.0 + 2.0 * model->dC0 * InstConst->dTotalVarFFx;
    coef_c = model->dC0 / temp;
    coef_b = model->dB0 / temp + 2.0 * coef_c * expectX - alpha;

    if (fabs(InstConst->dVarFFx) > 1e-10) {
      fixed_adjustment = Inst->dLnFwdFx - 0.5 * log(temp) -
                         0.5 * (InstConst->dVarFFx * varZ - covFFxZ * covFFxZ) /
                             InstConst->dVarFFx -
                         0.5 * pow(coef_b + covFFxZ / InstConst->dVarFFx, 2) *
                             InstConst->dVarFFx * temp;
    } else {
      fixed_adjustment = Inst->dLnFwdFx - 0.5 * log(temp);
    }

    /* Save results in the structure */
    InstConst->dQuadConst = fixed_adjustment;
    InstConst->dQuadLin = coef_b;
    InstConst->dQuadQuad = coef_c;
    if (varY == 0) {
      InstConst->dQuadMean = 0;
      InstConst->dQuadVar = InstConst->dVarFFx;
      if (InstConst->dQuadVar < 0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Negative Variance";
      }
    } else {
      InstConst->dQuadMean = covFFxZ / varZ;
      InstConst->dQuadVar =
          (varZ * InstConst->dVarFFx - covFFxZ * covFFxZ) / varZ;
      if (InstConst->dQuadVar < 0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Negative Variance";
      }
    }
    InstConst->dQuadStdY = sqrt(varZ);
  } else /* Pricing With a first conditioning on a minimum variance Linear
            varible*/
      if (NumParams->iMethod == 4) {
    varY =
        InstConst->dBetaDom * InstConst->dBetaDom * InstConst->dVarRatesDom +
        InstConst->dBetaFor * InstConst->dBetaFor * InstConst->dVarRatesFor -
        2 * InstConst->dBetaDom * InstConst->dBetaFor * InstConst->dCovarRates;
    covFFxY = InstConst->dBetaDom * InstConst->dCovarFFxDom -
              InstConst->dBetaFor * InstConst->dCovarFFxFor;

    alpha = NumParams->dAlpha;
    varZ = varY + alpha * alpha * InstConst->dVarFFx + 2.0 * alpha * covFFxY;
    covFFxZ = alpha * InstConst->dVarFFx + covFFxY;

    expectX = -InstConst->dCovarFFxDom * InstConst->dBetaDom;
    temp = 1.0 + 2.0 * model->dC0 * InstConst->dTotalVarFFx;
    coef_c = model->dC0 / temp;
    coef_b = model->dB0 / temp + 2.0 * coef_c * expectX - alpha;

    if (fabs(InstConst->dVarFFx) > 1e-10) {
      fixed_adjustment = Inst->dLnFwdFx - 0.5 * log(temp) -
                         0.5 * (InstConst->dVarFFx * varZ - covFFxZ * covFFxZ) /
                             InstConst->dVarFFx -
                         0.5 * pow(coef_b + covFFxZ / InstConst->dVarFFx, 2) *
                             InstConst->dVarFFx * temp;
    } else {
      fixed_adjustment = Inst->dLnFwdFx - 0.5 * log(temp);
    }

    /* Save results in the structure */
    InstConst->dQuadConst = fixed_adjustment;
    InstConst->dQuadLin = coef_b;
    InstConst->dQuadQuad = coef_c;
    if (varZ == 0) {
      InstConst->dQuadMean = 0;
      InstConst->dQuadVar = InstConst->dVarFFx;
      if (InstConst->dQuadVar < 0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Negative Variance";
      }
    } else {
      InstConst->dQuadMean = covFFxZ / varZ;
      InstConst->dQuadVar =
          (varZ * InstConst->dVarFFx - covFFxZ * covFFxZ) / varZ;
      if (InstConst->dQuadVar < 0) {
        err = "FxBetaDLM_Update_InstPrecalc_FromMoment: Negative Variance";
      }
    }
    InstConst->dQuadStdY = sqrt(varZ);
  } else {
    err =
        "FxBetaDLM_Update_InstPrecalc_FromMoment: Parameter Method is invalid";
    return err;
  }
  return err;
}

Err FxBetaDLM_Update_InstPrecalc(FxBetaDLM_FxOptInst *Inst,
                                 FxBetaDLM_model *model,
                                 FxBetaDLM_CalibConst *CalibConst,
                                 FxBetaDLM_InstPrecalc *InstConst) {
  Err err = NULL;

  InstConst->dVarRatesDom =
      CalibConst->cPrecalc->dVarRatesDom[Inst->iInstIndex];
  InstConst->dVarRatesFor =
      CalibConst->cPrecalc->dVarRatesFor[Inst->iInstIndex];
  InstConst->dCovarRates = CalibConst->cPrecalc->dCovarRates[Inst->iInstIndex];

  InstConst->dVarFFx = CalibConst->cPrecalc->dVarFFx[Inst->iInstIndex];
  InstConst->dTotalVarFFx = InstConst->dVarFFx;
  InstConst->dStdFFX = sqrt(InstConst->dVarFFx);

  InstConst->dCovarFFxDom =
      CalibConst->cPrecalc->dCovarFFxDom[Inst->iInstIndex];
  InstConst->dCovarFFxFor =
      CalibConst->cPrecalc->dCovarFFxFor[Inst->iInstIndex];

  if (fabs(model->dAlpha) > TINY) {
    InstConst->dVarRatesFor_down =
        CalibConst->cPrecalc->dVarRatesFor_down[Inst->iInstIndex];
    InstConst->dCovarRates_down =
        CalibConst->cPrecalc->dCovarRates_down[Inst->iInstIndex];

    InstConst->dVarFFx_down =
        CalibConst->cPrecalc->dVarFFx_down[Inst->iInstIndex];
    InstConst->dTotalVarFFx_down = InstConst->dVarFFx_down;
    InstConst->dStdFFX_down = sqrt(InstConst->dVarFFx_down);

    InstConst->dCovarFFxDom_down =
        CalibConst->cPrecalc->dCovarFFxDom_down[Inst->iInstIndex];
    InstConst->dCovarFFxFor_down =
        CalibConst->cPrecalc->dCovarFFxFor_down[Inst->iInstIndex];
  }

  err = FxBetaDLM_Update_InstPrecalc_FromMoment(
      Inst, model, CalibConst->cNumParams, InstConst);

  return err;
}

/* All the calibration functions */

Err FxBetaDLM_GetTargetVol(void *Inst_, void *GlobalConst, void *Model,
                           CALIBGEN_PARAMS CalibConsts, double *target) {
  FxBetaDLM_FxOptInst *Inst = (FxBetaDLM_FxOptInst *)(Inst_);

  *target = Inst->dPrice[0];

  return NULL;
}

Err FxBetaDLM_GetFirstGuessVol(void *Model_, void *GlobalConst, int index_param,
                               double target, double *vol1) {
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);

  /* get the new first guess */
  *vol1 = model->dSigmaFx[CalibConst->lSigIndex[index_param]];

  if (index_param > 0) {
    *vol1 += model->dSigmaFx[CalibConst->lSigIndex[index_param - 1]] -
             CalibConst->first_guess;

    if (*vol1 < 0.001) {
      *vol1 = 0.001;
    }
  }

  /* update the first guess */
  CalibConst->first_guess = model->dSigmaFx[CalibConst->lSigIndex[index_param]];

  return NULL;
}

Err FxBetaDLM_GetSecondGuessVol(void *Model_, void *GlobalConst, int vol_index,
                                double vol1, double price1, double target,
                                double *vol2) {
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);
  double t1, t2, dt;
  double cum_var_temp, cum_var;

  if (vol_index == 0) {
    t1 = 0.0;
    cum_var_temp = 0.0;
    cum_var = 0.0;
  } else {
    t1 = model->dPWTime[CalibConst->lSigIndex[vol_index - 1]];
    cum_var = CalibConst->cPrecalc->dVarX[vol_index - 1];
  }

  t2 = model->dPWTime[CalibConst->lSigIndex[vol_index]];

  dt = t2 - t1;

  cum_var_temp = cum_var + vol1 * vol1 * dt;

  *vol2 = (target * target / price1 / price1 * cum_var_temp - cum_var) / dt;

  if (*vol2 < 1.0E-08) {
    *vol2 = 1.0E-08;
  } else {
    *vol2 = sqrt(*vol2);
  }

  return NULL;
}

Err FxBetaDLM_GetLimitAndLastVol(void *Model_, CALIBGEN_PARAMS CalibParams,
                                 void *GlobalConst,

                                 int vol_index, double *last_vol,
                                 double *limit_down, double *limit_up) {
  Err err = NULL;

  double t1, t2, dt;
  int last_index, new_index;

  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);

  if (vol_index == 0) {
    *last_vol = 0.00001;
    *limit_down = 0.00001;
    *limit_up = 1000.0;
    return err;
  }

  last_index = CalibConst->lSigIndex[vol_index - 1];
  new_index = CalibConst->lSigIndex[vol_index];

  t1 = model->dPWTime[last_index];
  t2 = model->dPWTime[new_index];

  dt = t2 - t1;

  *last_vol = model->dSigmaFx[last_index];
  *limit_down = sqrt(CalibConst->cPrecalc->dVarX[vol_index - 1] *
                     CalibConst->dMinFact / t1);
  *limit_up = sqrt(CalibConst->cPrecalc->dVarX[vol_index - 1] /
                   CalibConst->dMaxFact / t1);

  return err;
}

Err FxBetaDLM_BumpVol(void *Model_, void *GlobalConst, int vol_index,
                      double vol) {
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);

  int i;
  int last_index, new_index, var_index;
  double varFFx, varFFxUp, varFFxDown;
  double t1, t2;
  Err err = NULL;

  if (vol_index == 0) {
    last_index = 0;
  } else {
    last_index = CalibConst->lSigIndex[vol_index - 1] + 1;
  }

  new_index = CalibConst->lSigIndex[vol_index];

  if (vol_index == 0) {
    var_index = 0;
    varFFx = 0.0;
    varFFxUp = 0.0;
    varFFxDown = 0.0;
    t1 = 0.0;
  } else {
    var_index = vol_index - 1;
    if (fabs(model->dAlpha) < TINY)
      varFFx = CalibConst->cPrecalc->dVarFFx[var_index];
    else {
      varFFx = CalibConst->cPrecalc->dVarFFx_mid[var_index];
      varFFxUp = CalibConst->cPrecalc->dVarFFx[var_index];
      varFFxDown = CalibConst->cPrecalc->dVarFFx_down[var_index];
    }
    t1 = model->dPWTime[CalibConst->lSigIndex[vol_index - 1]];
  }

  for (i = last_index; i <= new_index; i++) {
    model->dSigmaFx[i] = vol;
    if (fabs(model->dLambda) < TINY) {
      model->dSigmaFxUp[i] =
          vol * exp(model->dAlpha * sqrt(model->dPWTime[new_index]));
      model->dSigmaFxDown[i] =
          vol * exp(-model->dAlpha * sqrt(model->dPWTime[new_index]));
    } else {
      model->dSigmaFxUp[i] =
          vol *
          exp(model->dAlpha * sqrt((1.0 - exp(-2.0 * model->dLambda *
                                              model->dPWTime[new_index])) /
                                   (2.0 * model->dLambda)));
      model->dSigmaFxDown[i] =
          vol *
          exp(-model->dAlpha * sqrt((1.0 - exp(-2.0 * model->dLambda *
                                               model->dPWTime[new_index])) /
                                    (2.0 * model->dLambda)));
    }
    t2 = model->dPWTime[i];

    if (fabs(model->dAlpha) < TINY) {
      /* transform correlations */
      if (CalibConst->cNumParams->iMappingMethod == 1) {
        err = FxBetaDLM_correl_mapping(
            model->dSigmaFx[i], model->dB0, model->dC0, model->dTStar,
            model->dSigmaDom[i], model->dLambdaDom, model->dSigmaFor[i],
            model->dLambdaFor, model->dCorrDomFor[i], model->dCorrDomFxInput[i],
            model->dCorrForFxInput[i], !i ? 0.0 : model->dPWTime[i - 1],
            model->dPWTime[i], varFFx, &(model->dCorrDomFx[i]),
            &(model->dCorrForFx[i]));

      } else if (CalibConst->cNumParams->iMappingMethod == 2) {
        err = FxBetaDLM_correl_mapping_forward(
            CalibConst->cNumParams->dFFxCorrelTolerance, model->dPWTime[i],
            !i ? 0.0 : model->dPWTime[i - 1], model->dB0, model->dC0, varFFx,
            model->dTStar, model->dSigmaDom[i], model->dLambdaDom,
            model->dSigmaFor[i], model->dLambdaFor, model->dSigmaFx3F[i],
            model->dSigmaFx[i], model->dCorrDomFor[i],
            model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
            &(model->dCorrDomFx[i]), &(model->dCorrForFx[i]));

      } else {
        err = "Mapping Method not available";
        goto FREE_RETURN;
      }

      if (err)
        goto FREE_RETURN;

      /* update var ffx */
      varFFx += model->dSigmaFx[i] * model->dSigmaFx[i] * (t2 - t1);
      varFFx += -2.0 * model->dCorrForFx[i] * model->dSigmaFor[i] *
                    model->dSigmaFx[i] *
                    Etha_Func(model->dLambdaFor, model->dTStar, t1, t2) +
                2.0 * model->dCorrDomFx[i] * model->dSigmaDom[i] *
                    model->dSigmaFx[i] *
                    Etha_Func(model->dLambdaDom, model->dTStar, t1, t2);
      varFFx += model->dSigmaFor[i] * model->dSigmaFor[i] *
                    Psi_Func(model->dLambdaFor, model->dLambdaFor,
                             model->dTStar, t1, t2) +
                model->dSigmaDom[i] * model->dSigmaDom[i] *
                    Psi_Func(model->dLambdaDom, model->dLambdaDom,
                             model->dTStar, t1, t2) -
                2.0 * model->dCorrDomFor[i] * model->dSigmaDom[i] *
                    model->dSigmaFor[i] *
                    Psi_Func(model->dLambdaDom, model->dLambdaFor,
                             model->dTStar, t1, t2);
    } else {
      /* transform correlations */
      if (CalibConst->cNumParams->iMappingMethod == 1) {
        if (CalibConst->cNumParams->iMidCorrelation) {
          err = FxBetaDLM_correl_mapping(
              model->dSigmaFx[i], model->dB0, model->dC0, model->dTStar,
              model->dSigmaDom[i], model->dLambdaDom, model->dSigmaFor[i],
              model->dLambdaFor, model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dPWTime[i], varFFx,
              &(model->dCorrDomFx[i]), &(model->dCorrForFx[i]));

          if (err)
            goto FREE_RETURN;

          err = FxBetaDLM_correl_mapping(
              model->dSigmaFx[i], model->dB0, model->dC0, model->dTStar,
              model->dSigmaDom[i], model->dLambdaDom, model->dSigmaFor[i],
              model->dLambdaFor, model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dPWTime[i], varFFx,
              &(model->dCorrDomFxDown[i]), &(model->dCorrForFxDown[i]));

          if (err)
            goto FREE_RETURN;
        } else {
          err = FxBetaDLM_correl_mapping(
              model->dSigmaFxUp[i], model->dB0, model->dC0, model->dTStar,
              model->dSigmaDom[i], model->dLambdaDom, model->dSigmaFor[i],
              model->dLambdaFor, model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dPWTime[i], varFFxUp,
              &(model->dCorrDomFx[i]), &(model->dCorrForFx[i]));

          if (err)
            goto FREE_RETURN;

          err = FxBetaDLM_correl_mapping(
              model->dSigmaFxDown[i], model->dB0, model->dC0, model->dTStar,
              model->dSigmaDom[i], model->dLambdaDom, model->dSigmaFor[i],
              model->dLambdaFor, model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dPWTime[i], varFFxDown,
              &(model->dCorrDomFxDown[i]), &(model->dCorrForFxDown[i]));

          if (err)
            goto FREE_RETURN;
        }
      } else if (CalibConst->cNumParams->iMappingMethod == 2) {
        if (CalibConst->cNumParams->iMidCorrelation) {
          err = FxBetaDLM_correl_mapping_forward(
              CalibConst->cNumParams->dFFxCorrelTolerance, model->dPWTime[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dB0, model->dC0, varFFx,
              model->dTStar, model->dSigmaDom[i], model->dLambdaDom,
              model->dSigmaFor[i], model->dLambdaFor, model->dSigmaFx3F[i],
              model->dSigmaFx[i], model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              &(model->dCorrDomFx[i]), &(model->dCorrForFx[i]));

          if (err)
            goto FREE_RETURN;

          err = FxBetaDLM_correl_mapping_forward(
              CalibConst->cNumParams->dFFxCorrelTolerance, model->dPWTime[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dB0, model->dC0, varFFx,
              model->dTStar, model->dSigmaDom[i], model->dLambdaDom,
              model->dSigmaFor[i], model->dLambdaFor, model->dSigmaFx3F[i],
              model->dSigmaFx[i], model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              &(model->dCorrDomFxDown[i]), &(model->dCorrForFxDown[i]));

          if (err)
            goto FREE_RETURN;
        } else {
          err = FxBetaDLM_correl_mapping_forward(
              CalibConst->cNumParams->dFFxCorrelTolerance, model->dPWTime[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dB0, model->dC0,
              varFFxUp, model->dTStar, model->dSigmaDom[i], model->dLambdaDom,
              model->dSigmaFor[i], model->dLambdaFor, model->dSigmaFx3F[i],
              model->dSigmaFxUp[i], model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              &(model->dCorrDomFx[i]), &(model->dCorrForFx[i]));

          if (err)
            goto FREE_RETURN;

          err = FxBetaDLM_correl_mapping_forward(
              CalibConst->cNumParams->dFFxCorrelTolerance, model->dPWTime[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dB0, model->dC0,
              varFFxDown, model->dTStar, model->dSigmaDom[i], model->dLambdaDom,
              model->dSigmaFor[i], model->dLambdaFor, model->dSigmaFx3F[i],
              model->dSigmaFxDown[i], model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              &(model->dCorrDomFxDown[i]), &(model->dCorrForFxDown[i]));

          if (err)
            goto FREE_RETURN;
        }
      } else {
        err = "Mapping Method not available";
        goto FREE_RETURN;
      }

      /* update var ffx up*/
      varFFx += model->dSigmaFx[i] * model->dSigmaFx[i] * (t2 - t1);
      varFFx += -2.0 * model->dCorrForFx[i] * model->dSigmaFor[i] *
                    model->dSigmaFx[i] *
                    Etha_Func(model->dLambdaFor, model->dTStar, t1, t2) +
                2.0 * model->dCorrDomFx[i] * model->dSigmaDom[i] *
                    model->dSigmaFx[i] *
                    Etha_Func(model->dLambdaDom, model->dTStar, t1, t2);
      varFFx += model->dSigmaFor[i] * model->dSigmaFor[i] *
                    Psi_Func(model->dLambdaFor, model->dLambdaFor,
                             model->dTStar, t1, t2) +
                model->dSigmaDom[i] * model->dSigmaDom[i] *
                    Psi_Func(model->dLambdaDom, model->dLambdaDom,
                             model->dTStar, t1, t2) -
                2.0 * model->dCorrDomFor[i] * model->dSigmaDom[i] *
                    model->dSigmaFor[i] *
                    Psi_Func(model->dLambdaDom, model->dLambdaFor,
                             model->dTStar, t1, t2);

      varFFxUp += model->dSigmaFxUp[i] * model->dSigmaFxUp[i] * (t2 - t1);
      varFFxUp += -2.0 * model->dCorrForFx[i] * model->dSigmaFor[i] *
                      model->dSigmaFxUp[i] *
                      Etha_Func(model->dLambdaFor, model->dTStar, t1, t2) +
                  2.0 * model->dCorrDomFx[i] * model->dSigmaDom[i] *
                      model->dSigmaFxUp[i] *
                      Etha_Func(model->dLambdaDom, model->dTStar, t1, t2);
      varFFxUp += model->dSigmaFor[i] * model->dSigmaFor[i] *
                      Psi_Func(model->dLambdaFor, model->dLambdaFor,
                               model->dTStar, t1, t2) +
                  model->dSigmaDom[i] * model->dSigmaDom[i] *
                      Psi_Func(model->dLambdaDom, model->dLambdaDom,
                               model->dTStar, t1, t2) -
                  2.0 * model->dCorrDomFor[i] * model->dSigmaDom[i] *
                      model->dSigmaFor[i] *
                      Psi_Func(model->dLambdaDom, model->dLambdaFor,
                               model->dTStar, t1, t2);

      /* update var ffx down*/
      varFFxDown += model->dSigmaFxDown[i] * model->dSigmaFxDown[i] * (t2 - t1);
      varFFxDown += -2.0 * model->dCorrForFxDown[i] * model->dSigmaFor[i] *
                        model->dSigmaFxDown[i] *
                        Etha_Func(model->dLambdaFor, model->dTStar, t1, t2) +
                    2.0 * model->dCorrDomFxDown[i] * model->dSigmaDom[i] *
                        model->dSigmaFxDown[i] *
                        Etha_Func(model->dLambdaDom, model->dTStar, t1, t2);
      varFFxDown += model->dSigmaFor[i] * model->dSigmaFor[i] *
                        Psi_Func(model->dLambdaFor, model->dLambdaFor,
                                 model->dTStar, t1, t2) +
                    model->dSigmaDom[i] * model->dSigmaDom[i] *
                        Psi_Func(model->dLambdaDom, model->dLambdaDom,
                                 model->dTStar, t1, t2) -
                    2.0 * model->dCorrDomFor[i] * model->dSigmaDom[i] *
                        model->dSigmaFor[i] *
                        Psi_Func(model->dLambdaDom, model->dLambdaFor,
                                 model->dTStar, t1, t2);
    }

    t1 = t2;
  }

FREE_RETURN:
  return err;
}

Err FxBetaDLM_SetVol(void *Model_, CALIBGEN_PARAMS CalibConsts,
                     void *GlobalConst, int vol_index, double vol) {
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);

  int i;
  int last_index, new_index;
  double dt, t1, t2, cum_var;
  Err err = NULL;

  if (vol_index == 0) {
    last_index = 0;
    t1 = 0.0;
    cum_var = 0.0;
  } else {
    last_index = CalibConst->lSigIndex[vol_index - 1];
    t1 = model->dPWTime[last_index];
    cum_var = CalibConst->cPrecalc->dVarX[vol_index - 1];
  }

  new_index = CalibConst->lSigIndex[vol_index];
  t2 = model->dPWTime[new_index];
  dt = t2 - t1;

  for (i = last_index + 1; i <= new_index; i++) {
    model->dSigmaFx[i] = vol;
    if (fabs(model->dLambda) < TINY) {
      model->dSigmaFxUp[i] =
          vol * exp(model->dAlpha * sqrt(model->dPWTime[new_index]));
      model->dSigmaFxDown[i] =
          vol * exp(-model->dAlpha * sqrt(model->dPWTime[new_index]));
    } else {
      model->dSigmaFxUp[i] =
          vol *
          exp(model->dAlpha * sqrt((1.0 - exp(-2.0 * model->dLambda *
                                              model->dPWTime[new_index])) /
                                   (2.0 * model->dLambda)));
      model->dSigmaFxDown[i] =
          vol *
          exp(-model->dAlpha * sqrt((1.0 - exp(-2.0 * model->dLambda *
                                               model->dPWTime[new_index])) /
                                    (2.0 * model->dLambda)));
    }
  }

  CalibConst->cPrecalc->dVarX[vol_index] = cum_var + vol * vol * dt;

  err = FxBetaDLM_Update_ModelPrecalc(vol_index, model, CalibConst);

  if (err)
    return err;

  CalibConst->cPrecalc->iNbDone += 1;

  return NULL;
}

Err FxBetaDLM_ExtrapolVol(void *Model_, void *GlobalConst, int last_vol_index) {
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);

  int i, var_index, last_index;
  double varFFx, varFFxUp, varFFxDown, t1, t2;
  Err err = NULL;

  last_index = CalibConst->lSigIndex[last_vol_index];

  if (last_vol_index == 0) {
    var_index = 0;
    varFFx = 0.0;
    varFFxUp = 0.0;
    varFFxDown = 0.0;
    t1 = 0.0;
  } else {
    var_index = last_vol_index - 1;
    if (fabs(model->dAlpha) < TINY)
      varFFx = CalibConst->cPrecalc->dVarFFx[var_index];
    else {
      varFFxUp = CalibConst->cPrecalc->dVarFFx[var_index];
      varFFxDown = CalibConst->cPrecalc->dVarFFx_down[var_index];
    }
    t1 = model->dPWTime[CalibConst->lSigIndex[last_vol_index - 1]];
  }

  for (i = last_index + 1; i < model->iNbPWTime; i++) {
    t2 = model->dPWTime[i];

    model->dSigmaFx[i] = model->dSigmaFx[i - 1];
    model->dSigmaFxUp[i] = model->dSigmaFxUp[i - 1];
    model->dSigmaFxDown[i] = model->dSigmaFxDown[i - 1];

    if (fabs(model->dAlpha) < TINY) {
      /* transform correlations */
      if (CalibConst->cNumParams->iMappingMethod == 1) {
        err = FxBetaDLM_correl_mapping(
            model->dSigmaFx[i], model->dB0, model->dC0, model->dTStar,
            model->dSigmaDom[i], model->dLambdaDom, model->dSigmaFor[i],
            model->dLambdaFor, model->dCorrDomFor[i], model->dCorrDomFxInput[i],
            model->dCorrForFxInput[i], !i ? 0.0 : model->dPWTime[i - 1],
            model->dPWTime[i], varFFx, &(model->dCorrDomFx[i]),
            &(model->dCorrForFx[i]));

      } else if (CalibConst->cNumParams->iMappingMethod == 2) {
        err = FxBetaDLM_correl_mapping_forward(
            CalibConst->cNumParams->dFFxCorrelTolerance, model->dPWTime[i],
            !i ? 0.0 : model->dPWTime[i - 1], model->dB0, model->dC0, varFFx,
            model->dTStar, model->dSigmaDom[i], model->dLambdaDom,
            model->dSigmaFor[i], model->dLambdaFor, model->dSigmaFx3F[i],
            model->dSigmaFx[i], model->dCorrDomFor[i],
            model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
            &(model->dCorrDomFx[i]), &(model->dCorrForFx[i]));

      } else {
        err = "Mapping Method not available";
        goto FREE_RETURN;
      }

      if (err)
        goto FREE_RETURN;

      /* update var ffx */
      varFFx += model->dSigmaFx[i] * model->dSigmaFx[i] * (t2 - t1);
      varFFx += -2.0 * model->dCorrForFx[i] * model->dSigmaFor[i] *
                    model->dSigmaFx[i] *
                    Etha_Func(model->dLambdaFor, model->dTStar, t1, t2) +
                2.0 * model->dCorrDomFx[i] * model->dSigmaDom[i] *
                    model->dSigmaFx[i] *
                    Etha_Func(model->dLambdaDom, model->dTStar, t1, t2);
      varFFx += model->dSigmaFor[i] * model->dSigmaFor[i] *
                    Psi_Func(model->dLambdaFor, model->dLambdaFor,
                             model->dTStar, t1, t2) +
                model->dSigmaDom[i] * model->dSigmaDom[i] *
                    Psi_Func(model->dLambdaDom, model->dLambdaDom,
                             model->dTStar, t1, t2) -
                2.0 * model->dCorrDomFor[i] * model->dSigmaDom[i] *
                    model->dSigmaFor[i] *
                    Psi_Func(model->dLambdaDom, model->dLambdaFor,
                             model->dTStar, t1, t2);
    } else {
      /* transform correlations */
      if (CalibConst->cNumParams->iMappingMethod == 1) {
        if (CalibConst->cNumParams->iMidCorrelation) {
          err = FxBetaDLM_correl_mapping(
              model->dSigmaFx[i], model->dB0, model->dC0, model->dTStar,
              model->dSigmaDom[i], model->dLambdaDom, model->dSigmaFor[i],
              model->dLambdaFor, model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dPWTime[i], varFFx,
              &(model->dCorrDomFx[i]), &(model->dCorrForFx[i]));

          if (err)
            goto FREE_RETURN;

          err = FxBetaDLM_correl_mapping(
              model->dSigmaFx[i], model->dB0, model->dC0, model->dTStar,
              model->dSigmaDom[i], model->dLambdaDom, model->dSigmaFor[i],
              model->dLambdaFor, model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dPWTime[i], varFFx,
              &(model->dCorrDomFxDown[i]), &(model->dCorrForFxDown[i]));

          if (err)
            goto FREE_RETURN;
        } else {
          err = FxBetaDLM_correl_mapping(
              model->dSigmaFxUp[i], model->dB0, model->dC0, model->dTStar,
              model->dSigmaDom[i], model->dLambdaDom, model->dSigmaFor[i],
              model->dLambdaFor, model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dPWTime[i], varFFxUp,
              &(model->dCorrDomFx[i]), &(model->dCorrForFx[i]));

          if (err)
            goto FREE_RETURN;

          err = FxBetaDLM_correl_mapping(
              model->dSigmaFxDown[i], model->dB0, model->dC0, model->dTStar,
              model->dSigmaDom[i], model->dLambdaDom, model->dSigmaFor[i],
              model->dLambdaFor, model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dPWTime[i], varFFxDown,
              &(model->dCorrDomFxDown[i]), &(model->dCorrForFxDown[i]));

          if (err)
            goto FREE_RETURN;
        }
      } else if (CalibConst->cNumParams->iMappingMethod == 2) {
        if (CalibConst->cNumParams->iMidCorrelation) {
          err = FxBetaDLM_correl_mapping_forward(
              CalibConst->cNumParams->dFFxCorrelTolerance, model->dPWTime[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dB0, model->dC0, varFFx,
              model->dTStar, model->dSigmaDom[i], model->dLambdaDom,
              model->dSigmaFor[i], model->dLambdaFor, model->dSigmaFx3F[i],
              model->dSigmaFx[i], model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              &(model->dCorrDomFx[i]), &(model->dCorrForFx[i]));

          if (err)
            goto FREE_RETURN;

          err = FxBetaDLM_correl_mapping_forward(
              CalibConst->cNumParams->dFFxCorrelTolerance, model->dPWTime[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dB0, model->dC0, varFFx,
              model->dTStar, model->dSigmaDom[i], model->dLambdaDom,
              model->dSigmaFor[i], model->dLambdaFor, model->dSigmaFx3F[i],
              model->dSigmaFx[i], model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              &(model->dCorrForFxDown[i]), &(model->dCorrForFxDown[i]));

          if (err)
            goto FREE_RETURN;
        } else {
          err = FxBetaDLM_correl_mapping_forward(
              CalibConst->cNumParams->dFFxCorrelTolerance, model->dPWTime[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dB0, model->dC0,
              varFFxUp, model->dTStar, model->dSigmaDom[i], model->dLambdaDom,
              model->dSigmaFor[i], model->dLambdaFor, model->dSigmaFx3F[i],
              model->dSigmaFxUp[i], model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              &(model->dCorrDomFx[i]), &(model->dCorrForFx[i]));

          if (err)
            goto FREE_RETURN;

          err = FxBetaDLM_correl_mapping_forward(
              CalibConst->cNumParams->dFFxCorrelTolerance, model->dPWTime[i],
              !i ? 0.0 : model->dPWTime[i - 1], model->dB0, model->dC0,
              varFFxDown, model->dTStar, model->dSigmaDom[i], model->dLambdaDom,
              model->dSigmaFor[i], model->dLambdaFor, model->dSigmaFx3F[i],
              model->dSigmaFxDown[i], model->dCorrDomFor[i],
              model->dCorrDomFxInput[i], model->dCorrForFxInput[i],
              &(model->dCorrForFxDown[i]), &(model->dCorrForFxDown[i]));

          if (err)
            goto FREE_RETURN;
        }
      } else {
        err = "Mapping Method not available";
        goto FREE_RETURN;
      }

      /* update var ffx up*/
      varFFx += model->dSigmaFx[i] * model->dSigmaFx[i] * (t2 - t1);
      varFFx += -2.0 * model->dCorrForFx[i] * model->dSigmaFor[i] *
                    model->dSigmaFx[i] *
                    Etha_Func(model->dLambdaFor, model->dTStar, t1, t2) +
                2.0 * model->dCorrDomFx[i] * model->dSigmaDom[i] *
                    model->dSigmaFx[i] *
                    Etha_Func(model->dLambdaDom, model->dTStar, t1, t2);
      varFFx += model->dSigmaFor[i] * model->dSigmaFor[i] *
                    Psi_Func(model->dLambdaFor, model->dLambdaFor,
                             model->dTStar, t1, t2) +
                model->dSigmaDom[i] * model->dSigmaDom[i] *
                    Psi_Func(model->dLambdaDom, model->dLambdaDom,
                             model->dTStar, t1, t2) -
                2.0 * model->dCorrDomFor[i] * model->dSigmaDom[i] *
                    model->dSigmaFor[i] *
                    Psi_Func(model->dLambdaDom, model->dLambdaFor,
                             model->dTStar, t1, t2);

      varFFxUp += model->dSigmaFxUp[i] * model->dSigmaFxUp[i] * (t2 - t1);
      varFFxUp += -2.0 * model->dCorrForFx[i] * model->dSigmaFor[i] *
                      model->dSigmaFxUp[i] *
                      Etha_Func(model->dLambdaFor, model->dTStar, t1, t2) +
                  2.0 * model->dCorrDomFx[i] * model->dSigmaDom[i] *
                      model->dSigmaFxUp[i] *
                      Etha_Func(model->dLambdaDom, model->dTStar, t1, t2);
      varFFxUp += model->dSigmaFor[i] * model->dSigmaFor[i] *
                      Psi_Func(model->dLambdaFor, model->dLambdaFor,
                               model->dTStar, t1, t2) +
                  model->dSigmaDom[i] * model->dSigmaDom[i] *
                      Psi_Func(model->dLambdaDom, model->dLambdaDom,
                               model->dTStar, t1, t2) -
                  2.0 * model->dCorrDomFor[i] * model->dSigmaDom[i] *
                      model->dSigmaFor[i] *
                      Psi_Func(model->dLambdaDom, model->dLambdaFor,
                               model->dTStar, t1, t2);

      /* update var ffx down*/
      varFFxDown += model->dSigmaFxDown[i] * model->dSigmaFxDown[i] * (t2 - t1);
      varFFxDown += -2.0 * model->dCorrForFxDown[i] * model->dSigmaFor[i] *
                        model->dSigmaFxDown[i] *
                        Etha_Func(model->dLambdaFor, model->dTStar, t1, t2) +
                    2.0 * model->dCorrDomFxDown[i] * model->dSigmaDom[i] *
                        model->dSigmaFxDown[i] *
                        Etha_Func(model->dLambdaDom, model->dTStar, t1, t2);
      varFFxDown += model->dSigmaFor[i] * model->dSigmaFor[i] *
                        Psi_Func(model->dLambdaFor, model->dLambdaFor,
                                 model->dTStar, t1, t2) +
                    model->dSigmaDom[i] * model->dSigmaDom[i] *
                        Psi_Func(model->dLambdaDom, model->dLambdaDom,
                                 model->dTStar, t1, t2) -
                    2.0 * model->dCorrDomFor[i] * model->dSigmaDom[i] *
                        model->dSigmaFor[i] *
                        Psi_Func(model->dLambdaDom, model->dLambdaFor,
                                 model->dTStar, t1, t2);
    }

    t1 = t2;
  }

FREE_RETURN:
  return err;
}

Err FxBetaDLM_ExtrapolVol_WhenBetaCalib(void *Model_, void *GlobalConst,
                                        int last_vol_index) {
  Err err = NULL;

  return err;
}

Err FxBetaDLM_UpdateConstsAfterVol(void *Inst_, void *InstConst_,
                                   void *GlobalConst, void *Model_,
                                   CALIBGEN_PARAMS CalibConsts) {
  FxBetaDLM_FxOptInst *Inst = (FxBetaDLM_FxOptInst *)(Inst_);
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);
  FxBetaDLM_InstPrecalc *InstConst = (FxBetaDLM_InstPrecalc *)(InstConst_);

  Err err = NULL;

  err = FxBetaDLM_Update_ModelPrecalc(Inst->iInstIndex, model, CalibConst);

  if (err)
    return err;

  err = FxBetaDLM_Update_InstPrecalc(Inst, model, CalibConst, InstConst);

  if (err)
    return err;

  return err;
}

Err FxBetaDLM_PriceInstVol(void *Inst_, void *InstConst_, void *GlobalConst,
                           void *Model_, double *InstPrice) {
  FxBetaDLM_FxOptInst *Inst = (FxBetaDLM_FxOptInst *)(Inst_);
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);
  FxBetaDLM_InstPrecalc *InstConst = (FxBetaDLM_InstPrecalc *)(InstConst_);

  double Result[10];
  Err err = NULL;

  err = FxBetaDLM_Price_FxOptInst(Inst, model, InstConst, CalibConst->cHermite,
                                  CalibConst->cNumParams, &(Result[0]));

  if (err)
    return err;

  *InstPrice = Result[0];

  return NULL;
}

Err FxBetaDLM_GetTargetBeta(void *Inst_, void *GlobalConst, void *Model,
                            CALIBGEN_PARAMS CalibConsts, double *target) {
  FxBetaDLM_FxOptInst *Inst = (FxBetaDLM_FxOptInst *)(Inst_);

  if (Inst->iNbStrike > 2) {
    *target = Inst->dPrice[1] - Inst->dPrice[2];
  } else if (Inst->iNbStrike == 2) {
    *target = Inst->dPrice[1];
  } else {
    *target = 0.0;
  }

  return NULL;
}

Err FxBetaDLM_GetFirstGuessBeta(void *Model_, void *GlobalConst,
                                int index_param, double target, double *beta1) {
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);

  *beta1 = model->dC0;

  return NULL;
}

Err FxBetaDLM_GetSecondGuessBeta(void *Model_, void *GlobalConst,
                                 int beta_index, double beta1, double price1,
                                 double target, double *beta2) {
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);

  if (fabs(price1 - target) < CalibConst->CalibParamsForBeta->Precision) {
    *beta2 = beta1;
  } else {
    *beta2 = beta1 + 0.02;
  }

  return NULL;
}

Err FxBetaDLM_GetLimitAndLastBeta(void *Model_, CALIBGEN_PARAMS CalibParams,
                                  void *GlobalConst,

                                  int beta_index, double *last_beta,
                                  double *limit_down, double *limit_up) {
  Err err = NULL;

  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);

  *last_beta = model->dC0;

  if (CalibConst->cPrecalc->dVarFFx[CalibConst->iSmileIndex] > 0.0) {
    *limit_down =
        -1.0 / (2.0 * CalibConst->cPrecalc->dVarFFx[CalibConst->iSmileIndex]) *
        1.001;
  } else {
    *limit_down =
        -1.0 /
        (2.0 * 0.1 * 0.1 *
         model->dPWTime[CalibConst->lSigIndex[CalibConst->iSmileIndex]]) *
        1.001;
  }

  *limit_up = 100.0;

  return err;
}

Err FxBetaDLM_BumpBeta(void *Model_, void *GlobalConst, int beta_index,
                       double beta) {
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);

  model->dC0 = beta;

  return NULL;
}

Err FxBetaDLM_SetBeta(void *Model_, CALIBGEN_PARAMS CalibConsts,
                      void *GlobalConst, int beta_index, double beta) {
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);

  // double	shift;
  Err err = NULL;

  model->dC0 = beta;
  model->dEquiBeta = beta;

  /* Update the first guess for first not yet calibrated option */
  /*
  if (CalibConst->cPrecalc->iNbDone > 0)
  {
          shift =
  model->dSigmaFx[CalibConst->lSigIndex[CalibConst->cPrecalc->iNbDone - 1]]
                  - CalibConst->first_guess_for_smile;

          if (CalibConst->iNbOpt > CalibConst->cPrecalc->iNbDone)
          {
                  model->dSigmaFx[CalibConst->lSigIndex[CalibConst->cPrecalc->iNbDone]]
  += shift;
          }
  }

  CalibConst->cPrecalc->iNbDone = 0;
  */

  CalibConst->first_guess = CalibConst->first_guess_for_smile;

  /* Reset the classical extrapol function */
  CalibConst->CalibFunctionsForVol->ExtrapolParam =
      CalibConst->ExtrapolAfterSetBeta;

  err = CalibrateParamTS(CalibConst->cPrecalc->iNbDone, CalibConst->iNbOpt - 1,
                         CalibConst->cAllInst, CalibConst->cAllInstPrecalc,
                         CalibConst, model, CalibConst->CalibParams,
                         CalibConst->CalibFunctionsForVol);

  if (err)
    return err;

  return err;
}

Err FxBetaDLM_ExtrapolBeta(void *Model_, void *GlobalConst,
                           int last_beta_index) {
  return NULL;
}

Err FxBetaDLM_UpdateConstsAfterBeta(void *Inst_, void *InstConst_,
                                    void *GlobalConst, void *Model_,
                                    CALIBGEN_PARAMS CalibConsts) {
  FxBetaDLM_FxOptInst *Inst = (FxBetaDLM_FxOptInst *)(Inst_);
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);
  FxBetaDLM_InstPrecalc *InstConst = (FxBetaDLM_InstPrecalc *)(InstConst_);

  Err err = NULL;

  CalibConst->cPrecalc->iNbDone = 0;

  err = CalibrateParamTS(0, CalibConst->iSmileIndex, CalibConst->cAllInst,
                         CalibConst->cAllInstPrecalc, CalibConst, model,
                         CalibConst->CalibParams,
                         CalibConst->CalibFunctionsForVol);

  if (err)
    return err;

  err = FxBetaDLM_Update_InstPrecalc(Inst, model, CalibConst, InstConst);

  if (err)
    return err;

  return NULL;
}

Err FxBetaDLM_PriceInstBeta(void *Inst_, void *InstConst_, void *GlobalConst,
                            void *Model_, double *InstPrice) {
  FxBetaDLM_FxOptInst *Inst = (FxBetaDLM_FxOptInst *)(Inst_);
  FxBetaDLM_model *model = (FxBetaDLM_model *)(Model_);
  FxBetaDLM_CalibConst *CalibConst = (FxBetaDLM_CalibConst *)(GlobalConst);
  FxBetaDLM_InstPrecalc *InstConst = (FxBetaDLM_InstPrecalc *)(InstConst_);

  double Result[10];
  Err err = NULL;

  err = FxBetaDLM_Price_FxOptInst(Inst, model, InstConst, CalibConst->cHermite,
                                  CalibConst->cNumParams, &(Result[0]));

  if (err)
    return err;

  if (Inst->iNbStrike > 2) {
    *InstPrice = Result[1] - Result[2];
  } else if (Inst->iNbStrike == 2) {
    *InstPrice = Result[1];
  } else {
    *InstPrice = 0.0;
  }

  return NULL;
}