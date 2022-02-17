
#include "DiagCalibDLMNew.h"
#include "DiagCalibDLM.h"
#include "DiagCalibGen.h"
#include "Fx3FUtils.h"
#include "cpdcalib.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"

void free_LGM_model(LGM_MODEL model) {
  if (model->dPWTime)
    free(model->dPWTime);
  if (model->dSigma)
    free(model->dSigma);
  if (model->dLambda)
    free(model->dLambda);
  if (model->dLambdaTime)
    free(model->dLambdaTime);
  if (model->dInitLambda)
    free(model->dInitLambda);
  if (model->dNewLambda)
    free(model->dNewLambda);
  if (model->dNewLambda2)
    free(model->dNewLambda2);
}

Err Initialise_Factors_LGM(CALIBCPNSCHEDULEDLM CalibCpnSchedule,
                           CALIBEXESCHEDULEDLM sCalibExeSchedule,
                           LGM_MODEL model, LGM_MODELFACTORS sModelFactors) {
  sModelFactors->dCpn_G1 = calloc(CalibCpnSchedule->iNCpn, sizeof(double));
  sModelFactors->dEx_G1 = calloc(sCalibExeSchedule->iNExe, sizeof(double));
  sModelFactors->dExpFact1 = calloc(sCalibExeSchedule->iNExe, sizeof(double));
  sModelFactors->dZeta1 = calloc(CalibCpnSchedule->iNCpn, sizeof(double));

  if (!sModelFactors->dCpn_G1 || !sModelFactors->dExpFact1 ||
      !sModelFactors->dEx_G1 || !sModelFactors->dZeta1) {
    return "Memory allocation faillure in Initialise_Factors_LGM";
  }

  if (model->iNbFactor >= 2) {
    sModelFactors->dCpn_G2 = calloc(CalibCpnSchedule->iNCpn, sizeof(double));
    sModelFactors->dEx_G2 = calloc(sCalibExeSchedule->iNExe, sizeof(double));
    sModelFactors->dExpFact2 = calloc(sCalibExeSchedule->iNExe, sizeof(double));
    sModelFactors->dExpFact12 =
        calloc(sCalibExeSchedule->iNExe, sizeof(double));
    sModelFactors->dZeta2 = calloc(sCalibExeSchedule->iNExe, sizeof(double));
    sModelFactors->dZeta12 = calloc(sCalibExeSchedule->iNExe, sizeof(double));

    if (!sModelFactors->dCpn_G2 || !sModelFactors->dExpFact2 ||
        !sModelFactors->dExpFact12 || !sModelFactors->dEx_G2 ||
        !sModelFactors->dZeta2 || !sModelFactors->dZeta12) {
      return "Memory allocation faillure in Initialise_Factors_LGM";
    }
  }

  return NULL;
}

void Free_Factors_LGM(LGM_MODELFACTORS sModelFactors) {
  if (sModelFactors->dCpn_G1)
    free(sModelFactors->dCpn_G1);
  if (sModelFactors->dEx_G1)
    free(sModelFactors->dEx_G1);
  if (sModelFactors->dExpFact1)
    free(sModelFactors->dExpFact1);
  if (sModelFactors->dZeta1)
    free(sModelFactors->dZeta1);
  if (sModelFactors->dCpn_G2)
    free(sModelFactors->dCpn_G2);
  if (sModelFactors->dEx_G2)
    free(sModelFactors->dEx_G2);
  if (sModelFactors->dExpFact2)
    free(sModelFactors->dExpFact2);
  if (sModelFactors->dExpFact12)
    free(sModelFactors->dExpFact12);
  if (sModelFactors->dZeta2)
    free(sModelFactors->dZeta2);
  if (sModelFactors->dZeta12)
    free(sModelFactors->dZeta12);
}

Err Initialise_AllParams_LGM(
    LGM_MODEL model,

    long *lSigLongIndex, CALIBCPNSCHEDULEDLM CalibLongCpnSchedule,
    ALLCALININSTRUMENTSDLM AllCalibLongInst, CPD_DIAG_CALIB_PARAM paraml,

    long *lSigShortIndex, CALIBCPNSCHEDULEDLM CalibShortCpnSchedule,
    ALLCALININSTRUMENTSDLM AllCalibShortInst, CPD_DIAG_CALIB_PARAM params,

    int fix_lambda,

    LGM_ALLPARAMS AllParams) {
  Err err = NULL;

  AllParams->CalibParams = calloc(1, sizeof(CALIBGEN_Params));
  AllParams->CalibParamsLambda = calloc(1, sizeof(CALIBGEN_Params));
  AllParams->sLongFactors = calloc(1, sizeof(LGM_ModelFactors));
  AllParams->sShortFactors = calloc(1, sizeof(LGM_ModelFactors));

  if (!AllParams->CalibParams || !AllParams->CalibParamsLambda ||
      !AllParams->sLongFactors || !AllParams->sShortFactors) {
    err = "Memory allocation faillure in Initialise_AllParams";
    return err;
  }

  AllParams->long_param = paraml;
  AllParams->short_param = params;

  err = Initialise_CalibParams(paraml->use_jumps, paraml->precision,
                               paraml->nb_iter_max, 0, 1, 0, 0, 0.0,
                               AllParams->CalibParams);

  if (err)
    goto FREE_RETURN;

  err = Initialise_CalibParams(params->use_jumps, params->precision,
                               params->nb_iter_max, 1, !fix_lambda, 0, 1, 0.001,
                               AllParams->CalibParamsLambda);

  if (err)
    goto FREE_RETURN;

  AllParams->lSigLongIndex = lSigLongIndex;
  AllParams->iNbLongInst = AllCalibLongInst->iNbInst;

  err = Initialise_Factors_LGM(CalibLongCpnSchedule,
                               AllCalibLongInst->sCalibExeSchedule, model,
                               AllParams->sLongFactors);

  if (err)
    goto FREE_RETURN;

  if (!fix_lambda) {
    lSigShortIndex = lSigShortIndex;
    AllParams->iNbShortInst = AllCalibShortInst->iNbInst;

    err = Initialise_Factors_LGM(CalibShortCpnSchedule,
                                 AllCalibShortInst->sCalibExeSchedule, model,
                                 AllParams->sShortFactors);

    if (err)
      goto FREE_RETURN;
  }

  AllParams->CalibCpnLongSchedule = AllCalibLongInst->sCalibCpnSchedule;
  AllParams->CalibExeLongSchedule = AllCalibLongInst->sCalibExeSchedule;
  AllParams->CalibCpnShortSchedule = AllCalibShortInst->sCalibCpnSchedule;
  AllParams->CalibExeShortSchedule = AllCalibShortInst->sCalibExeSchedule;

FREE_RETURN:

  return err;
}

void Free_AllParams_LGM(LGM_ALLPARAMS AllParams) {
  if (AllParams->CalibParams) {
    free(AllParams->CalibParams);
  }

  if (AllParams->CalibParamsLambda) {
    free(AllParams->CalibParamsLambda);
  }

  if (AllParams->sLongFactors) {
    Free_Factors_LGM(AllParams->sLongFactors);
    free(AllParams->sLongFactors);
  }

  if (AllParams->sShortFactors) {
    Free_Factors_LGM(AllParams->sShortFactors);
    free(AllParams->sShortFactors);
  }
}

void Calculate_Lambda_Factors(CALIBCPNSCHEDULEDLM sSchedule,
                              CALIBEXESCHEDULEDLM sCalibExeSchedule,
                              LGM_MODEL model, LGM_MODELFACTORS sModelFactors) {
  int i;
  double t1, t2;

  static_lgmsetupG_tauts(model->iNbLam, model->dLambdaTime, model->dNewLambda,
                         sSchedule->iNCpn, sSchedule->dCpnTime,
                         sModelFactors->dCpn_G1, sCalibExeSchedule->iNExe,
                         sCalibExeSchedule->dExeTimes, sModelFactors->dEx_G1);

  if (model->iNbFactor == 2) {
    static_lgmsetupG_tauts(
        model->iNbLam, model->dLambdaTime, model->dNewLambda2, sSchedule->iNCpn,
        sSchedule->dCpnTime, sModelFactors->dCpn_G2, sCalibExeSchedule->iNExe,
        sCalibExeSchedule->dExeTimes, sModelFactors->dEx_G2);
  }

  t1 = 0.0;

  for (i = 0; i < sCalibExeSchedule->iNExe; i++) {
    t2 = sCalibExeSchedule->dExeTimes[i];

    sModelFactors->dExpFact1[i] = export_lgmcalcexpfact_tauts(
        t1, t2, model->iNbLam, model->dLambdaTime, model->dLambda, 0.0);

    if (model->iNbFactor == 2) {
      sModelFactors->dExpFact2[i] =
          export_lgmcalcexpfact_tauts(t1, t2, model->iNbLam, model->dLambdaTime,
                                      model->dLambda, model->dGamma) /
          sModelFactors->dExpFact1[i];

      sModelFactors->dExpFact12[i] =
          export_lgmcalcexpfact_tauts(t1, t2, model->iNbLam, model->dLambdaTime,
                                      model->dLambda, 0.5 * model->dGamma) /
          sModelFactors->dExpFact1[i];
    }
  }
}

/*	Calibrate lgm: main function */
/*	New version: calibrates not necessarily to diagonal
                with lambda calibration */
Err cpd_calib_diagonal_dlm_new(
    /*	Market */
    char *yc_name,        /*	Name of the yield curve */
    char *vol_curve_name, /*	Name of the market vol curve */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),

    /* Get Cash Vol ref */
    char *vol_ref_rate_name,

    /* Long Instruments */
    char *instr_long_freq, /*	Frequency and basis of instruments */
    char *instr_long_basis,
    char *long_ref_rate_name, /*	Name of the reference rate */

    int num_ex_datesl, /*	Long Exercise dates */
    long *ex_datel_,   /*	Supposed to be sorted */
    int *cal_datel,    /*	1: use ex_date as calibration date  , 0: don't */
    char **end_tenorl, /*	Tenors of the underlying instruments */
    long end_datel,    /*	End date for diagonal */
    double *strikel_,  /*	Strikes */

    CPD_DIAG_CALIB_PARAM paraml,

    /* Short Instruments */
    char *instr_short_freq, /*	Frequency and basis of instruments */
    char *instr_short_basis,
    char *short_ref_rate_name, /*	Name of the reference rate */

    int num_ex_datess, /*	Short Exercise dates */
    long *ex_dates_,   /*	Supposed to be sorted */
    int *cal_dates,    /*	1: use ex_date as calibration date  , 0: don't */
    char **end_tenors, /*	Tenors of the underlying instruments */
    long end_dates,    /*	End date for diagonal */
    double *strikes_,  /*	Strikes */
    double *weights_,  /*	Weights on secondary instruments */

    CPD_DIAG_CALIB_PARAM params,

    /*	Model */
    int fix_lambda,
    int one_factor_equi, /*	Calibrate 2 Factor to 1 Factor price */
    int nlam,            /*	Lambda TS: may NOT be changed in the process */
    double lam_time[], double lam[], double lam_shift[],
    int nfactor,  /*	Number of factors */
    double alpha, /*	Alpha  , Gamma  , Rho (2F only) */
    double gamma, double rho,

    /*	Output */
    int *num_sig, /*	Answer */
    double **sig_time, double **sig,

    /*	Parameters */
    DIAG_CALIB_LM_PARAMS lm_params,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA
        inst_data) /*	NULL = don't save calibration instrument data */
{
  Err err = NULL;

  long today;
  SrtCurvePtr yc_ptr;
  int i, j;

  long ex_datel[MAX_INST], start_datel[MAX_INST], theo_end_datel[MAX_INST],
      act_end_datel[MAX_INST], ex_dates[MAX_INST], start_dates[MAX_INST],
      theo_end_dates[MAX_INST], act_end_dates[MAX_INST];

  double strikel[MAX_INST], ex_timel[MAX_INST], strikes[MAX_INST],
      weightl[MAX_INST], weights[MAX_INST], ex_times[MAX_INST];

  CalibCpnScheduleDLM *CalibCpnLongSchedule = NULL,
                      *CalibCpnShortSchedule = NULL;

  CalibExeScheduleDLM *CalibExeLongSchedule = NULL,
                      *CalibExeShortSchedule = NULL;

  AllCalibInstrumentsDLM *AllCalibLongInst = NULL, *AllCalibShortInst = NULL,
                         *AllCalibShortInstOne = NULL;

  void **AllLongInstruments = NULL, **AllShortInstruments = NULL;

  LGM_AllParams *AllParams = NULL;

  LGM_model *model = NULL;

  long *lSigLongIndex = NULL, *lSigShortIndex = NULL;
  double time0 = 100;

  long total_und_long, total_und_short;
  int shift_lam;
  clock_t time1, time2;

  time1 = clock();

  /*	Copy data so as not to change the original */
  *sig_time = NULL;
  *sig = NULL;

  /* Initialisations */
  CalibCpnLongSchedule = calloc(1, sizeof(CalibCpnScheduleDLM));
  CalibCpnShortSchedule = calloc(1, sizeof(CalibCpnScheduleDLM));
  CalibExeLongSchedule = calloc(1, sizeof(CalibExeScheduleDLM));
  CalibExeShortSchedule = calloc(1, sizeof(CalibExeScheduleDLM));
  AllCalibLongInst = calloc(1, sizeof(AllCalibInstrumentsDLM));
  AllCalibShortInst = calloc(1, sizeof(AllCalibInstrumentsDLM));
  AllCalibShortInstOne = calloc(1, sizeof(AllCalibInstrumentsDLM));
  AllParams = calloc(1, sizeof(LGM_AllParams));
  model = calloc(1, sizeof(LGM_model));

  if (!CalibCpnLongSchedule || !CalibCpnShortSchedule ||
      !CalibExeLongSchedule || !CalibExeShortSchedule || !AllCalibLongInst ||
      !AllCalibShortInst || !AllParams || !AllCalibShortInstOne || !model) {
    err = "Memory allocation faillure in cpd_calib_diagonal_LGM_dlm";
    goto FREE_RETURN;
  }

  /* Get today */
  yc_ptr = lookup_curve(yc_name);

  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }

  today = get_today_from_curve(yc_ptr);

  /* Initialisation */
  memset(start_datel, 0, MAX_INST * sizeof(long));
  memset(start_dates, 0, MAX_INST * sizeof(long));

  /* Long instruments */
  AllocateCalibExeSchedule(&(ex_datel[0]), &(ex_timel[0]), &(start_datel[0]),
                           &(theo_end_datel[0]), &(act_end_datel[0]),
                           &(strikel[0]), NULL, NULL, &(weightl[0]),
                           CalibExeLongSchedule);

  err = Construct_CalibSchedule(
      yc_name, today, instr_long_freq, instr_long_basis, num_ex_datesl,
      ex_datel_, cal_datel, end_tenorl, end_datel, strikel_, NULL, NULL, NULL,
      CalibCpnLongSchedule, CalibExeLongSchedule);

  if (err)
    goto FREE_RETURN;

  err = Reduct_ExeSchedule(CalibCpnLongSchedule, CalibExeLongSchedule,
                           cal_datel, paraml);

  if (err)
    goto FREE_RETURN;

  err = AllocateAllCalibInst(CalibExeLongSchedule->iNExe, 1, AllCalibLongInst);

  if (err)
    goto FREE_RETURN;

  AllLongInstruments =
      (void **)calloc(CalibExeLongSchedule->iNExe, sizeof(void *));

  if (!AllLongInstruments) {
    err = "Memory allocation faillure in cpd_calib_diagonal_LGM_new_dlm";
    goto FREE_RETURN;
  }

  for (i = 0; i < CalibExeLongSchedule->iNExe; i++) {
    AllLongInstruments[i] = (void *)(&AllCalibLongInst->sCalibInst[i]);
  }

  err = Calculate_CalibInst(today, yc_name, vol_curve_name, get_cash_vol,
                            vol_ref_rate_name, instr_long_freq,
                            instr_long_basis, long_ref_rate_name, 1, paraml,
                            CalibCpnLongSchedule, CalibExeLongSchedule,
                            AllCalibLongInst, AllCalibLongInst);

  if (err)
    goto FREE_RETURN;

  /* Short instruments */
  if (!fix_lambda) {
    /* allocation */
    AllocateCalibExeSchedule(&(ex_dates[0]), &(ex_times[0]), &(start_dates[0]),
                             &(theo_end_dates[0]), &(act_end_dates[0]),
                             &(strikes[0]), NULL, NULL, &(weights[0]),
                             CalibExeShortSchedule);

    err = Construct_CalibSchedule(
        yc_name, today, instr_short_freq, instr_short_basis, num_ex_datess,
        ex_dates_, cal_dates, end_tenors, end_dates, strikes_, NULL, NULL,
        weights_, CalibCpnShortSchedule, CalibExeShortSchedule);

    if (err)
      goto FREE_RETURN;

    err = Reduct_ExeSchedule(CalibCpnShortSchedule, CalibExeShortSchedule,
                             cal_dates, params);

    if (err)
      goto FREE_RETURN;

    err = AllocateAllCalibInst(CalibExeShortSchedule->iNExe, 1,
                               AllCalibShortInst);

    if (err)
      goto FREE_RETURN;

    AllShortInstruments =
        (void **)calloc(CalibExeShortSchedule->iNExe, sizeof(void *));

    if (!AllShortInstruments) {
      err = "Memory allocation faillure in cpd_calib_diagonal_LGM_new_dlm";
      goto FREE_RETURN;
    }

    for (i = 0; i < CalibExeShortSchedule->iNExe; i++) {
      AllShortInstruments[i] = (void *)(&AllCalibShortInst->sCalibInst[i]);
    }

    err = Calculate_CalibInst(today, yc_name, vol_curve_name, get_cash_vol,
                              vol_ref_rate_name, instr_short_freq,
                              instr_short_basis, short_ref_rate_name, 1, params,
                              CalibCpnShortSchedule, CalibExeShortSchedule,
                              AllCalibShortInst, AllCalibLongInst);

    if (err)
      goto FREE_RETURN;
  } else {
    AllShortInstruments = (void **)calloc(1, sizeof(void *));

    if (!AllShortInstruments) {
      err = "Memory allocation faillure in cpd_calib_diagonal_LGM_new_dlm";
      goto FREE_RETURN;
    }
  }

  /*	2.)	Fill the model */
  model->lToday = today;
  model->iNbFactor = nfactor;
  model->dAlpha = alpha;
  model->dGamma = gamma;
  model->dRho = rho;

  model->iNbPWTime = CalibExeLongSchedule->iNExe;
  model->dPWTime = calloc(model->iNbPWTime, sizeof(double));
  model->iNbLam = nlam;
  model->dLambdaTime = calloc(model->iNbLam, sizeof(double));
  model->dInitLambda = calloc(model->iNbLam, sizeof(double));
  model->dNewLambda = calloc(model->iNbLam, sizeof(double));
  model->dNewLambda2 = calloc(model->iNbLam, sizeof(double));

  if (!model->dPWTime || !model->dLambdaTime || !model->dInitLambda ||
      !model->dNewLambda || !model->dNewLambda2) {
    err = "Memory allocation faillure in cpd_calib_diagonal_dlm_new";
    goto FREE_RETURN;
  }

  model->iNbLam = nlam;
  memcpy(model->dLambdaTime, lam_time, nlam * sizeof(double));
  memcpy(model->dInitLambda, lam, nlam * sizeof(double));
  memcpy(model->dNewLambda, lam, nlam * sizeof(double));

  if (model->iNbFactor == 2) {
    for (i = 0; i < nlam; i++) {
      model->dNewLambda2[i] = model->dNewLambda[i] + model->dGamma;
    }
  }

  /* Init with Primary Exe Times */
  memcpy(model->dPWTime, CalibExeLongSchedule->dExeTimes,
         CalibExeLongSchedule->iNExe * sizeof(double));

  if (model->iNbLam > 1) {
    /* Add the TS of Lambda */
    num_f_concat_vector(&model->iNbPWTime, &model->dPWTime, nlam, lam_time);
  }

  if (!fix_lambda) {
    /* Add the Secondary Exe Times */
    num_f_concat_vector(&model->iNbPWTime, &model->dPWTime,
                        CalibExeShortSchedule->iNExe,
                        CalibExeShortSchedule->dExeTimes);
  }

  num_f_sort_vector(model->iNbPWTime, model->dPWTime);
  num_f_unique_vector(&model->iNbPWTime, model->dPWTime);

  model->dSigma = calloc(model->iNbPWTime, sizeof(double));
  model->dLambda = calloc(model->iNbPWTime, sizeof(double));

  lSigLongIndex = lvector(0, CalibExeLongSchedule->iNExe - 1);

  if (!fix_lambda) {
    lSigShortIndex = lvector(0, CalibExeShortSchedule->iNExe - 1);
  }

  if (!model->dSigma || !model->dLambda || !lSigLongIndex ||
      (!lSigShortIndex && !fix_lambda)) {
    err = "Memory allocation faillure (2) in cpd_calibSV_approx";
    goto FREE_RETURN;
  }

  for (i = 0; i < model->iNbPWTime; i++) {
    j = Get_Index(model->dPWTime[i], lam_time, nlam);
    model->dLambda[i] = lam[j];
  }

  for (i = 0; i < CalibExeLongSchedule->iNExe; i++) {
    lSigLongIndex[i] = Get_Index(CalibExeLongSchedule->dExeTimes[i],
                                 model->dPWTime, model->iNbPWTime);
  }

  if (!fix_lambda) {
    for (i = 0; i < CalibExeShortSchedule->iNExe; i++) {
      lSigShortIndex[i] = Get_Index(CalibExeShortSchedule->dExeTimes[i],
                                    model->dPWTime, model->iNbPWTime);
    }
  }

  /* Check Lambda Shift */
  if (params && fabs(params->lambda_shift) > 1.0E-08) {
    shift_lam = 1;

    if (fix_lambda) {
      /* We shift before calibration */
      for (i = 0; i < model->iNbPWTime; i++) {
        model->dLambda[i] += params->lambda_shift;
      }

      for (i = 0; i < model->iNbLam; i++) {
        model->dInitLambda[i] += params->lambda_shift;
        model->dNewLambda[i] += params->lambda_shift;

        if (model->iNbFactor == 2)
          model->dNewLambda2[i] += params->lambda_shift;
      }

      shift_lam = 0;
    }
  } else {
    shift_lam = 0;
  }

  err = Initialise_AllParams_LGM(model, lSigLongIndex, CalibCpnLongSchedule,
                                 AllCalibLongInst, paraml, lSigShortIndex,
                                 CalibCpnShortSchedule, AllCalibShortInst,
                                 params, fix_lambda, AllParams);

  if (err)
    goto FREE_RETURN;

  if (!fix_lambda) {
    /* Set the Lambda Sensitivity */
    total_und_long = 0;

    for (i = 0; i < AllCalibLongInst->iNbInst; i++) {
      total_und_long +=
          CalibCpnLongSchedule
              ->lCpnDate[AllCalibLongInst->sCalibInst[i].iEndCpn] -
          CalibCpnLongSchedule
              ->lCpnDate[AllCalibLongInst->sCalibInst[i].iStartCpn];
    }

    total_und_long /= AllCalibLongInst->iNbInst;

    total_und_short = 0;

    for (i = 0; i < AllCalibShortInst->iNbInst; i++) {
      total_und_short +=
          CalibCpnShortSchedule
              ->lCpnDate[AllCalibShortInst->sCalibInst[i].iEndCpn] -
          CalibCpnShortSchedule
              ->lCpnDate[AllCalibShortInst->sCalibInst[i].iStartCpn];
    }

    total_und_short /= AllCalibShortInst->iNbInst;

    if (total_und_long > total_und_short) {
      AllParams->sens_lambda = 1.0;
    } else {
      AllParams->sens_lambda = -1.0;
    }
  }

  /* Set the Calibration Vol Function */
  AllParams->CalibFunctionsForVol.GetTarget = GetTargetVol_LGM;
  AllParams->CalibFunctionsForVol.BumpParam = BumpVol_LGM;
  AllParams->CalibFunctionsForVol.ExtrapolParam = ExtrapolVol_LGM;
  AllParams->CalibFunctionsForVol.GetFirstGuess = GetFirstGuessVol_LGM;
  AllParams->CalibFunctionsForVol.GetLimitAndLastParam = GetLimitAndLastVol_LGM;
  AllParams->CalibFunctionsForVol.GetSecondGuess = GetSecondGuessVol_LGM;
  AllParams->CalibFunctionsForVol.PriceInst = PriceInstVol_LGM;
  AllParams->CalibFunctionsForVol.SetParam = SetVol_LGM;
  AllParams->CalibFunctionsForVol.UpdateConstsAfterParam =
      UpdateParamsAfterVol_LGM;

  AllParams->CalibFunctionsForLambda.GetTarget = GetTargetLambda_LGM;
  AllParams->CalibFunctionsForLambda.BumpParam = BumpLambda_LGM;
  AllParams->CalibFunctionsForLambda.ExtrapolParam = ExtrapolLambda_LGM;
  AllParams->CalibFunctionsForLambda.GetFirstGuess = GetFirstGuessLambda_LGM;
  AllParams->CalibFunctionsForLambda.GetLimitAndLastParam =
      GetLimitAndLastLambda_LGM;
  AllParams->CalibFunctionsForLambda.GetSecondGuess = GetSecondGuessLambda_LGM;
  AllParams->CalibFunctionsForLambda.PriceInst = PriceInstLambda_LGM;
  AllParams->CalibFunctionsForLambda.SetParam = SetLambda_LGM;
  AllParams->CalibFunctionsForLambda.UpdateConstsAfterParam =
      UpdateParamsAfterLambda_LGM;

  err = CalibrateParamTS(0, 0, AllShortInstruments, AllShortInstruments,
                         AllParams, model, AllParams->CalibParamsLambda,
                         &(AllParams->CalibFunctionsForLambda));

  if (err)
    goto FREE_RETURN;

  time2 = clock();
  smessage("Calibration time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /*	5.)	Create the result */
  *num_sig = model->iNbPWTime;
  *sig_time = calloc(*num_sig, sizeof(double));
  *sig = calloc(*num_sig, sizeof(double));

  if (!(*sig_time) || !(*sig)) {
    err = "Memory allocation faillure in LGMSVCalibApprox";
    goto FREE_RETURN;
  }

  for (i = 0; i < *num_sig; i++) {
    (*sig_time)[i] = model->dPWTime[i];
    (*sig)[i] = model->dSigma[i];
  }

  /*	6.)	Save instrument data if required */
  if (inst_data) {
    num_ex_datesl = CalibExeLongSchedule->iNExe;

    inst_data->num_inst = num_ex_datesl;

    inst_data->exer_dates_long = (long *)calloc(num_ex_datesl, sizeof(long));
    inst_data->start_dates = (long *)calloc(num_ex_datesl, sizeof(long));
    inst_data->end_dates = (long *)calloc(num_ex_datesl, sizeof(long));
    inst_data->long_strikes = (double *)calloc(num_ex_datesl, sizeof(double));
    inst_data->market_prices_long =
        (double *)calloc(num_ex_datesl, sizeof(double));

    if (!fix_lambda) {
      num_ex_datess = CalibExeShortSchedule->iNExe;

      inst_data->num_insts = num_ex_datess;
      inst_data->exer_dates_short = (long *)calloc(num_ex_datess, sizeof(long));
      inst_data->start_datess = (long *)calloc(num_ex_datess, sizeof(long));
      inst_data->end_datess = (long *)calloc(num_ex_datess, sizeof(long));
      inst_data->short_strikes =
          (double *)calloc(num_ex_datess, sizeof(double));
      inst_data->short_weights =
          (double *)calloc(num_ex_datess, sizeof(double));
      inst_data->market_prices_short =
          (double *)calloc(num_ex_datess, sizeof(double));

      if (!inst_data->exer_dates_long || !inst_data->start_dates ||
          !inst_data->end_dates || !inst_data->long_strikes ||
          !inst_data->market_prices_long || !inst_data->exer_dates_short ||
          !inst_data->start_datess || !inst_data->end_datess ||
          !inst_data->short_weights || !inst_data->short_strikes ||
          !inst_data->market_prices_short) {
        err = "Allocation error (4) in cpd_calib_diagonal";
        goto FREE_RETURN;
      }
    } else {
      inst_data->num_insts = 0;
      if (!inst_data->exer_dates_long || !inst_data->start_dates ||
          !inst_data->end_dates || !inst_data->long_strikes ||
          !inst_data->market_prices_long) {
        err = "Allocation error (4) in cpd_calib_diagonal";
        goto FREE_RETURN;
      }
    }

    for (i = 0; i < num_ex_datesl; i++) {
      inst_data->exer_dates_long[i] =
          (long)(CalibExeLongSchedule->dExeTimes[i] * DAYS_IN_YEAR + 0.5) +
          today;
      inst_data->start_dates[i] =
          CalibCpnLongSchedule->lCpnDate[CalibExeLongSchedule->iStartCpn[i]];
      inst_data->end_dates[i] =
          CalibCpnLongSchedule->lCpnDate[CalibExeLongSchedule->iEndCpn[i]];
      inst_data->long_strikes[i] = AllCalibLongInst->sCalibInst[i].dStrike[0];
      inst_data->market_prices_long[i] = CalibExeLongSchedule->dPrices[i];
    }

    if (!fix_lambda) {
      for (i = 0; i < num_ex_datess; i++) {
        inst_data->exer_dates_short[i] =
            (long)(CalibExeShortSchedule->dExeTimes[i] * DAYS_IN_YEAR + 0.5) +
            today;
        inst_data->start_datess[i] =
            CalibCpnShortSchedule
                ->lCpnDate[CalibExeShortSchedule->iStartCpn[i]];
        inst_data->end_datess[i] =
            CalibCpnShortSchedule->lCpnDate[CalibExeShortSchedule->iEndCpn[i]];
        inst_data->short_strikes[i] =
            AllCalibShortInst->sCalibInst[i].dStrike[0];
        inst_data->short_weights[i] = 1.0;
        inst_data->market_prices_short[i] = CalibExeShortSchedule->dPrices[i];
      }
    }
  }

FREE_RETURN:

  if (err) {
    if (*sig_time)
      free(*sig_time);
    *sig_time = NULL;

    if (*sig)
      free(*sig);
    *sig = NULL;

    if (inst_data) {
      cpd_free_calib_inst_data(inst_data);
    }
  }

  if (lSigLongIndex)
    free_lvector(lSigLongIndex, 0, CalibExeLongSchedule->iNExe - 1);
  if (lSigShortIndex)
    free_lvector(lSigShortIndex, 0, CalibExeShortSchedule->iNExe - 1);

  if (CalibCpnLongSchedule) {
    free(CalibCpnLongSchedule);
  }

  if (CalibCpnShortSchedule) {
    free(CalibCpnShortSchedule);
  }

  if (CalibExeLongSchedule) {
    free(CalibExeLongSchedule);
  }

  if (CalibExeShortSchedule) {
    free(CalibExeShortSchedule);
  }

  if (AllCalibLongInst) {
    FreeAllCalibInst(AllCalibLongInst);
    free(AllCalibLongInst);
  }

  if (AllLongInstruments) {
    free(AllLongInstruments);
  }

  if (AllCalibShortInst) {
    FreeAllCalibInst(AllCalibShortInst);
    free(AllCalibShortInst);
  }

  if (AllShortInstruments) {
    free(AllShortInstruments);
  }

  if (AllCalibShortInstOne) {
    free(AllCalibShortInstOne);
  }

  if (AllParams) {
    Free_AllParams_LGM(AllParams);
    free(AllParams);
  }

  if (model) {
    free_LGM_model(model);
    free(model);
  }

  return err;
}

/*	Calibrate lgm: main function */
/*	New version: calibrates not necessarily to diagonal
                with lambda calibration */

Err GetTargetVol_LGM(void *Inst_, void *GlobalParam, void *Model,
                     CALIBGEN_PARAMS CalibConsts, double *target) {
  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);
  CALIBINSTRUMENTDLM Inst = (CALIBINSTRUMENTDLM)Inst_;

  *target = Inst->dPrice[0];

  if (AllParams->long_param->vega_prec) {
    *target /= Inst->dVega[0];
  }

  return NULL;
}

Err GetFirstGuessVol_LGM(void *Model, void *GlobalParam, int vol_index,
                         double target, double *vol1) {
  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);

  *vol1 = pow(model->dSigma[AllParams->lSigLongIndex[vol_index]], 2);

  return NULL;
}

Err GetSecondGuessVol_LGM(void *Model, void *GlobalParam, int vol_index,
                          double var1, double price1, double target,
                          double *var2) {
  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);

  double t1, t2, dt;
  double cum_var_temp, cum_var_lgm;

  if (vol_index == 0) {
    t1 = 0.0;
    cum_var_temp = 0.0;
    cum_var_lgm = 0.0;
  } else {
    t1 = model->dPWTime[AllParams->lSigLongIndex[vol_index - 1]];
    cum_var_lgm = AllParams->cum_var_lgm;
  }

  t2 = model->dPWTime[AllParams->lSigLongIndex[vol_index]];

  dt = t2 - t1;

  cum_var_temp = cum_var_lgm + var1 * dt;

  *var2 = (target * target / price1 / price1 * cum_var_temp - cum_var_lgm) / dt;

  if (*var2 < 0.0) {
    *var2 = 1.0E-08;
  }

  return NULL;
}

Err BumpVol_LGM(void *Model, void *GlobalParam, int vol_index, double var) {
  Err err = NULL;

  int i;
  int last_index, new_index;
  double vol;

  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);

  if (vol_index == 0) {
    last_index = 0;
  } else {
    last_index = AllParams->lSigLongIndex[vol_index - 1] + 1;
  }

  new_index = AllParams->lSigLongIndex[vol_index];

  vol = sqrt(var);

  for (i = last_index; i <= new_index; i++) {
    model->dSigma[i] = vol;
  }

  return err;
}

Err SetVol_LGM(void *Model, CALIBGEN_PARAMS CalibConsts, void *GlobalParam,
               int vol_index, double var) {
  Err err = NULL;

  int i;
  int last_index, new_index;
  double t1, t2, dt;
  double vol, cum_var_lgm;

  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);

  if (vol_index == 0) {
    last_index = -1;
    t1 = 0.0;
    cum_var_lgm = 0.0;
  } else {
    last_index = AllParams->lSigLongIndex[vol_index - 1];
    t1 = model->dPWTime[last_index];
    cum_var_lgm = AllParams->cum_var_lgm;
  }

  new_index = AllParams->lSigLongIndex[vol_index];
  t2 = model->dPWTime[new_index];

  dt = t2 - t1;

  vol = sqrt(var);

  for (i = last_index + 1; i <= new_index; i++) {
    model->dSigma[i] = vol;
  }

  AllParams->cum_var_lgm = cum_var_lgm + var * dt;

  AllParams->index_vol++;

  return err;
}

Err GetLimitAndLastVol_LGM(void *Model, CALIBGEN_PARAMS CalibParams,
                           void *GlobalParam,

                           int vol_index, double *last_vol, double *limit_down,
                           double *limit_up) {
  Err err = NULL;

  double t1, t2, dt;
  int last_index, new_index;
  double sig_limit_down_lgm, sig_limit_up_lgm;

  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);

  if (vol_index == 0) {
    *last_vol = 0.00001;
    *limit_down = 0.000000001;
    *limit_up = 10000.0;
    return err;
  }

  last_index = AllParams->lSigLongIndex[vol_index - 1];
  new_index = AllParams->lSigLongIndex[vol_index];

  t1 = model->dPWTime[last_index];
  t2 = model->dPWTime[new_index];
  dt = t2 - t1;

  sig_limit_down_lgm =
      sqrt(AllParams->cum_var_lgm * AllParams->long_param->min_fact / t1);
  sig_limit_up_lgm =
      sqrt(AllParams->cum_var_lgm / AllParams->long_param->max_fact / t1);

  *last_vol = pow(model->dSigma[last_index], 2);
  *limit_down = pow(sig_limit_down_lgm, 2);
  *limit_up = pow(sig_limit_up_lgm, 2);

  return err;
}

Err ExtrapolVol_LGM(void *Model, void *GlobalParam, int last_vol_index) {
  Err err = NULL;
  int last_index;
  int i;

  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);

  last_index = AllParams->lSigLongIndex[last_vol_index];

  for (i = last_index + 1; i < model->iNbPWTime; i++) {
    model->dSigma[i] = model->dSigma[last_index];
  }

  return err;
}

Err UpdateParamsAfterVol_LGM(void *Inst_, void *InstParam, void *GlobalParam,
                             void *Model, CALIBGEN_PARAMS CalibConsts) {
  Err err = NULL;

  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);
  double last_zeta1, last_zeta2, last_zeta12, var;

  /* Update the Zetas */
  if (AllParams->index_vol == 0) {
    last_zeta1 = 0.0;

    if (model->iNbFactor == 2) {
      last_zeta2 = 0.0;
      last_zeta12 = 0.0;
    }
  } else {
    last_zeta1 = AllParams->sLongFactors->dZeta1[AllParams->index_vol - 1];

    if (model->iNbFactor == 2) {
      last_zeta2 = AllParams->sLongFactors->dZeta2[AllParams->index_vol - 1];
      last_zeta12 = AllParams->sLongFactors->dZeta12[AllParams->index_vol - 1];
    }
  }

  var = model->dSigma[AllParams->lSigLongIndex[AllParams->index_vol]];
  var *= var;

  AllParams->sLongFactors->dZeta1[AllParams->index_vol] =
      last_zeta1 +
      var * AllParams->sLongFactors->dExpFact1[AllParams->index_vol];

  if (model->iNbFactor == 2) {
    AllParams->sLongFactors->dZeta2[AllParams->index_vol] =
        last_zeta2 +
        var * AllParams->sLongFactors->dExpFact2[AllParams->index_vol];
    AllParams->sLongFactors->dZeta12[AllParams->index_vol] =
        last_zeta12 +
        var * AllParams->sLongFactors->dExpFact12[AllParams->index_vol];
  }

  return err;
}

Err PriceInstVol_LGM(void *Inst_, void *InstParam, void *GlobalParam,
                     void *Model, double *InstPrice) {
  Err err = NULL;

  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);
  CALIBINSTRUMENTDLM Inst = (CALIBINSTRUMENTDLM)Inst_;

  if (model->iNbFactor == 1) {
    *InstPrice =
        lgmopval1F(Inst->iNbCoupon, Inst->dCpnValues,
                   &AllParams->sLongFactors->dCpn_G1[Inst->iStartCpn],
                   AllParams->sLongFactors->dZeta1[AllParams->index_vol],
                   AllParams->sLongFactors->dEx_G1[AllParams->index_vol]) *
        Inst->dWeight;
  } else {
    *InstPrice =
        lgmopval2F(Inst->iNbCoupon, Inst->dCpnValues,
                   &AllParams->sLongFactors->dCpn_G1[Inst->iStartCpn],
                   &AllParams->sLongFactors->dCpn_G2[Inst->iStartCpn],
                   AllParams->sLongFactors->dZeta1[AllParams->index_vol],
                   AllParams->sLongFactors->dZeta2[AllParams->index_vol],
                   AllParams->sLongFactors->dZeta12[AllParams->index_vol],
                   AllParams->sLongFactors->dEx_G1[AllParams->index_vol],
                   AllParams->sLongFactors->dEx_G2[AllParams->index_vol]) *
        Inst->dWeight;
  }

  if (AllParams->long_param->vega_prec) {
    *InstPrice /= Inst->dVega[0];
  }

  return err;
}

Err GetTargetLambda_LGM(void *Inst_, void *GlobalConst, void *Model,

                        CALIBGEN_PARAMS CalibConsts, double *target) {
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalConst);
  CalibInstrumentDLM *NextInst;
  CALIBINSTRUMENTDLM Inst = (CALIBINSTRUMENTDLM)Inst_;
  double total;

  NextInst = Inst;
  total = 0.0;

  while (NextInst != NULL) {
    total += NextInst->dPrice[0] * NextInst->dWeight;
    NextInst = NextInst->NextInst;
  }

  *target = total;

  return NULL;
}

Err GetFirstGuessLambda_LGM(void *Model, void *GlobalConst, int index_lambda,
                            double target, double *lambda1) {
  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalConst);

  *lambda1 = 0.0;

  return NULL;
}

Err GetSecondGuessLambda_LGM(void *Model, void *GlobalParam, int index_lambda,
                             double lambda1, double price1, double target,
                             double *lambda2) {
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);

  if (target > price1) {
    *lambda2 = 0.01;
  } else {
    *lambda2 = -0.01;
  }

  return NULL;
}

Err GetLimitAndLastLambda_LGM(void *Model, CALIBGEN_PARAMS CalibConsts,
                              void *GlobalParam, int index_lambda,
                              double *last_lambda, double *limit_down,
                              double *limit_up) {
  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);
  int i;

  *last_lambda = 0.0;

  *limit_down = 999;
  *limit_up = -999;

  for (i = 0; i < model->iNbLam; i++) {
    *limit_down = DMIN(
        AllParams->short_param->lambda_min - model->dNewLambda[i], *limit_down);
    *limit_up = DMAX(AllParams->short_param->lambda_max - model->dNewLambda[i],
                     *limit_up);
  }

  return NULL;
}

Err BumpLambda_LGM(void *Model, void *GlobalParam, int index_lambda,
                   double lambda) {
  Err err = NULL;

  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);
  int i;

  for (i = 0; i < model->iNbLam; i++) {
    model->dNewLambda[i] = model->dInitLambda[i] + lambda;

    if (model->iNbFactor == 2) {
      model->dNewLambda2[i] = model->dNewLambda[i] + model->dGamma;
    }
  }

  AllParams->index_vol = 0;

  return err;
}

Err SetLambda_LGM(void *Model, CALIBGEN_PARAMS CalibConsts, void *GlobalParam,
                  int index_lambda, double lambda) {
  Err err = NULL;

  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);
  int i;

  if (CalibConsts->do_calib) {
    for (i = 0; i < model->iNbLam; i++) {
      model->dNewLambda[i] = model->dInitLambda[i] + lambda;

      if (model->iNbFactor == 2) {
        model->dNewLambda2[i] = model->dNewLambda[i] + model->dGamma;
      }
    }
  }

  AllParams->index_vol = 0;

  return err;
}

Err PriceInstLambda_LGM(void *Inst_, void *InstParam, void *GlobalParam,
                        void *Model, double *InstPrice) {
  Err err = NULL;

  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);
  CALIBINSTRUMENTDLM Inst = (CALIBINSTRUMENTDLM)Inst_;

  CalibInstrumentDLM *NextInst;
  double total;
  int i;

  total = 0.0;
  NextInst = Inst;
  i = 0;

  while (NextInst != NULL) {
    if (model->iNbFactor == 1) {
      total += lgmopval1F(Inst->iNbCoupon, Inst->dCpnValues,
                          &AllParams->sShortFactors->dCpn_G1[Inst->iStartCpn],
                          AllParams->sShortFactors->dZeta1[i],
                          AllParams->sShortFactors->dEx_G1[i]) *
               Inst->dWeight;
    } else {
      total += lgmopval2F(Inst->iNbCoupon, Inst->dCpnValues,
                          &AllParams->sShortFactors->dCpn_G1[Inst->iStartCpn],
                          &AllParams->sShortFactors->dCpn_G2[Inst->iStartCpn],
                          AllParams->sShortFactors->dZeta1[i],
                          AllParams->sShortFactors->dZeta2[i],
                          AllParams->sShortFactors->dZeta12[i],
                          AllParams->sShortFactors->dEx_G1[i],
                          AllParams->sShortFactors->dEx_G2[i]) *
               Inst->dWeight;
    }

    NextInst = NextInst->NextInst;
    i++;
  }

  *InstPrice = total;

  return err;
}

Err ExtrapolLambda_LGM(void *Model, void *GlobalParam, int last_smile_index) {
  return NULL;
}

Err UpdateParamsAfterLambda_LGM(void *Inst_, void *InstParam, void *GlobalParam,
                                void *Model, CALIBGEN_PARAMS CalibConsts) {
  Err err = NULL;
  int i, j, last_index, nb_inst;
  double time0 = 100;
  double ex_zeta[MAX_CPN];
  double temp, fact;

  LGM_MODEL model = (LGM_MODEL)(Model);
  LGM_ALLPARAMS AllParams = (LGM_ALLPARAMS)(GlobalParam);
  CalibInstrumentDLM *NextInst;
  CALIBINSTRUMENTDLM Inst = (CALIBINSTRUMENTDLM)Inst_;

  /* Setup the Initial G for Long Instruments */
  Calculate_Lambda_Factors(AllParams->CalibCpnLongSchedule,
                           AllParams->CalibExeLongSchedule, model,
                           AllParams->sLongFactors);

  if (CalibConsts->do_calib) {
    /* Setup the Initial G for Short Instruments */
    Calculate_Lambda_Factors(AllParams->CalibCpnShortSchedule,
                             AllParams->CalibExeShortSchedule, model,
                             AllParams->sShortFactors);
  }

  /* Recalibrate the Volatility */
  err = CalibrateParamTS(
      0, nb_inst - 1, AllParams->AllLongInst, AllParams->AllLongInst, AllParams,
      model, AllParams->CalibParams, &(AllParams->CalibFunctionsForVol));

  if (CalibConsts->do_calib) {
    /* Update Zeta for Short */
    static_interpolate_zeta(AllParams->CalibExeLongSchedule->iNExe,
                            AllParams->CalibExeLongSchedule->dExeTimes,
                            AllParams->sLongFactors->dZeta1, model->iNbLam,
                            model->dLambdaTime, model->dLambda,
                            AllParams->CalibExeShortSchedule->iNExe,
                            AllParams->CalibExeShortSchedule->dExeTimes,
                            AllParams->sShortFactors->dZeta1);

    if (model->iNbFactor == 2) {
      static_lgmcalczeta2zeta12_tauts(
          AllParams->CalibExeShortSchedule->iNExe,
          AllParams->CalibExeShortSchedule->dExeTimes,
          AllParams->sShortFactors->dZeta1, model->iNbLam, model->dLambdaTime,
          model->dLambda, model->dAlpha, model->dGamma, model->dRho,
          AllParams->sShortFactors->dZeta2, AllParams->sShortFactors->dZeta12);
    }
  }

FREE_RETURN:

  return err;
}