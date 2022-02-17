#include "srt_h_all.h"
#include "srt_h_lgmtypes.h"
//#include "opfnctns.h"
#include "diagCalibDLM.h"
#include "srt_h_allFx3F.h"
#include "swp_h_swap_pricing.h"

#include "srtaccess.h"

#include "TargetNote.h"
#include "TargetNoteCalib.h"
#include "TargetNoteProdStruct.h"
#include "math.h"

/* --------------------------------------------------------------------------------
 */
/*
        Calibrate to caplets.  If the lambda value is not input  , calibrate to
   diagonal swaptions.  The LGM2F parameters must be input.
*/
/* --------------------------------------------------------------------------------
 */

void TARN_calib_setInst(
    TARN_Struct *tarn, char *szTenor,
    int iStrikeType, /* 0: ATM  ,  1: CASH  ,  2: SWAP  ,  3: STD  ,  4:  KO */
    DiagCalibInstStruct *ptrInst) {
  /* Variable declaration */
  int i;

  /* set known values */
  ptrInst->num_ex_dates = tarn->deal.nCouponDates;
  ptrInst->end_date = tarn->deal.lvFltrEndDates[tarn->deal.nCouponDates - 1];

  /* set all the struct pointers null */
  ptrInst->cal_date = 0;
  ptrInst->end_tenor = 0;
  ptrInst->strike = 0;

  /* Allocate memory */
  ptrInst->cal_date = (int *)calloc(tarn->deal.nCouponDates, sizeof(int));
  ptrInst->strike = (double *)calloc(tarn->deal.nCouponDates, sizeof(double));
  ptrInst->strikeP1 = (double *)calloc(tarn->deal.nCouponDates, sizeof(double));
  ptrInst->strikeM1 = (double *)calloc(tarn->deal.nCouponDates, sizeof(double));
  ptrInst->end_tenor = (char **)calloc(tarn->deal.nCouponDates, sizeof(char *));

  /* set the strikes: if it is both capped and floored we calibrate ATM */
  for (i = 0; i < tarn->deal.nCouponDates; i++) {
    if (fabs(tarn->deal.dvGearing[i]) <= 1.0e-8 ||
        (tarn->deal.ivIsFloored[i] && tarn->deal.ivIsCapped[i]))
      ptrInst->strike[i] = 0.0;
    else
      ptrInst->strike[i] =
          fabs(tarn->deal.dvCoupon[i] / tarn->deal.dvGearing[i]);

    if (fabs(tarn->deal.dvGearing[i]) <= 1.0e-8 && tarn->calibration.skip_start)
      ptrInst->cal_date[i] = 0;
    else
      ptrInst->cal_date[i] = 1;

    ptrInst->strikeP1[i] = 1.0;
    ptrInst->strikeM1[i] = -1.0;
    ptrInst->end_tenor[i] = calloc(10, sizeof(char));
    strcpy(ptrInst->end_tenor[i], szTenor);
  }
}

void TARN_calib_freeInst(DiagCalibInstStruct *ptrInst) {
  int i;

  if (ptrInst->cal_date)
    free(ptrInst->cal_date);
  if (ptrInst->strike)
    free(ptrInst->strike);
  if (ptrInst->strikeP1)
    free(ptrInst->strikeP1);
  if (ptrInst->strikeM1)
    free(ptrInst->strikeM1);
  if (ptrInst->end_tenor) {
    for (i = 0; i < ptrInst->num_ex_dates; i++)
      if (ptrInst->end_tenor[i])
        free(ptrInst->end_tenor[i]);
    free(ptrInst->end_tenor);
  }
}

char *TARN_set_calib_param(TARN_Struct *tarn, TARN_AUX *aux,
                           TARN_CALIB_AUX *calib_aux) {
  /* variable declaration */
  char *err = 0;

  /* set the LM parameters */
  diag_calib_lm_params_set_default_param(&calib_aux->LM_params);
  calib_aux->LM_params.nb_iter = tarn->calibration.nb_iter_LM;
  calib_aux->LM_params.precision = tarn->calibration.precision_LM;

  /* set the numerical values */
  lgmsv_app_set_default_params_struct(&calib_aux->NumerParams);
  calib_aux->NumerParams.iNbX = tarn->calibration.iNbX;
  calib_aux->NumerParams.iNbX = tarn->calibration.iNbX;
  calib_aux->NumerParams.dParam1 = tarn->calibration.iNbSigmaXLeft;
  calib_aux->NumerParams.dParam2 = tarn->calibration.iNbSigmaXRight;
  calib_aux->NumerParams.dIntegParam = tarn->calibration.dIntegParam;
  calib_aux->NumerParams.iIntegMethod = tarn->calibration.iIntegMethod;
  calib_aux->NumerParams.dVolLimit = tarn->calibration.dVolLimit;
  calib_aux->NumerParams.iCalibLGM = tarn->calibration.iCalibLGM;
  calib_aux->NumerParams.dMinStd = tarn->calibration.dMinStd;
  calib_aux->NumerParams.dMaxStd = tarn->calibration.dMaxStd;

  /* set up the strikes and tenors */
  TARN_calib_setInst(tarn, tarn->calibration.szPrimTenor,
                     tarn->calibration.prim_param.strike_type,
                     &calib_aux->Primary);

  TARN_calib_setInst(tarn, tarn->calibration.szSecTenor,
                     tarn->calibration.sec_param.strike_type,
                     &calib_aux->Secondary);

  // added to change the end date of the swaption calibration to be always the
  // end date of the deal
  //	calib_aux->Secondary.end_date = tarn->deal.lvFundEndDates[
  //tarn->deal.nFundDates- 1 ];

  /* if the KO calibration is requested  , generate the probs and set the
   * required instruments */
  if (tarn->calibration.sec_param.strike_type == 4) {
    if (err = TARN_calc_KO_prob_2F(tarn, aux, calib_aux))
      goto FREE_RETURN;
    tarn->calibration.prim_param.strike_type = 2;
  }

/* return */
FREE_RETURN:
  return err;
}

/* ------------------------------------------------------------------------------------------------------------------
 */
/* The TARN_lgm structure */
void set_one_factor(TARN_lgm *pCal) {
  pCal->alpha = 0.0001;
  pCal->gamma = 0.0;
  pCal->rho = 0.0;
  pCal->nFactors = 1;
}

void init_tarn_lgm(TARN_lgm *pCal, int nFactors, TARN_Struct *tarn) {
  int i;

  pCal->dmVol = 0;
  pCal->sig_time = 0;
  pCal->sig = 0;
  pCal->dmTau = 0;
  pCal->smile_datas = 0;
  pCal->sigma_dates = 0;
  pCal->smile_dates = 0;
  pCal->tau_dates = 0;
  pCal->num_tau = 0;
  pCal->tau = 0;
  pCal->tstar = tarn->model.tstar;

  pCal->nlam = tarn->model.nlam;
  pCal->lam = calloc(tarn->model.nlam, sizeof(double));
  pCal->lam_time = calloc(tarn->model.nlam, sizeof(double));
  for (i = 0; i < pCal->nlam; i++) {
    pCal->lam[i] = tarn->model.lam[i];
    pCal->lam_time[i] = tarn->model.lam_time[i];
  }
  pCal->alpha = tarn->model.alpha;
  pCal->gamma = tarn->model.gamma;
  pCal->rho = tarn->model.rho;
  pCal->nFactors = 2;

  pCal->calalpha = 0;
  pCal->callambdaeps = 0;
  pCal->calrho = 0;
  pCal->calrho2 = 0;

  /* if the model is 1F overwrite the correlation parameters */
  if (nFactors == 1)
    set_one_factor(pCal);
}

void free_tarn_lgm(TARN_lgm *pCal) {
  free_and_zero(&pCal->lam);
  free_and_zero(&pCal->lam_time);
  free_and_zero(&pCal->calalpha);
  free_and_zero(&pCal->callambdaeps);
  free_and_zero(&pCal->calrho);
  free_and_zero(&pCal->calrho2);
  free_and_zero(&pCal->sig);
  free_and_zero(&pCal->sig_time);
  free_and_zero(&pCal->sigma_dates);
  free_and_zero(&pCal->smile_dates);
  free_and_zero(&pCal->tau_dates);
  free_and_zero(&pCal->tau);
  if (pCal->dmVol)
    free_dmatrix(pCal->dmVol, 0, 1, 0, pCal->num_sig - 1);
  if (pCal->dmTau)
    free_dmatrix(pCal->dmTau, 0, 1, 0, pCal->num_tau - 1);
  if (pCal->smile_datas)
    free_dmatrix(pCal->smile_datas, 0, pCal->smile_col, 0, pCal->num_smile - 1);
}

void alloc_tarn_lgm(TARN_Struct *tarn, TARN_lgm *pCal) {
  int i;

  pCal->dmVol = dmatrix(0, 1, 0, pCal->num_sig - 1);
  for (i = 0; i < pCal->num_sig; i++) {
    pCal->dmVol[0][i] = pCal->sigma_dates
                            ? pCal->sigma_dates[i]
                            : tarn->market.lToday + 365.0 * pCal->sig_time[i];
    pCal->dmVol[1][i] = pCal->sig[i];
  }
  pCal->num_tau = pCal->num_tau == 0 ? pCal->nlam : pCal->num_tau;
  pCal->dmTau = dmatrix(0, 1, 0, pCal->num_tau - 1);
  for (i = 0; i < pCal->nlam; i++) {
    pCal->dmTau[0][i] = pCal->tau_dates
                            ? pCal->tau_dates[i]
                            : tarn->market.lToday + 365.0 * pCal->lam_time[i];
    pCal->dmTau[1][i] = pCal->tau ? pCal->tau[i] : 1.0 / pCal->lam[i];
  }

  if (pCal->calalpha) {
    if (pCal->nFactors == 1)
      pCal->smile_col = 4;
    else
      pCal->smile_col = 5;
    pCal->num_smile = pCal->num_sig;
    pCal->smile_datas = dmatrix(0, pCal->smile_col, 0, pCal->num_smile - 1);
    for (i = 0; i < pCal->num_smile; i++) {
      pCal->smile_datas[0][i] =
          pCal->smile_dates ? pCal->smile_dates[i]
                            : tarn->market.lToday + 365.0 * pCal->sig_time[i];
      pCal->smile_datas[1][i] = pCal->calalpha[i];
      pCal->smile_datas[2][i] = pCal->callambdaeps[i];
      pCal->smile_datas[3][i] = pCal->calrho[i];
      if (pCal->nFactors == 2)
        pCal->smile_datas[4][i] = pCal->calrho2[i];
    }
  }
}

/* The TARN_lgm structure */
/* ------------------------------------------------------------------------------------------------------------------
 */

char *TARN_calib_LGM2F(
    /* Market and deal info */
    char *szLGM2FUnd, char *szTARN_LGM2F_UND, TARN_Struct *tarn, int nFactors,
    TARN_CALIB_AUX *calib_aux,
    /* output */
    CPD_CALIB_INST_DATA inst_data) {
  /* Variable declarations */
  char *err = 0;
  TARN_lgm cal, *pCal = &cal;
  init_tarn_lgm(pCal, nFactors, tarn);

  /* calibration */
  if (!strcmp(szLGM2FUnd, "CAL")) {
    /* call the calibration routine */
    if (err = cpd_calib_diagonal_dlm(
            /* Market info */
            tarn->market.szYieldCurve, tarn->market.szVolCurve,
            tarn->market.getCashVol, tarn->market.szVolCurveRefRate,
            /* primary info */
            tarn->calibration.szPrimFreq, tarn->calibration.szPrimBasis,
            tarn->calibration.szPrimRefRate, tarn->deal.nCouponDates,
            tarn->deal.lvFltrFixingDates, calib_aux->Primary.cal_date,
            calib_aux->Primary.end_tenor, calib_aux->Primary.end_date,
            calib_aux->Primary.strike, &(tarn->calibration.prim_param),
            /* secondary info */
            tarn->calibration.szSecFreq, tarn->calibration.szSecBasis,
            tarn->calibration.szSecRefRate, tarn->deal.nCouponDates,
            tarn->deal.lvFltrFixingDates, calib_aux->Secondary.cal_date,
            calib_aux->Secondary.end_tenor, calib_aux->Secondary.end_date,
            calib_aux->Secondary.strike, 0, /*calib_aux->Swap.weights  ,*/
            &(tarn->calibration.sec_param),
            /* Model info */
            tarn->calibration.fix_lambda, 0, /*	Force 2F lambda */
            cal.nlam,                        /*	Force single value of lambda */
            cal.lam_time, cal.lam, NULL, 2,  /*	Force 2F */
            cal.alpha, cal.gamma, cal.rho,

            0, NULL, NULL, 0, NULL, NULL,
            /*	Output */
            &cal.num_sig, &cal.sig_time, &cal.sig,

            /*	Parameters */
            &calib_aux->LM_params,
            /*	Calibration instrument data */
            inst_data))
      goto FREE_RETURN;
  } else {
    /* Get the term structure */
    double tau;
    if (err = Get_LGM2F_TermStructure(szLGM2FUnd, &cal.sig_time, &cal.sig,
                                      &cal.num_sig, &tau, &cal.alpha,
                                      &cal.gamma, &cal.rho))
      goto FREE_RETURN;
    cal.lam[0] = 1 / tau;
    cal.nlam = 1;
    cal.lam_time[0] = 10.0;
  }

  /* Copy the data into the underlying */
  alloc_tarn_lgm(tarn, pCal);

  /* Instantiate the underlying */
  if (err = SrtInitIRUnd(szTARN_LGM2F_UND, tarn->market.szYieldCurve, "LGM2F",
                         cal.num_sig, 2, cal.dmVol, cal.nlam, 2, cal.dmTau, 0.0,
                         cal.alpha, cal.gamma, cal.rho, 0.0, 0.0, 0.0, 0, 0, 0))
    goto FREE_RETURN;

/* Error handling */
FREE_RETURN:
  free_tarn_lgm(pCal);

  /* return */
  return err;
}

/* calibrate the SV model */
char *TARN_calib_LGMSV(
    /* Market and deal info */
    char *szLGMSVUnd, char *szTARN_LGMSV_UND, TARN_Struct *tarn, int nFactors,
    TARN_CALIB_AUX *calib_aux,
    /* output */
    CPD_CALIB_INST_DATA inst_data) {
  /* variable declaration */
  char *err = 0;

  TARN_lgm cal;
  init_tarn_lgm(&cal, nFactors, tarn);

  if (!strcmp(szLGMSVUnd, "CAL")) {
    /* Calibrate */
    if (err = cpd_calib_diagonal_LGMSV_new_dlm(
            tarn->market.szYieldCurve, tarn->market.szVolCurve,
            tarn->market.getCashVol, tarn->market.szVolCurveRefRate,
            /* primary info */
            tarn->calibration.szPrimFreq, tarn->calibration.szPrimBasis,
            tarn->calibration.szPrimRefRate, tarn->deal.nCouponDates,
            tarn->deal.lvFltrFixingDates, calib_aux->Primary.cal_date,
            calib_aux->Primary.end_tenor, calib_aux->Primary.end_date,
            calib_aux->Primary.strike, calib_aux->Primary.strikeP1,
            calib_aux->Primary.strikeM1, &tarn->calibration.prim_param,
            /* secondary info */
            tarn->calibration.szSecFreq, tarn->calibration.szSecBasis,
            tarn->calibration.szSecRefRate, tarn->deal.nCouponDates,
            tarn->deal.lvFltrFixingDates, calib_aux->Secondary.cal_date,
            calib_aux->Secondary.end_tenor, calib_aux->Secondary.end_date,
            calib_aux->Secondary.strike, calib_aux->Secondary.strikeP1,
            calib_aux->Secondary.strikeM1, NULL, &tarn->calibration.sec_param,
            /* Calibration */
            &tarn->calibration.lgmsv_calib_params, cal.nFactors, cal.lam,
            cal.nFactors == 1 ? tarn->model.nsmilepar1F : tarn->model.nsmilepar,
            cal.nFactors == 1 ? tarn->model.smilepartime1F
                              : tarn->model.smilepartime,
            cal.nFactors == 1 ? tarn->model.alphaepsts1F
                              : tarn->model.alphaepsts,
            cal.nFactors == 1 ? tarn->model.ldaepsts1F : tarn->model.ldaepsts,
            cal.nFactors == 1 ? tarn->model.rhoepsts1F : tarn->model.rhoepsts,
            tarn->model.tstar, tarn->model.alpha, tarn->model.gamma,
            tarn->model.rho, tarn->model.rho2epsts, &(calib_aux->NumerParams),
            &cal.num_sig, &cal.sig_time, &cal.sig, &cal.calalpha,
            &cal.callambdaeps, &cal.calrho, &cal.calrho2,
            &(calib_aux->LM_params), inst_data))
      goto FREE_RETURN;
  } else {
    /* Get the terms structure */
    if (err = irm_sv_get_term_struct_date(
            szLGMSVUnd, &cal.nFactors, &cal.num_sig, &cal.sigma_dates, &cal.sig,
            &cal.num_tau, &cal.tau_dates, &cal.tau, &cal.alpha, &cal.gamma,
            &cal.rho, &cal.num_smile, &cal.smile_dates, &cal.calalpha,
            &cal.callambdaeps, &cal.calrho, &cal.calrho2, &cal.tstar))
      goto FREE_RETURN;
  }

  /* convert the output to matrices for creating an underlying */
  alloc_tarn_lgm(tarn, &cal);

  /* Initialize the new und */
  if (err = SrtInitIRMSVUnd(
          szTARN_LGMSV_UND, tarn->market.szYieldCurve, cal.nFactors, cal.dmVol,
          cal.num_sig, 2, cal.dmTau, cal.num_tau, 2, cal.alpha, cal.gamma,
          cal.rho, cal.smile_datas, cal.num_smile, cal.smile_col, cal.tstar))
    goto FREE_RETURN;

FREE_RETURN:
  free_tarn_lgm(&cal);

  return err;
}

char *TARN_calib(TARN_Struct *tarn, TARN_AUX *aux) {
  /* Variable declarations */
  char *err = 0;
  TARN_CALIB_AUX CALIB_AUX, *calib_aux = &CALIB_AUX;

  /* set up the calibraton parameters */
  if (err = TARN_set_calib_param(tarn, aux, calib_aux))
    goto FREE_RETURN;

  /* LGM2F calibration */
  if (tarn->pricing.iPrice2FCV) {
    if (err = TARN_calib_LGM2F(tarn->market.szLGM2FUnd,
                               tarn->output.szTARN_LGM2F_UND, tarn, 2,
                               calib_aux, &tarn->output.inst_data_lgm2F))
      goto FREE_RETURN;

    /* LGM1F calibration */
    if (err = TARN_calib_LGM2F(tarn->market.szLGM1FUnd,
                               tarn->output.szTARN_LGM1F_UND, tarn, 1,
                               calib_aux, &tarn->output.inst_data_lgm1F))
      goto FREE_RETURN;

    /* LGM1FSV calibration */
    if (err = TARN_calib_LGMSV(tarn->market.szLGM1FSVUnd,
                               tarn->output.szTARN_LGM1FSV_UND, tarn, 1,
                               calib_aux, &tarn->output.inst_data_lgm1FSV))
      goto FREE_RETURN;
  }

  /* LGMSV calibration */
  if (tarn->pricing.iPriceSV)
    if (err = TARN_calib_LGMSV(tarn->market.szLGMSVUnd,
                               tarn->output.szTARN_LGMSV_UND, tarn, 2,
                               calib_aux, &tarn->output.inst_data_lgmSV))
      goto FREE_RETURN;

/* Error handling */
FREE_RETURN:

  return err;
}

/* --------------------------------------------------------------------------------
 */
/* Calculate the knock-out probs */
char *TARN_calc_KO_prob_2F(TARN_Struct *tarn, TARN_AUX *aux,
                           TARN_CALIB_AUX *calib_aux) {
  /* Variable declarations */
  char *err = 0;
  //	TARN_CALIB_AUX CALIB_AUX_TEMP  , *calib_aux_temp = &CALIB_AUX_TEMP;
  char szTempUnd[128] = "TEMP_UND_2F";
  double dCumulProb = 0.0, *dvProbs = 0;
  int i;
  long *lvIndex = 0;

  /* set up the calibraton parameters */
  //	if ( err = TARN_set_calib_param( tarn  , aux  , calib_aux_temp )  )
  //		goto FREE_RETURN;

  /* LGM2F calibration */
  if (err = TARN_calib_LGM2F("CAL", szTempUnd, tarn, 2, calib_aux, 0))
    goto FREE_RETURN;

  /* run a Monte-Carlo to calculate the KO probs */
  if (err = TargetNoteMC_KO_2F(szTempUnd, tarn, aux, tarn->output.dmKO_prob))
    goto FREE_RETURN;

  /* sort the output probs */
  dvProbs = calloc(tarn->deal.nCouponDates, sizeof(double));
  lvIndex = calloc(tarn->deal.nCouponDates, sizeof(long));
  for (i = 0; i < tarn->deal.nCouponDates; i++) {
    //		calib_aux->Swap.weights[i] = tarn->output.dmKO_prob[i][0];
    calib_aux->Secondary.cal_date[i] = 0;
    dvProbs[i] = tarn->output.dmKO_prob[i][0];
    lvIndex[i] = i;
  }

  if (err = indexx_dl(dvProbs, lvIndex, tarn->deal.nCouponDates))
    goto FREE_RETURN;

  /* find the KOs with the highest probs and set them as calibration instruments
   */
  for (i = tarn->deal.nCouponDates - 1; i >= 0 && dCumulProb < 0.9; i--) {
    dCumulProb += dvProbs[lvIndex[i]];
    calib_aux->Secondary.cal_date[lvIndex[i]] = 1;
  }

/* Error handling */
FREE_RETURN:
  free_and_zero(&dvProbs);
  free_and_zero(&lvIndex);
  return err;
}
