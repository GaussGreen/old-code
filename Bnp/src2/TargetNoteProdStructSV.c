#include "LGMSVUtil.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"

#include "TargetNoteProdStructSV.h"

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Memory management for the SV parameters */

char *init_TARN_MC_AUX_SV(TARN_MC_AUX_SV **mc_aux_sv) {
  char *err;

  /* allocate the memory */
  if (!(*mc_aux_sv = (TARN_MC_AUX_SV *)malloc(sizeof(TARN_MC_AUX_SV))))
    return "Memory allocation error in TARN_MC_AUX_SV";

  /* LGM SV parameters */
  (*mc_aux_sv)->sigma = 0;
  (*mc_aux_sv)->alpha = 0;
  (*mc_aux_sv)->rho = 0;
  (*mc_aux_sv)->lameps = 0;
  (*mc_aux_sv)->lvleps = 0;
  /*  parameters */
  (*mc_aux_sv)->dLGMalpha = 0;
  (*mc_aux_sv)->dLGMgamma = 0;
  (*mc_aux_sv)->dLGMrho = 0;
  (*mc_aux_sv)->dRho2 = 0;
  /*  derived quantities */
  (*mc_aux_sv)->logdff_star = 0;
  (*mc_aux_sv)->gam1_star = 0;
  (*mc_aux_sv)->gam2_star = 0;
  (*mc_aux_sv)->gam1_2_star = 0;
  (*mc_aux_sv)->gam2_2_star = 0;
  (*mc_aux_sv)->gam12_star = 0;

  /*	Initialisation of the LGMSV model */
  init_NULL_LGMSV_model(&(*mc_aux_sv)->model);
  if (err = Fill_lgmSV_defaultParam(&(*mc_aux_sv)->Params))
    goto FREE_RETURN;

FREE_RETURN:
  return err;
}

void free_TARN_MC_AUX_SV(TARN_MC_AUX_SV *mc_aux_sv) {
  /* LGM SV parameters */
  free_and_zero(&mc_aux_sv->sigma);
  free_and_zero(&mc_aux_sv->alpha);
  free_and_zero(&mc_aux_sv->rho);
  free_and_zero(&mc_aux_sv->lameps);
  free_and_zero(&mc_aux_sv->lvleps);
  /*  parameters */
  free_and_zero(&mc_aux_sv->dLGMalpha);
  free_and_zero(&mc_aux_sv->dLGMgamma);
  free_and_zero(&mc_aux_sv->dLGMrho);
  free_and_zero(&mc_aux_sv->dRho2);
  /* derived quantities */
  free_and_zero(&mc_aux_sv->logdff_star);
  free_and_zero(&mc_aux_sv->gam1_star);
  free_and_zero(&mc_aux_sv->gam2_star);
  free_and_zero(&mc_aux_sv->gam1_2_star);
  free_and_zero(&mc_aux_sv->gam2_2_star);
  free_and_zero(&mc_aux_sv->gam12_star);
  /* model */
  free_LGMSV_model(&mc_aux_sv->model);
}

char *alloc_TARN_MC_AUX_SV(TARN_Struct *tarn, TARN_AUX *aux,
                           TARN_MC_AUX *mc_aux) {
  /* variable declaration */
  char *err;
  int iIsTwoF;
  TARN_MC_AUX_SV *mc_aux_sv;

  /* allocate the memory */
  if (err = init_TARN_MC_AUX_SV(&mc_aux_sv))
    goto FREE_RETURN;
  mc_aux->ptrModel = (TARN_MC_AUX_SV *)mc_aux_sv;

  /* Get the LGMSV model */
  if (err = Get_LGMSV_model(mc_aux->und_name, &mc_aux_sv->model))
    goto FREE_RETURN;

  /* Find the model type */
  iIsTwoF = (mc_aux_sv->model.iOne2F == 2);
  mc_aux->tnModel = iIsTwoF ? TN_2F_SV : TN_1F_SV;

  /*	Next  , set the event dates at each of the floater fixing dates and
   * allocate the time steps */
  if (err = set_TARN_TimeSteps_SV(tarn, aux, mc_aux))
    goto FREE_RETURN;

  /* Allocate the memory for the event variables */
  err = "Memory allocation error";
  mc_aux_sv->logdff_star = (double *)calloc(mc_aux->nEvents, sizeof(double));
  mc_aux_sv->gam1_star = (double *)calloc(mc_aux->nEvents, sizeof(double));
  mc_aux_sv->gam1_2_star = (double *)calloc(mc_aux->nEvents, sizeof(double));
  if (!mc_aux_sv->logdff_star || !mc_aux_sv->gam1_star ||
      !mc_aux_sv->gam1_2_star)
    goto FREE_RETURN;

  if (iIsTwoF) {
    mc_aux_sv->gam2_star = (double *)calloc(mc_aux->nEvents, sizeof(double));
    mc_aux_sv->gam2_2_star = (double *)calloc(mc_aux->nEvents, sizeof(double));
    mc_aux_sv->gam12_star = (double *)calloc(mc_aux->nEvents, sizeof(double));
    if (!mc_aux_sv->gam2_star || !mc_aux_sv->gam2_2_star ||
        !mc_aux_sv->gam12_star)
      goto FREE_RETURN;
  }

  /*	Allocate memory for the time step variables variables */
  mc_aux->is_event = (int *)calloc(mc_aux->nStps, sizeof(int));
  mc_aux->date = (double *)calloc(mc_aux->nStps, sizeof(double));
  mc_aux_sv->sigma = (double *)calloc(mc_aux->nStps, sizeof(double));
  mc_aux_sv->alpha = (double *)calloc(mc_aux->nStps, sizeof(double));
  mc_aux_sv->rho = (double *)calloc(mc_aux->nStps, sizeof(double));
  mc_aux_sv->lameps = (double *)calloc(mc_aux->nStps, sizeof(double));
  mc_aux_sv->lvleps = (double *)calloc(mc_aux->nStps, sizeof(double));
  if (!mc_aux->is_event || !mc_aux->date || !mc_aux_sv->sigma ||
      !mc_aux_sv->alpha || !mc_aux_sv->rho || !mc_aux_sv->lameps ||
      !mc_aux_sv->lvleps)
    goto FREE_RETURN;

  if (iIsTwoF) {
    mc_aux_sv->dLGMalpha = (double *)calloc(mc_aux->nStps, sizeof(double));
    mc_aux_sv->dLGMgamma = (double *)calloc(mc_aux->nStps, sizeof(double));
    mc_aux_sv->dLGMrho = (double *)calloc(mc_aux->nStps, sizeof(double));
    mc_aux_sv->dRho2 = (double *)calloc(mc_aux->nStps, sizeof(double));
    if (!mc_aux_sv->dLGMalpha || !mc_aux_sv->dLGMgamma || !mc_aux_sv->dLGMrho ||
        !mc_aux_sv->dRho2)
      goto FREE_RETURN;
  }
  err = 0;

FREE_RETURN:
  return err;
}

/* Memory management for the SV parameters */
/* -----------------------------------------------------------------------------------------------------------------
 */
/* Get the time steps for the event dates */
/* DO NOT include the zero time step as an event */

char *set_TARN_TimeSteps_SV(TARN_Struct *tarn, TARN_AUX *aux,
                            TARN_MC_AUX *mc_aux) {
  /* variable declaration */
  int i;
  char *err = 0;
  TARN_MC_AUX_SV *mc_aux_sv = (TARN_MC_AUX_SV *)mc_aux->ptrModel;

  /* calculate the number of event dates */
  mc_aux->nEvents = tarn->deal.nCouponDates - aux->i1stCpn;

  /*	Allocate memory for event dates */
  mc_aux->time = (double *)calloc(mc_aux->nEvents, sizeof(double));
  mc_aux->evt_tms = (double *)calloc(mc_aux->nEvents, sizeof(double));
  mc_aux->evt_dts = (double *)calloc(mc_aux->nEvents, sizeof(double));
  if (!mc_aux->time || !mc_aux->evt_tms || !mc_aux->evt_dts) {
    err = "Memory allocation error";
    goto FREE_RETURN;
  }

  /* Loop over the coupon fixing dates setting the dates and times */
  for (i = 0; i < mc_aux->nEvents; i++) {
    mc_aux->time[i] =
        (tarn->deal.lvFltrFixingDates[i + aux->i1stCpn] - tarn->market.lToday) /
        365.0;
    mc_aux->evt_tms[i] = mc_aux->time[i];
    mc_aux->evt_dts[i] = tarn->deal.lvFltrFixingDates[i + aux->i1stCpn];
  }

  /* add in the sigma dates */
  i = 0;
  mc_aux->nStps = mc_aux->nEvents;
  while (i < mc_aux_sv->model.iNbPWTime &&
         mc_aux_sv->model.dPWTime[i] < mc_aux->time[mc_aux->nEvents - 1]) {
    num_f_add_number(&(mc_aux->nStps), &(mc_aux->time),
                     mc_aux_sv->model.dPWTime[i]);
    i++;
  }

  /*	Make sure that the time vector has at least the required number of
   * points */
  if (err = fill_time_vector(&mc_aux->time, &(mc_aux->nStps), 0, NULL, 0, NULL,
                             tarn->pricing.nStepT))
    goto FREE_RETURN;

FREE_RETURN:
  return err;
}

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Set reconstruction parameters */
void setLGM1FSV_star(double dTimeEvent, LGMSV_model *model,
                     double *out_gam1_star, double *out_gam1_2_star) {
  *out_gam1_star = (1.0 - exp(model->dLambdaX * (model->dTStar - dTimeEvent))) /
                   model->dLambdaX;
  *out_gam1_2_star = 0.5 * (*out_gam1_star) * (*out_gam1_star);
}

void setLGM1FSV_event(double dTimeEvent, double dTimeDF, LGMSV_model *model,
                      double gam1_star, double gam1_2_star, double *out_gam1,
                      double *out_gam1_2) {
  double temp1;
  temp1 =
      (1.0 - exp(model->dLambdaX * (model->dTStar - dTimeEvent - dTimeDF))) /
      model->dLambdaX;
  *out_gam1 = temp1 - gam1_star;
  *out_gam1_2 = 0.5 * temp1 * temp1 - gam1_2_star;
}

void setLGM2FSV_star(double dTimeEvent, LGMSV_model *model,
                     double *out_gam1_star, double *out_gam1_2_star,
                     double *out_gam2_star, double *out_gam2_2_star,
                     double *out_gam12_star) {
  *out_gam1_star = (1.0 - exp(model->dLambdaX * (model->dTStar - dTimeEvent))) /
                   model->dLambdaX;
  *out_gam2_star =
      (1.0 - exp(model->dLambdaX2 * (model->dTStar - dTimeEvent))) /
      model->dLambdaX2;
  *out_gam1_2_star = 0.5 * (*out_gam1_star) * (*out_gam1_star);
  *out_gam2_2_star = 0.5 * (*out_gam2_star) * (*out_gam2_star);
  *out_gam12_star = (*out_gam1_star) * (*out_gam2_star);
}

void setLGM2FSV_event(double dTimeEvent, double dTimeDF, LGMSV_model *model,
                      double gam1_star, double gam1_2_star, double gam2_star,
                      double gam2_2_star, double gam12_star, double *out_gam1,
                      double *out_gam1_2, double *out_gam2, double *out_gam2_2,
                      double *out_gam12) {
  double temp1, temp2;

  temp1 =
      (1.0 - exp(model->dLambdaX * (model->dTStar - dTimeEvent - dTimeDF))) /
      model->dLambdaX;
  temp2 =
      (1.0 - exp(model->dLambdaX2 * (model->dTStar - dTimeEvent - dTimeDF))) /
      model->dLambdaX2;
  *out_gam1 = temp1 - gam1_star;
  *out_gam2 = temp2 - gam2_star;
  *out_gam1_2 = 0.5 * temp1 * temp1 - gam1_2_star;
  *out_gam2_2 = 0.5 * temp2 * temp2 - gam2_2_star;
  *out_gam12 = temp1 * temp2 - gam12_star;
}

setLGMSV_event(int iEvt, /* MC event */
               int iFund, TARN_MC_AUX *mc_aux, TARN_EVENT *event) {
  /* Variable declaration */
  TARN_MC_AUX_SV *mc_aux_sv = (TARN_MC_AUX_SV *)mc_aux->ptrModel;
  TARN_EVENT_SV *event_sv = (TARN_EVENT_SV *)event->ptrModel;

  /* set the 1F values */
  if (mc_aux->tnModel == TN_1F_SV) {
    setLGM1FSV_event(
        mc_aux->evt_tms[iEvt], event->fund_start_tms[iFund], &mc_aux_sv->model,
        mc_aux_sv->gam1_star[iEvt], mc_aux_sv->gam1_2_star[iEvt],
        &event_sv->fund_start_gam1[iFund], &event_sv->fund_start_gam1_2[iFund]);

    setLGM1FSV_event(
        mc_aux->evt_tms[iEvt], event->fund_end_tms[iFund], &mc_aux_sv->model,
        mc_aux_sv->gam1_star[iEvt], mc_aux_sv->gam1_2_star[iEvt],
        &event_sv->fund_end_gam1[iFund], &event_sv->fund_end_gam1_2[iFund]);
  }
  /* set the 2F values if necessary */
  else {
    setLGM2FSV_event(
        mc_aux->evt_tms[iEvt], event->fund_start_tms[iFund], &mc_aux_sv->model,
        mc_aux_sv->gam1_star[iEvt], mc_aux_sv->gam1_2_star[iEvt],
        mc_aux_sv->gam2_star[iEvt], mc_aux_sv->gam2_2_star[iEvt],
        mc_aux_sv->gam12_star[iEvt], &event_sv->fund_start_gam1[iFund],
        &event_sv->fund_start_gam1_2[iFund], &event_sv->fund_start_gam2[iFund],
        &event_sv->fund_start_gam2_2[iFund],
        &event_sv->fund_start_gam12[iFund]);

    setLGM2FSV_event(
        mc_aux->evt_tms[iEvt], event->fund_end_tms[iFund], &mc_aux_sv->model,
        mc_aux_sv->gam1_star[iEvt], mc_aux_sv->gam1_2_star[iEvt],
        mc_aux_sv->gam2_star[iEvt], mc_aux_sv->gam2_2_star[iEvt],
        mc_aux_sv->gam12_star[iEvt], &event_sv->fund_end_gam1[iFund],
        &event_sv->fund_end_gam1_2[iFund], &event_sv->fund_end_gam2[iFund],
        &event_sv->fund_end_gam2_2[iFund], &event_sv->fund_end_gam12[iFund]);
  }
}

/* Set reconstruction parameters */
/* -----------------------------------------------------------------------------------------------------------------
 */

char *check_ptr(int n, void *v, ...) {
  va_list ap;
  int i;
  char *err = "Memory allocation error";

  va_start(ap, v);
  for (i = 0; i < 0; i++) {
    if (!va_arg(ap, void *))
      goto FREE_RETURN;
  }
  err = 0;

FREE_RETURN:
  va_end(ap);
  return err;
}

char *
fillTargetNoteFundingSV(int iEvtMC, /* MC event date---not deal event date */
                        TARN_Struct *tarn, TARN_AUX *aux, TARN_MC_AUX *mc_aux,
                        TARN_EVENT *event) {
  /* Variable declaration */
  int iRelFund, iFund, iCpn, iEvtDeal;
  long lEvent;
  char *err = 0;
  long *lvFundEndIndex = aux->ivFundEndIndex,
       *lvFundStartIndex = aux->ivFundStartIndex;
  TARN_MC_AUX_SV *mc_aux_sv = 0;
  TARN_EVENT_SV *event_sv = (TARN_EVENT_SV *)event->ptrModel;

  /* add the skipped coupons to get the event index */
  iEvtDeal = iEvtMC + aux->i1stCpn;

  /* Variable initializations */
  lEvent = tarn->deal.lvFltrFixingDates[iEvtDeal];
  iCpn = iEvtDeal + 1;

  /* Set up the target note */
  if (iCpn >= tarn->deal.nCouponDates)
    event->num_fund_df = 0;
  else
    event->num_fund_df = lvFundEndIndex[iCpn] - lvFundStartIndex[iCpn] + 1;
  if (event->num_fund_df > 0) {
    event->fund_start_tms = dvector(
        0, event->num_fund_df -
               1); // xStr.evt[i].evt->dft[fund_idx]; // discount factor times
    event->fund_start_dts = dvector(
        0, event->num_fund_df -
               1); // xStr.evt[i].evt->dfd[fund_idx]; // discount factor dates
    event->fund_start_logdff = dvector(0, event->num_fund_df - 1);
    event->fund_start = dvector(0, event->num_fund_df - 1);
    event->fund_end_tms = dvector(
        0, event->num_fund_df -
               1); // xStr.evt[i].evt->dft[fund_idx]; // discount factor times
    event->fund_end_dts = dvector(
        0, event->num_fund_df -
               1); // xStr.evt[i].evt->dfd[fund_idx]; // discount factor dates
    event->fund_end_logdff = dvector(0, event->num_fund_df - 1);
    event->fund_end = dvector(0, event->num_fund_df - 1);
    event->fund_cvgMargin = dvector(0, event->num_fund_df - 1);
    if (err = check_ptr(9, event->fund_start_tms, event->fund_start_dts,
                        event->fund_start, event->fund_end_tms,
                        event->fund_end_dts, event->fund_end_logdff,
                        event->fund_end, event->fund_cvgMargin))
      goto FREE_RETURN;

    event_sv->fund_start_gam1 = dvector(0, event->num_fund_df - 1);
    event_sv->fund_start_gam1_2 = dvector(0, event->num_fund_df - 1);
    event_sv->fund_end_gam1 = dvector(0, event->num_fund_df - 1);
    event_sv->fund_end_gam1_2 = dvector(0, event->num_fund_df - 1);

    if (mc_aux->tnModel == TN_2F_SV) {
      event_sv->fund_start_gam2 = dvector(0, event->num_fund_df - 1);
      event_sv->fund_start_gam2_2 = dvector(0, event->num_fund_df - 1);
      event_sv->fund_start_gam12 = dvector(0, event->num_fund_df - 1);
      event_sv->fund_end_gam1 = dvector(0, event->num_fund_df - 1);
      event_sv->fund_end_gam2 = dvector(0, event->num_fund_df - 1);
      event_sv->fund_end_gam1_2 = dvector(0, event->num_fund_df - 1);
      event_sv->fund_end_gam2_2 = dvector(0, event->num_fund_df - 1);
      event_sv->fund_end_gam12 = dvector(0, event->num_fund_df - 1);
    }

    for (iRelFund = 0; iRelFund < event->num_fund_df; iRelFund++) {
      iFund = iRelFund + lvFundStartIndex[iCpn];
      /* start of the float period quantities */
      event->fund_start_dts[iRelFund] = tarn->deal.lvFundStartDates[iFund];
      event->fund_start_tms[iRelFund] =
          (event->fund_start_dts[iRelFund] - lEvent) / 365.0;
      event->fund_start_logdff[iRelFund] = log(swp_f_df(
          lEvent, event->fund_start_dts[iRelFund], tarn->market.szYieldCurve));

      event->fund_end_dts[iRelFund] = tarn->deal.lvFundEndDates[iFund];
      event->fund_end_tms[iRelFund] =
          (event->fund_end_dts[iRelFund] - lEvent) / 365.0;
      event->fund_end_logdff[iRelFund] = log(swp_f_df(
          lEvent, event->fund_end_dts[iRelFund], tarn->market.szYieldCurve));

      event->fund_cvgMargin[iRelFund] =
          1.0 - aux->dvFundCvg[iFund] * (tarn->deal.dvFundMargin[iFund] +
                                         tarn->deal.dvFundSpread[iFund]);

      setLGMSV_event(iEvtMC, iRelFund, mc_aux, event);
    }
  }

FREE_RETURN:
  return 0;
}

char *fillTargetNoteCouponSV(int iEvt, /* MC event date */
                             TARN_AUX *aux, TARN_MC_AUX *mc_aux,
                             TARN_EVENT *event) {
  /* Variable decalaration */
  TARN_MC_AUX_SV *mc_aux_sv = (TARN_MC_AUX_SV *)mc_aux->ptrModel;
  TARN_EVENT_SV *event_sv = (TARN_EVENT_SV *)event->ptrModel;

  /* Find out what the model is and set the model structure pointer */
  if (mc_aux->tnModel == TN_2F_SV) {

    setLGM2FSV_event(mc_aux->evt_tms[iEvt], event->fixed_pay_tms,
                     &mc_aux_sv->model, mc_aux_sv->gam1_star[iEvt],
                     mc_aux_sv->gam1_2_star[iEvt], mc_aux_sv->gam2_star[iEvt],
                     mc_aux_sv->gam2_2_star[iEvt], mc_aux_sv->gam12_star[iEvt],
                     &event_sv->fixed_pay_gam1, &event_sv->fixed_pay_gam1_2,
                     &event_sv->fixed_pay_gam2, &event_sv->fixed_pay_gam2_2,
                     &event_sv->fixed_pay_gam12);

    setLGM2FSV_event(mc_aux->evt_tms[iEvt], event->fltr_start_tms,
                     &mc_aux_sv->model, mc_aux_sv->gam1_star[iEvt],
                     mc_aux_sv->gam1_2_star[iEvt], mc_aux_sv->gam2_star[iEvt],
                     mc_aux_sv->gam2_2_star[iEvt], mc_aux_sv->gam12_star[iEvt],
                     &event_sv->fltr_start_gam1, &event_sv->fltr_start_gam1_2,
                     &event_sv->fltr_start_gam2, &event_sv->fltr_start_gam2_2,
                     &event_sv->fltr_start_gam12);

    setLGM2FSV_event(mc_aux->evt_tms[iEvt], event->fltr_end_tms,
                     &mc_aux_sv->model, mc_aux_sv->gam1_star[iEvt],
                     mc_aux_sv->gam1_2_star[iEvt], mc_aux_sv->gam2_star[iEvt],
                     mc_aux_sv->gam2_2_star[iEvt], mc_aux_sv->gam12_star[iEvt],
                     &event_sv->fltr_end_gam1, &event_sv->fltr_end_gam1_2,
                     &event_sv->fltr_end_gam2, &event_sv->fltr_end_gam2_2,
                     &event_sv->fltr_end_gam12);
  } else {

    setLGM1FSV_event(mc_aux->evt_tms[iEvt], event->fixed_pay_tms,
                     &mc_aux_sv->model, mc_aux_sv->gam1_star[iEvt],
                     mc_aux_sv->gam1_2_star[iEvt], &event_sv->fixed_pay_gam1,
                     &event_sv->fixed_pay_gam1_2);

    setLGM1FSV_event(mc_aux->evt_tms[iEvt], event->fltr_start_tms,
                     &mc_aux_sv->model, mc_aux_sv->gam1_star[iEvt],
                     mc_aux_sv->gam1_2_star[iEvt], &event_sv->fltr_start_gam1,
                     &event_sv->fltr_start_gam1_2);

    setLGM1FSV_event(mc_aux->evt_tms[iEvt], event->fltr_end_tms,
                     &mc_aux_sv->model, mc_aux_sv->gam1_star[iEvt],
                     mc_aux_sv->gam1_2_star[iEvt], &event_sv->fltr_end_gam1,
                     &event_sv->fltr_end_gam1_2);
  }

  return 0;
}
