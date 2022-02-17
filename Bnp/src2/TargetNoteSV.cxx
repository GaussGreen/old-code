/*

TargetNoteSV.cxx

Monte-Carlo pricer for the TargetNote in the LGMSV framework

*/

#include "BGMEval.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "SrtAccess.h"
#include "math.h"
#include "srt_h_all.h"

// SV files
#include "LGMSVGrfn.h"
#include "LGMSVMC.h"
#include "LGMSVUtil.h"
#include "LGMSVpde.h"
#include "LgmSVClosedForm.h"

// TargetNote files
#include "TargetNoteProdStruct.h"
#include "TargetNoteProdStructSV.h"
#include "TargetNoteSV.h"

// Extern global variable
extern TargetNote TARGET_NOTE;

/*
        nEvents is the number of event dates
        nStepT is the requested number of time steps
        out_nStp is the number of time steps
*/

/*

Payoff function for SV

*/

char *TargetNotePayoffSV1(long path_index, double evt_date, double evt_time,
                          void *func_parm,

                          /* Model data	*/
                          double ft, double psi, double v,

                          /* Vector of results to be updated */
                          int num_col,
                          double *res, // res[0] = fund        , res[1] = coupon
                          int *stop_path) {
  /* Variable declaration */
  char *err;
  int l;
  TARN_EVENT *event = (TARN_EVENT *)func_parm;
  TARN_EVENT_SV *event_sv = (TARN_EVENT_SV *)event->ptrModel;

  /* Initialize the return values */
  res[0] = 0.0;
  res[1] = 0.0;
  *stop_path = 0;

  /* see if we need to initialize anything
          if ( event->bIsFirst )
                  resetTargetNote();
  */
  /* Check that there is anything to do
          if ( ! event->do_dom )
                  return 0;
  */
  /* Reconstruct the discount factors */
  /* coupon pay */
  event->fixed_pay =
      exp(event->fixed_pay_logdff - event_sv->fixed_pay_gam1 * ft -
          event_sv->fixed_pay_gam1_2 * psi);
  /* floater */
  event->fltr_start =
      exp(event->fltr_start_logdff - event_sv->fltr_start_gam1 * ft -
          event_sv->fltr_start_gam1_2 * psi);
  event->fltr_end = exp(event->fltr_end_logdff - event_sv->fltr_end_gam1 * ft -
                        event_sv->fltr_end_gam1_2 * psi);

  /* Call the coupon calculator */
  if (err = TargetNoteCoupon(event, res, stop_path))
    return err;

  /* Check to see if it has knocked out */
  if (*stop_path)
    return 0;

  /* Reconstruct the start float discount factors */
  for (l = 0; l < event->num_fund_df; l++) {
    event->fund_start[l] =
        exp(event->fund_start_logdff[l] - event_sv->fund_start_gam1[l] * ft -
            event_sv->fund_start_gam1_2[l] * psi);
    event->fund_end[l] =
        exp(event->fund_end_logdff[l] - event_sv->fund_end_gam1[l] * ft -
            event_sv->fund_end_gam1_2[l] * psi);
  }

  /* Call the funding calculator */
  return TargetNoteFunding(event, res);
}

Err TargetNotePayoffSV2(long path_index, double evt_date, double evt_time,
                        void *func_parm,

                        double ft1, double ft2, double psi1, double psi2,
                        double psi12, double v,

                        /* Vector of results to be updated */
                        int nprod,
                        /* Result	*/
                        double *res, int *stop_path) {
  /* Variable declaration */
  char *err;
  int l;
  /* Get the event */
  TARN_EVENT *event = (TARN_EVENT *)func_parm;
  TARN_EVENT_SV *event_sv = (TARN_EVENT_SV *)event->ptrModel;

  /* Initialize return values */
  res[0] = 0.0;
  res[1] = 0.0;
  *stop_path = 0;

  /* see if we need to initialize anything
          if ( event->bIsFirst )
                  resetTargetNote();
  */
  /* Check that there is anything to do
          if ( ! event->do_dom )
                  return 0;
  */
  /* Reconstruct the discount factors */
  /* coupon pay */
  event->fixed_pay = exp(
      event->fixed_pay_logdff - event_sv->fixed_pay_gam1 * ft1 -
      event_sv->fixed_pay_gam2 * ft2 - event_sv->fixed_pay_gam1_2 * psi1 -
      event_sv->fixed_pay_gam2_2 * psi2 - event_sv->fixed_pay_gam12 * psi12);
  /* floater */
  event->fltr_start = exp(
      event->fltr_start_logdff - event_sv->fltr_start_gam1 * ft1 -
      event_sv->fltr_start_gam2 * ft2 - event_sv->fltr_start_gam1_2 * psi1 -
      event_sv->fltr_start_gam2_2 * psi2 - event_sv->fltr_start_gam12 * psi12);
  event->fltr_end =
      exp(event->fltr_end_logdff - event_sv->fltr_end_gam1 * ft1 -
          event_sv->fltr_end_gam2 * ft2 - event_sv->fltr_end_gam1_2 * psi1 -
          event_sv->fltr_end_gam2_2 * psi2 - event_sv->fltr_end_gam12 * psi12);

  /* Call the coupon calculator */
  if (err = TargetNoteCoupon(event, res, stop_path))
    return err;

  /* Check to see if it has knocked out */
  if (*stop_path)
    return 0;

  /* Reconstruct the start float discount factors */
  for (l = 0; l < event->num_fund_df; l++) {
    event->fund_start[l] =
        exp(event->fund_start_logdff[l] - event_sv->fund_start_gam1[l] * ft1 -
            event_sv->fund_start_gam2[l] * ft2 -
            event_sv->fund_start_gam1_2[l] * psi1 -
            event_sv->fund_start_gam2_2[l] * psi2 -
            event_sv->fund_start_gam12[l] * psi12);
    event->fund_end[l] = exp(
        event->fund_end_logdff[l] - event_sv->fund_end_gam1[l] * ft1 -
        event_sv->fund_end_gam2[l] * ft2 - event_sv->fund_end_gam1_2[l] * psi1 -
        event_sv->fund_end_gam2_2[l] * psi2 -
        event_sv->fund_end_gam12[l] * psi12);
  }

  /* Call the funding calculator */
  return TargetNoteFunding(event, res);
}

/* set up the event structures for pricing the TargetNote */
char *set_TARN_EVENTS_SV(TARN_Struct *tarn, TARN_AUX *aux,
                         TARN_MC_AUX *mc_aux) {
  /* Variable declarations */
  int i, j, index;
  char *err = 0;
  TARN_MC_AUX_SV *mc_aux_sv = (TARN_MC_AUX_SV *)mc_aux->ptrModel;
  TARN_EVENT *event = 0;
  TARN_EVENT_SV *event_sv = 0;

  /*	Allocate memory for the product structure */
  if (!(mc_aux->void_prm = (void **)calloc(mc_aux->nStps, sizeof(void *)))) {
    err = "Memory allocation error in TargetNote";
    goto FREE_RETURN;
  }

  /* Set a counter for the nunmber of event dates */
  j = mc_aux->nEvents - 1;

  /* loop over the total number of time steps */
  for (i = mc_aux->nStps - 1; i >= 0; i--) {
    /* calculate the date */
    mc_aux->date[i] =
        ((long)(tarn->market.lToday + mc_aux->time[i] * 365.0 + 1.0E-08)) * 1.0;

    /* Check that we are not at the last date and calculate the model parameters
     */
    if (i < mc_aux->nStps - 1) {
      index = Get_Index(mc_aux->time[i + 1], mc_aux_sv->model.dPWTime,
                        mc_aux_sv->model.iNbPWTime);
      mc_aux_sv->sigma[i + 1] = mc_aux_sv->model.dSigma[index];
      mc_aux_sv->alpha[i + 1] = mc_aux_sv->model.dAlpha[index];
      mc_aux_sv->rho[i + 1] = mc_aux_sv->model.dRho[index];
      mc_aux_sv->lameps[i + 1] = mc_aux_sv->model.dLambdaEps[index];
      mc_aux_sv->lvleps[i + 1] = mc_aux_sv->model.dLvlEps[index];
      if (mc_aux->tnModel == TN_2F_SV) {
        mc_aux_sv->dLGMalpha[i + 1] = mc_aux_sv->model.dLGMAlpha[index];
        mc_aux_sv->dLGMgamma[i + 1] = mc_aux_sv->model.dLGMGamma;
        mc_aux_sv->dLGMrho[i + 1] = mc_aux_sv->model.dLGMRho[index];
        mc_aux_sv->dRho2[i + 1] = mc_aux_sv->model.dRho2[index];
      }
    }

    /* check if this is an event */
    if (j >= 0 && fabs(mc_aux->time[i] - mc_aux->evt_tms[j]) < 1.0E-08) {

      /* Allocate memory for the local structure */
      if (err = alloc_TARN_EVENT(&event, mc_aux))
        goto FREE_RETURN;
      event_sv = (TARN_EVENT_SV *)event->ptrModel;

      /* set flags */
      //			event->bIsFirst = ! j;
      //			event->do_dom = 1;
      event->iIsFinal = (i == mc_aux->nStps - 1);

      /* DF(t        , T*) reconstruction */
      mc_aux_sv->logdff_star[j] =
          log(swp_f_df(mc_aux->evt_dts[j], mc_aux_sv->model.lTStarDate,
                       tarn->market.szYieldCurve));
      mc_aux_sv->gam1_star[j] =
          (1.0 - exp(mc_aux_sv->model.dLambdaX *
                     (mc_aux_sv->model.dTStar - mc_aux->evt_tms[j]))) /
          mc_aux_sv->model.dLambdaX;
      mc_aux_sv->gam1_2_star[j] =
          0.5 * mc_aux_sv->gam1_star[j] * mc_aux_sv->gam1_star[j];
      if (mc_aux->tnModel == TN_2F_SV) {
        mc_aux_sv->logdff_star[j] =
            log(swp_f_df(mc_aux->date[i], mc_aux_sv->model.lTStarDate,
                         tarn->market.szYieldCurve));
        mc_aux_sv->gam1_star[j] =
            (1.0 - exp(mc_aux_sv->model.dLambdaX *
                       (mc_aux_sv->model.dTStar - mc_aux->evt_tms[j]))) /
            mc_aux_sv->model.dLambdaX;
        mc_aux_sv->gam2_star[j] =
            (1.0 - exp(mc_aux_sv->model.dLambdaX2 *
                       (mc_aux_sv->model.dTStar - mc_aux->evt_tms[j]))) /
            mc_aux_sv->model.dLambdaX2;
        mc_aux_sv->gam12_star[j] =
            mc_aux_sv->gam1_star[j] * mc_aux_sv->gam2_star[j];
        mc_aux_sv->gam1_2_star[j] =
            0.5 * mc_aux_sv->gam1_star[j] * mc_aux_sv->gam1_star[j];
        mc_aux_sv->gam2_2_star[j] =
            0.5 * mc_aux_sv->gam2_star[j] * mc_aux_sv->gam2_star[j];
      }

      /* Fill the legs */
      if (err = fillTargetNoteFundingSV(j, tarn, aux, mc_aux, event))
        goto FREE_RETURN;
      if (err =
              fillTargetNoteCoupon(j + aux->i1stCpn, tarn, aux, mc_aux, event))
        goto FREE_RETURN;
      if (err = fillTargetNoteCouponSV(j, aux, mc_aux, event))
        goto FREE_RETURN;

      /* Calculate the numeraire values */

      mc_aux->is_event[i] = 1;
      mc_aux->void_prm[i] = (void *)event;
      j--;
    }
    /* if is not an event date */
    else {
      mc_aux->is_event[i] = 0;
      mc_aux->void_prm[i] = NULL;
    }
  }

FREE_RETURN:
  return err;
}

// -------------------------------------------------------------------------------------------------------------------
// //
//
// Wrapper routine for pricing a TargetNote in LGMSV
//
char *TargetNoteMC_SV(char *szUnd, TARN_Struct *tarn, TARN_AUX *aux,
                      double **prod_val) {
  /* -----------------------------------------------------------------------------------------------------------------
   */
  /* Variable declaration */

  /* Dates and DF */
  double df;
  //	int iNum_Fund_DF;
  //	double dFundPV;

  /* Target Note */
  TARN_MC_AUX MC_AUX, *mc_aux = &MC_AUX;
  TARN_MC_AUX_SV *mc_aux_sv;

  /* Counters */
  int i;
  Err err = NULL;
  int num_col = 2;
  /* Variable declaration */
  /* -----------------------------------------------------------------------------------------------------------------
   */

  /* set up the MC auxiliary structure */
  init_TARN_MC_AUX(mc_aux, szUnd);
  if (err = alloc_TARN_MC_AUX(tarn, aux, mc_aux))
    goto FREE_RETURN;

  /* Set the target note values */
  initTargetNote(tarn);

  /* set up the event structures */
  if (err = set_TARN_EVENTS_SV(tarn, aux, mc_aux))
    goto FREE_RETURN;
  mc_aux_sv = (TARN_MC_AUX_SV *)mc_aux->ptrModel;

  /* Call the MC */
  if (mc_aux->tnModel == TN_2F_SV) {
    const int nUseMemOpt = 1;

    if (nUseMemOpt) {

      if (err = lgmSV2F_mc_balsam_optim_mem /*_rev*/ (
              /*	Time Information  */
              mc_aux->nStps, mc_aux->nEvents, mc_aux->time, mc_aux->date,
              tarn->pricing.num_paths,
              /*	Model data Information	*/
              mc_aux_sv->model.dLambdaX, mc_aux_sv->model.dLambdaX2,

              mc_aux_sv->sigma, mc_aux_sv->dLGMalpha, mc_aux_sv->dLGMrho,

              mc_aux_sv->alpha, mc_aux_sv->lameps, mc_aux_sv->lvleps,
              mc_aux_sv->rho, mc_aux_sv->dRho2,

              /* Parameters for DF(t        ,T*) reconstruction */
              mc_aux_sv->logdff_star, mc_aux_sv->gam1_star,
              mc_aux_sv->gam2_star, mc_aux_sv->gam1_2_star,
              mc_aux_sv->gam2_2_star, mc_aux_sv->gam12_star,
              /* Parameters */
              &mc_aux_sv->Params, mc_aux->void_prm, mc_aux->is_event,
              /* for Optimisation of exercise boundary */
              0, NULL, NULL,
              /*	Initialisation function to be called at the beggining of
                 each path or NULL if none */
              resetTargetNote,
              /*	Payoff function */
              TargetNotePayoffSV2,
              /*	Result */
              num_col, prod_val))
        goto FREE_RETURN;
    } else {

      if (err = lgmSV2F_mc_balsam(
              /*	Time Information  */
              mc_aux->nStps, mc_aux->nEvents, mc_aux->time, mc_aux->date,
              tarn->pricing.num_paths,
              /*	Model data Information	*/
              mc_aux_sv->model.dLambdaX, mc_aux_sv->model.dLambdaX2,

              mc_aux_sv->sigma, mc_aux_sv->dLGMalpha, mc_aux_sv->dLGMrho,

              mc_aux_sv->alpha, mc_aux_sv->lameps, mc_aux_sv->lvleps,
              mc_aux_sv->rho, mc_aux_sv->dRho2,

              /* Parameters for DF(t        ,T*) reconstruction */
              mc_aux_sv->logdff_star, mc_aux_sv->gam1_star,
              mc_aux_sv->gam2_star, mc_aux_sv->gam1_2_star,
              mc_aux_sv->gam2_2_star, mc_aux_sv->gam12_star,
              /* Parameters */
              mc_aux_sv->Params, mc_aux->void_prm, mc_aux->is_event,
              /* for Optimisation of exercise boundary */
              0, NULL, NULL,
              /*	Initialisation function to be called at the beggining of
                 each path or NULL if none */
              resetTargetNote,
              /*	Payoff function */
              TargetNotePayoffSV2,
              /*	Result */
              num_col, prod_val))
        goto FREE_RETURN;
    }

  } else {
    if (err = lgmSV_mc_balsam(mc_aux->nStps, mc_aux->nEvents, mc_aux->time,
                              mc_aux->date, tarn->pricing.num_paths,
                              mc_aux_sv->model.dLambdaX,

                              mc_aux_sv->sigma,

                              mc_aux_sv->alpha, mc_aux_sv->lameps,
                              mc_aux_sv->lvleps, mc_aux_sv->rho,

                              mc_aux_sv->logdff_star, mc_aux_sv->gam1_star,
                              mc_aux_sv->gam1_2_star,

                              mc_aux_sv->Params,

                              mc_aux->void_prm, mc_aux->is_event, 0, NULL, NULL,

                              resetTargetNote,

                              TargetNotePayoffSV1, num_col, prod_val))
      goto FREE_RETURN;
  }

  /* calculate the initial value of the numeraire and hence discount */
  df = swp_f_df(tarn->market.lToday, mc_aux_sv->model.lTStarDate,
                tarn->market.szYieldCurve);
  for (i = 0; i < num_col; i++) {
    (prod_val)[i][0] *= df;
    (prod_val)[i][1] *= df;
  }

  /*	Add PV of Past
          iNum_Fund_DF = aux->ivFundEndIndex[0] - aux->ivFundStartIndex[0] + 1;
          dFundPV = 0.0;
          for ( i=0; i<iNum_Fund_DF; i++ )
          {
                  double dFund_Start_DF = swp_f_df(tarn->market.lToday        ,
     tarn->deal.lvFundStartDates[i]        , tarn->market.szYieldCurve); double
     dFund_End_DF = swp_f_df(tarn->market.lToday        ,
     tarn->deal.lvFundEndDates[i] , tarn->market.szYieldCurve); double
     dFund_cvgMargin = 1.0 - coverage( tarn->deal.lvFundStartDates[i]        ,
     tarn->deal.lvFundEndDates[i]        , tarn->deal.basisFund )
                                                                          * (
     tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i] ); dFundPV +=
     dFund_Start_DF - dFund_cvgMargin * dFund_End_DF;
          }
          (prod_val)[0][0] += dFundPV;*/
  prod_val[0][0] += aux->dHistFundPV;
  prod_val[1][0] += aux->dHistCpnPV;

FREE_RETURN:
  free_TARN_MC_AUX(mc_aux);

  return err;
}
