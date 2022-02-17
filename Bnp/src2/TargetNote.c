// ------------------------------------------------------------------------------------------------------------------
// //
//
// TargetNote.c
//
//
#include "math.h"

#include "BGMEval.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "LGM2FMC.h"
#include "LGM2Fgrfn.h"
#include "LGM2Fpde.h"
#include "SrtAccess.h"
#include "TargetNote.h"
#include "TargetNoteProdStruct.h"
#include "srt_h_all.h"

// Extern global variable
extern TargetNote TARGET_NOTE;

char *TARN_TimeSteps_2F(TARN_Struct *tarn, TARN_AUX *aux, TARN_MC_AUX *aux_mc) {
  // Variable declaration
  int iEvt, iCpn;

  /* calculate the number of events and time steps */
  aux_mc->nEvents = tarn->deal.nCouponDates - aux->i1stCpn;
  aux_mc->nStps = aux_mc->nEvents + !aux->iIsTodayFixing;

  /* allocate the memory for the dates  , times and event dates and times */
  if (!(aux_mc->date = (double *)calloc(aux_mc->nStps, sizeof(double))))
    return "Memory allocation error in TARN_TimeSteps_2F";
  if (!(aux_mc->time = (double *)calloc(aux_mc->nStps, sizeof(double))))
    return "Memory allocation error in TARN_TimeSteps_2F";
  /*
          if ( ! (aux_mc->evt_dts = (double*) calloc(aux_mc->nEvents  , sizeof
     (double)))  ) return "Memory allocation error in TARN_TimeSteps_2F"; if ( !
     (aux_mc->evt_tms = (double*) calloc(aux_mc->nEvents  , sizeof (double)))  )
                  return "Memory allocation error in TARN_TimeSteps_2F";
  */

  /* set the first event index */
  iCpn = aux->i1stCpn;
  iEvt = 0;

  /* if today is not an event date  , we need to add in the zero date */
  if (!aux->iIsTodayFixing) {
    aux_mc->time[iEvt] = 0.0;
    aux_mc->date[iEvt] = tarn->market.lToday; // today + DAYS_IN_YEAR * time[i];
    /*		aux_mc->evt_dts[iEvt] = tarn->market.lToday;	//today + DAYS_IN_YEAR *
       time[i]; aux_mc->evt_tms[iEvt] = 0.0; aux_mc->is_event[iEvt] = 0;
    */
    ++iEvt;
  }

  /* Loop over the coupon fixing dates setting the dates and times */
  while (iCpn < tarn->deal.nCouponDates) {
    aux_mc->time[iEvt] =
        (tarn->deal.lvFltrFixingDates[iCpn] - tarn->market.lToday) / 365.0;
    aux_mc->date[iEvt] =
        tarn->deal.lvFltrFixingDates[iCpn]; // today + DAYS_IN_YEAR * time[i];
    /*		aux_mc->evt_dts[iEvt] = tarn->deal.lvFltrFixingDates[iCpn];	//today +
       DAYS_IN_YEAR * time[i]; aux_mc->evt_tms[iEvt] = (
       tarn->deal.lvFltrFixingDates[iCpn] - tarn->market.lToday ) /365.0;
    */
    ++iEvt;
    ++iCpn;
  }

  return 0;
}

/* ------------------------------------------------------------------------------------------------------------------
//

Payoff function for the TargetNote

*/
char *TargetNotePayoff2F(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double R1D, double R2D,
    /* Results */
    int num_col,
    double *res, // res[0] = fund  , res[1] = coupon
    int *stop_path) {
  /* Variable declaration */
  TARN_EVENT *event;
  TARN_EVENT_2F *event_2f;

  int l;
  double dFloatPV, dFloater, dCoupon, dNotionalRepay, dNextCumulCoupon;
  char *err = 0;

  /* Initialize return values */
  res[0] = 0.0;
  res[1] = 0.0;
  *stop_path = 0;

  /* Get the event */
  event = (TARN_EVENT *)func_parm;
  event_2f = (TARN_EVENT_2F *)event->ptrModel;

  /* Check that there is anything to do
          if ( ! event->do_dom )
                  return 0;
  */
  /* Reconstruct the coupon pay discount factors */
  event->fixed_pay =
      exp(event->fixed_pay_logdff + event_2f->fixed_pay_gam12 -
          event_2f->fixed_pay_gam1 * R1D + -event_2f->fixed_pay_gam2 * R2D);

  /* Reconstruct the floater discount factors */
  event->fltr_start =
      exp(event->fltr_start_logdff + event_2f->fltr_start_gam12 -
          event_2f->fltr_start_gam1 * R1D + -event_2f->fltr_start_gam2 * R2D);
  event->fltr_end =
      exp(event->fltr_end_logdff + event_2f->fltr_end_gam12 -
          event_2f->fltr_end_gam1 * R1D + -event_2f->fltr_end_gam2 * R2D);

  /* Calculate the coupon */
  dFloater = (event->fltr_start / event->fltr_end - 1.0) / event->fltr_cvg +
             event->fltr_spread;
  dCoupon = event->fixed_coupon + event->fltr_gearing * dFloater;

  /* Check to see if we need to floor or cap it */
  if (event->iIsFloored && dCoupon < event->dFloor)
    dCoupon = event->dFloor;
  if (event->iIsCapped && dCoupon > event->dCap)
    dCoupon = event->dFloor;

  /* Calculate the notional to be repaid */
  dNextCumulCoupon = TARGET_NOTE.m_dCumulCoupon + dCoupon * event->cumul_cvg;
  dNotionalRepay =
      1.0 - (TARGET_NOTE.m_dTarget - dNextCumulCoupon) / TARGET_NOTE.m_dSpread;
  dNotionalRepay = dNotionalRepay > 1.0 ? 1.0 : dNotionalRepay;
  dNotionalRepay = dNotionalRepay < TARGET_NOTE.m_dCouponNotionalRepaid
                       ? TARGET_NOTE.m_dCouponNotionalRepaid
                       : dNotionalRepay;

  /* Call the relevant coupon calculation function */
  switch (TARGET_NOTE.m_couponType) {
  case TN_ZC:
    err = TargetNoteCoupon_ZC(event, dCoupon, dNotionalRepay, res, stop_path);
    break;
  case TN_SWAP:
    err = TargetNoteCoupon_SWAP(event, dCoupon, dNotionalRepay, res, stop_path);
    break;
  case TN_SWAPKO:
    err =
        TargetNoteCoupon_SWAPKO(event, dCoupon, dNotionalRepay, res, stop_path);
    break;
  case TN_FINAL:
    err = TargetNoteCoupon_FINAL(event, dCoupon, dFloater, dNotionalRepay, res,
                                 stop_path);
    break;
  default:
    return "Unknown coupon type";
  }
  if (err)
    return err;
  if (*stop_path)
    return 0;

  /* Add the coupon to the cumulative value */
  TARGET_NOTE.m_dCumulCoupon += dCoupon * event->cumul_cvg;

  /* Reconstruct the start float discount factors */
  for (l = 0; l < event->num_fund_df; l++) {
    event->fund_start[l] =
        exp(event->fund_start_logdff[l] + event_2f->fund_start_gam12[l] -
            event_2f->fund_start_gam1[l] * R1D +
            -event_2f->fund_start_gam2[l] * R2D);
    event->fund_end[l] = exp(
        event->fund_end_logdff[l] + event_2f->fund_end_gam12[l] -
        event_2f->fund_end_gam1[l] * R1D + -event_2f->fund_end_gam2[l] * R2D);
  }

  /* Calculate the funding */
  dFloatPV = 0.0;
  for (l = 0; l < event->num_fund_df; l++)
    dFloatPV +=
        event->fund_start[l] - event->fund_cvgMargin[l] * event->fund_end[l];
  res[0] = dFloatPV * TARGET_NOTE.m_dFundingNotional;

  /* Return */
  return 0;
}

/* ------------------------------------------------------------------------------------------------------------------
//

Payoff function for the TargetNote

*/
char *TargetNoteKO2F(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double R1D, double R2D,
    /* Results */
    int num_col,
    double *res, // res[i] = 0/1
    int *stop_path) {
  /* Variable declaration */
  TARN_EVENT *event;
  TARN_EVENT_2F *event_2f;

  double dFloater, dCoupon;
  char *err = 0;

  /* Get the event */
  event = (TARN_EVENT *)func_parm;
  event_2f = (TARN_EVENT_2F *)event->ptrModel;

  /* Initialize return values */
  memset(res, 0, sizeof(double) * event->nEventDates);
  /*	for ( l=0; l<event->nEventDates; l++ )
                  res[l] = 0.0;
  */
  *stop_path = 0;

  /* Reconstruct the coupon pay discount factors */
  event->fixed_pay =
      exp(event->fixed_pay_logdff + event_2f->fixed_pay_gam12 -
          event_2f->fixed_pay_gam1 * R1D + -event_2f->fixed_pay_gam2 * R2D);

  /* Reconstruct the floater discount factors */
  event->fltr_start =
      exp(event->fltr_start_logdff + event_2f->fltr_start_gam12 -
          event_2f->fltr_start_gam1 * R1D + -event_2f->fltr_start_gam2 * R2D);
  event->fltr_end =
      exp(event->fltr_end_logdff + event_2f->fltr_end_gam12 -
          event_2f->fltr_end_gam1 * R1D + -event_2f->fltr_end_gam2 * R2D);

  /* Calculate the coupon */
  dFloater = (event->fltr_start / event->fltr_end - 1.0) / event->fltr_cvg +
             event->fltr_spread;
  dCoupon = event->fixed_coupon + event->fltr_gearing * dFloater;

  /* Check to see if we need to floor or cap it */
  if (event->iIsFloored && dCoupon < event->dFloor)
    dCoupon = event->dFloor;
  if (event->iIsCapped && dCoupon > event->dCap)
    dCoupon = event->dFloor;

  /* Add the coupon to the cumulative value */
  TARGET_NOTE.m_dCumulCoupon += dCoupon * event->cumul_cvg;

  /* check if it has knocked out */
  if (TARGET_NOTE.m_dCumulCoupon >= TARGET_NOTE.m_dTarget) {

    res[event->iEventDate] =
        exp(event->numeraire_logdff + event_2f->numeraire_gam12 -
            event_2f->numeraire_gam1 * R1D + -event_2f->numeraire_gam2 * R2D);
    *stop_path = 1;
  }

  /* Return */
  return 0;
}

char *fillTargetNoteFunding2F(int iEvt, int iStp, TARN_Struct *tarn,
                              TARN_AUX *aux, TARN_MC_AUX *mc_aux,
                              TARN_EVENT *event) {
  /* Variable declaration */
  int k, iCpn;
  char *dom_yc;
  double lam_dom, lam_dom2, *dom_phi1, *dom_phi2, *dom_phi12, *date;
  long lEvent;
  long *lvFundEndIndex = aux->ivFundEndIndex,
       *lvFundStartIndex = aux->ivFundStartIndex;
  TARN_MC_AUX_2F *mkt2FPtr = (TARN_MC_AUX_2F *)mc_aux->ptrModel;
  TARN_EVENT_2F *event_2f = (TARN_EVENT_2F *)event->ptrModel;

  //	init_TARN_EVENT_2F( event_2f );

  /* Variable initializations */
  lam_dom = mkt2FPtr->lam;
  lam_dom2 = mkt2FPtr->lam2;
  dom_phi1 = mkt2FPtr->phi1;
  dom_phi2 = mkt2FPtr->phi2;
  dom_phi12 = mkt2FPtr->phi12;

  dom_yc = tarn->market.szYieldCurve;
  date = mc_aux->evt_dts;
  lEvent = tarn->deal.lvFltrFixingDates[iEvt];
  iCpn = iEvt + 1;

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

    event_2f->fund_start_gam1 = dvector(0, event->num_fund_df - 1);
    event_2f->fund_start_gam2 = dvector(0, event->num_fund_df - 1);
    event_2f->fund_start_gam12 = dvector(0, event->num_fund_df - 1);

    event->fund_end_tms = dvector(
        0, event->num_fund_df -
               1); // xStr.evt[i].evt->dft[fund_idx]; // discount factor times
    event->fund_end_dts = dvector(
        0, event->num_fund_df -
               1); // xStr.evt[i].evt->dfd[fund_idx]; // discount factor dates
    event->fund_end_logdff = dvector(0, event->num_fund_df - 1);
    event->fund_end = dvector(0, event->num_fund_df - 1);

    event_2f->fund_end_gam1 = dvector(0, event->num_fund_df - 1);
    event_2f->fund_end_gam2 = dvector(0, event->num_fund_df - 1);
    event_2f->fund_end_gam12 = dvector(0, event->num_fund_df - 1);

    event->fund_cvgMargin = dvector(0, event->num_fund_df - 1);

    /* Add error handling in the even of memory problems
                    if (!event->dom_dff || !event->dom_gam1 || !event->dom_gam2
       || !event->dom_gam12)
                    {
                            err = "Memory allocation error (8) in
       SrtGrfn5DFXMc"; goto FREE_RETURN;
                    }
    */

    for (k = 0; k < event->num_fund_df; k++) {
      int iFund = k + lvFundStartIndex[iCpn];
      /* start of the float period quantities */
      event->fund_start_dts[k] = tarn->deal.lvFundStartDates[iFund];
      event->fund_start_tms[k] = (event->fund_start_dts[k] - lEvent) / 365.0;
      event->fund_start_logdff[k] =
          log(swp_f_df(lEvent, event->fund_start_dts[k], (char *)dom_yc));

      event_2f->fund_start_gam1[k] =
          (1.0 - exp(-lam_dom * event->fund_start_tms[k])) / lam_dom;
      event_2f->fund_start_gam2[k] =
          (1.0 - exp(-lam_dom2 * event->fund_start_tms[k])) / lam_dom2;
      event_2f->fund_start_gam12[k] =
          -0.5 * (event_2f->fund_start_gam1[k] * event_2f->fund_start_gam1[k] *
                      dom_phi1[iStp] +
                  event_2f->fund_start_gam2[k] * event_2f->fund_start_gam2[k] *
                      dom_phi2[iStp]) -
          event_2f->fund_start_gam1[k] * event_2f->fund_start_gam2[k] *
              dom_phi12[iStp];
      /* end of the float period quantities */
      event->fund_end_dts[k] = tarn->deal.lvFundEndDates[iFund];
      event->fund_end_tms[k] = (event->fund_end_dts[k] - lEvent) / 365.0;
      event->fund_end_logdff[k] =
          log(swp_f_df(lEvent, event->fund_end_dts[k], (char *)dom_yc));
      event_2f->fund_end_gam1[k] =
          (1.0 - exp(-lam_dom * event->fund_end_tms[k])) / lam_dom;
      event_2f->fund_end_gam2[k] =
          (1.0 - exp(-lam_dom2 * event->fund_end_tms[k])) / lam_dom2;

      event_2f->fund_end_gam12[k] =
          -0.5 * (event_2f->fund_end_gam1[k] * event_2f->fund_end_gam1[k] *
                      dom_phi1[iStp] +
                  event_2f->fund_end_gam2[k] * event_2f->fund_end_gam2[k] *
                      dom_phi2[iStp]) -
          event_2f->fund_end_gam1[k] * event_2f->fund_end_gam2[k] *
              dom_phi12[iStp];
      /* period quantities */
      event->fund_cvgMargin[k] =
          1.0 - aux->dvFundCvg[iFund] * (tarn->deal.dvFundMargin[iFund] +
                                         tarn->deal.dvFundSpread[iFund]);
    }

    event->ptrModel = (void *)event_2f;
  }

  return 0;
}

/* --------------------------------------------------------------------------------------------------------------------

        Set up the coupon leg

*/
char *fillTargetNoteCoupon(int i, /* actual index */
                           TARN_Struct *tarn, TARN_AUX *aux,
                           TARN_MC_AUX *mc_aux, TARN_EVENT *event) {
  /* Variable declaration */
  long lEvent;
  char *dom_yc;
  double dDFdate;

  /* add the number of skipped periods to get the actual event */
  //	i += aux->i1stCpn;

  /* Variable initializations */
  dom_yc = tarn->market.szYieldCurve;
  lEvent = tarn->deal.lvFltrFixingDates[i];
  dDFdate = (double)lEvent;

  /* Set up the target note */

  /* payment details */
  event->fixed_pay_dts = tarn->deal.lvCpnPayDates[i];
  event->fixed_pay_tms = (event->fixed_pay_dts - lEvent) / 365.0;
  event->fixed_pay_logdff =
      log(swp_f_df(dDFdate, event->fixed_pay_dts, (char *)dom_yc));

  /* period quantities */
  event->fixed_cvg = aux->dvCouponCvg[i];
  event->fixed_coupon = tarn->deal.dvCoupon[i];
  event->cumul_cvg = aux->dvCumulCvg[i];

  /* floater details:  start date */
  event->fltr_start_dts = tarn->deal.lvFltrStartDates[i];
  event->fltr_start_tms = (event->fltr_start_dts - lEvent) / 365.0;
  event->fltr_start_logdff =
      log(swp_f_df(dDFdate, event->fltr_start_dts, (char *)dom_yc));

  /* floater details:  end date */
  event->fltr_end_dts = tarn->deal.lvFltrEndDates[i];
  event->fltr_end_tms = (event->fltr_end_dts - lEvent) / 365.0;
  event->fltr_end_logdff =
      log(swp_f_df(dDFdate, event->fltr_end_dts, (char *)dom_yc));
  event->fltr_spread = tarn->deal.dvFltrSpread[i];

  /* period quantities */
  event->fltr_cvg = coverage((long)event->fltr_start_dts,
                             (long)event->fltr_end_dts, tarn->deal.basisFltr);
  event->fltr_gearing = tarn->deal.dvGearing[i];

  // cap/floor details
  event->iIsCapped = tarn->deal.ivIsCapped[i];
  event->dCap = tarn->deal.dvCap[i];
  event->iIsFloored = tarn->deal.ivIsFloored[i];
  event->dFloor = tarn->deal.dvFloor[i];

  /* copy the fee */
  event->dKnockOutFee = aux->dvKnockOutFee[i];

  return 0;
}

char *fillTargetNoteCoupon2F(int iStp, /* time step  , not event index */
                             TARN_MC_AUX *mc_aux, TARN_EVENT *event) {
  /* Variable declaration */
  double lam_dom, lam_dom2, *dom_phi1, *dom_phi2, *dom_phi12;
  TARN_EVENT_2F *event_2f = (TARN_EVENT_2F *)event->ptrModel;
  TARN_MC_AUX_2F *mc_aux_2f = (TARN_MC_AUX_2F *)mc_aux->ptrModel;

  /* Variable initializations */
  lam_dom = mc_aux_2f->lam;
  lam_dom2 = mc_aux_2f->lam2;
  dom_phi1 = mc_aux_2f->phi1;
  dom_phi2 = mc_aux_2f->phi2;
  dom_phi12 = mc_aux_2f->phi12;

  /* Set up the target note */

  /* payment details */
  event_2f->fixed_pay_gam1 =
      (1.0 - exp(-lam_dom * event->fixed_pay_tms)) / lam_dom;
  event_2f->fixed_pay_gam2 =
      (1.0 - exp(-lam_dom2 * event->fixed_pay_tms)) / lam_dom2;

  event_2f->fixed_pay_gam12 =
      -0.5 * (event_2f->fixed_pay_gam1 * event_2f->fixed_pay_gam1 *
                  dom_phi1[iStp] +
              event_2f->fixed_pay_gam2 * event_2f->fixed_pay_gam2 *
                  dom_phi2[iStp]) -
      event_2f->fixed_pay_gam1 * event_2f->fixed_pay_gam2 * dom_phi12[iStp];

  /* floater details:  start date */
  event_2f->fltr_start_gam1 =
      (1.0 - exp(-lam_dom * event->fltr_start_tms)) / lam_dom;
  event_2f->fltr_start_gam2 =
      (1.0 - exp(-lam_dom2 * event->fltr_start_tms)) / lam_dom2;
  event_2f->fltr_start_gam12 =
      -0.5 * (event_2f->fltr_start_gam1 * event_2f->fltr_start_gam1 *
                  dom_phi1[iStp] +
              event_2f->fltr_start_gam2 * event_2f->fltr_start_gam2 *
                  dom_phi2[iStp]) -
      event_2f->fltr_start_gam1 * event_2f->fltr_start_gam2 * dom_phi12[iStp];
  /* floater details:  end date */
  event_2f->fltr_end_gam1 =
      (1.0 - exp(-lam_dom * event->fltr_end_tms)) / lam_dom;
  event_2f->fltr_end_gam2 =
      (1.0 - exp(-lam_dom2 * event->fltr_end_tms)) / lam_dom2;
  event_2f->fltr_end_gam12 =
      -0.5 *
          (event_2f->fltr_end_gam1 * event_2f->fltr_end_gam1 * dom_phi1[iStp] +
           event_2f->fltr_end_gam2 * event_2f->fltr_end_gam2 * dom_phi2[iStp]) -
      event_2f->fltr_end_gam1 * event_2f->fltr_end_gam2 * dom_phi12[iStp];

  return 0;
}

/* ------------------------------------------------------------------------------------------------------------------
 */
/* ------------------------------------------------------------------------------------------------------------------
 */
/* ------------------------------------------------------------------------------------------------------------------
 */
/* set up the event structures for pricing the TargetNote */
char *set_TARN_EVENTS_2F(TARN_Struct *tarn, TARN_AUX *aux,
                         TARN_MC_AUX *mc_aux) {
  /* Variable declaration */
  int iStp, i;
  char *err = 0;
  TARN_EVENT *event = 0;
  TARN_EVENT_2F *event_2f = 0;
  TARN_MC_AUX_2F *mc_aux_2f = (TARN_MC_AUX_2F *)mc_aux->ptrModel;
  double DF = 0.0;

  /*	Allocate memory for the product structure */
  if (!(mc_aux->void_prm = (void **)calloc(mc_aux->nStps, sizeof(void *)))) {
    err = "Memory allocation error in TargetNote";
    goto FREE_RETURN;
  }

  /* check to see if we have added an empty event date */
  if (!aux->iIsTodayFixing) {
    mc_aux->void_prm[0] = 0;
    iStp = 1;
  } else
    iStp = 0;

  // Set the period information
  for (i = aux->i1stCpn; i < tarn->deal.nCouponDates; i++) {
    /* allocate memory for the event structure */
    if (err = alloc_TARN_EVENT(&event, mc_aux))
      goto FREE_RETURN;

    /* set the number of dates */
    event->iEventDate = i;
    event->nEventDates = tarn->deal.nCouponDates;

    /* Fill the coupon leg with market info */
    if (err = fillTargetNoteCoupon(i, tarn, aux, mc_aux, event))
      goto FREE_RETURN;

    /* Fill model dependent info in the leg */
    event_2f = (TARN_EVENT_2F *)event->ptrModel;
    init_TARN_EVENT_2F(event_2f);
    if (err = fillTargetNoteCoupon2F(iStp, mc_aux, event))
      goto FREE_RETURN;
    if (err = fillTargetNoteFunding2F(i, iStp, tarn, aux, mc_aux, event))
      goto FREE_RETURN;

    /* Fill the numeraire info */
    event_2f = (TARN_EVENT_2F *)event->ptrModel;
    if (err = (*tarn->market.getDF)(tarn->market.szYieldCurve,
                                    tarn->deal.lvFltrFixingDates[i],
                                    mc_aux_2f->pay_date, &DF))
      goto FREE_RETURN;
    event->numeraire_logdff = log(DF);
    event_2f->numeraire_gam1 =
        (1.0 -
         exp(-mc_aux_2f->lam *
             (mc_aux_2f->pay_date - tarn->deal.lvFltrFixingDates[i]) / 365.0)) /
        mc_aux_2f->lam;
    event_2f->numeraire_gam2 =
        (1.0 -
         exp(-mc_aux_2f->lam2 *
             (mc_aux_2f->pay_date - tarn->deal.lvFltrFixingDates[i]) / 365.0)) /
        mc_aux_2f->lam2;
    event_2f->numeraire_gam12 =
        -0.5 * (event_2f->numeraire_gam1 * event_2f->numeraire_gam1 *
                    mc_aux_2f->phi1[iStp] +
                event_2f->numeraire_gam2 * event_2f->numeraire_gam2 *
                    mc_aux_2f->phi2[iStp]) -
        event_2f->numeraire_gam1 * event_2f->numeraire_gam2 *
            mc_aux_2f->phi12[iStp];

    /* check if we are at the end */
    event->iIsFinal = (i == tarn->deal.nCouponDates - 1);

    /* save the local data */
    mc_aux->void_prm[iStp++] = (void *)event;
  }

FREE_RETURN:
  return 0;
}

/* ------------------------------------------------------------------------------------------------------------------
 */
/* ------------------------------------------------------------------------------------------------------------------
 */
/* ------------------------------------------------------------------------------------------------------------------
 */
/* LGM2F MonteCarlo routine for pricing the TargetNote */

char *TargetNoteMC_2F(char *szUnd, TARN_Struct *tarn, TARN_AUX *aux,
                      double **prod_val) {
  /* -----------------------------------------------------------------------------------------------------------------
   */
  /* Variable declaration */

  /* Dates and DF */
  //	double dFund_Start_DF  , dFund_End_DF  , dFund_cvgMargin  , dFundPV;
  //	int iNum_Fund_DF;
  double df;

  /* Target Note */
  TARN_MC_AUX MC_AUX, *mc_aux = &MC_AUX;
  TARN_MC_AUX_2F *mc_aux_2f;
  TARGET_NOTE_MODEL tarn_model = TN_2F;

  /* Book-keeping */
  int i, num_col = 2;
  clock_t t1, t2;
  Err err = 0;

  /* End of variable declaration */
  /* -----------------------------------------------------------------------------------------------------------------
   */

  /* Start the clock */
  t1 = clock();

  /* set up the MC aux structure */
  init_TARN_MC_AUX(mc_aux, szUnd);
  if (err = alloc_TARN_MC_AUX(tarn, aux, mc_aux))
    goto FREE_RETURN;

  /* Set the target note values */
  initTargetNote(tarn);

  /* set up the event structures */
  if (err = set_TARN_EVENTS_2F(tarn, aux, mc_aux))
    goto FREE_RETURN;

  /* check the time */
  t2 = clock();
  smessage("Phase 1 -preprocessing  , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

  mc_aux_2f = (TARN_MC_AUX_2F *)mc_aux->ptrModel;
  err = mc_main_lgm2f(
      /*	Time data */
      tarn->pricing.num_paths, num_col, mc_aux->time, mc_aux->date,
      mc_aux->nStps, 0, mc_aux_2f->fwd1, mc_aux_2f->fwd2, mc_aux_2f->exp1,
      mc_aux_2f->exp2, mc_aux_2f->phi1, mc_aux_2f->phi2, mc_aux_2f->phi12,
      mc_aux_2f->gam1_fwd, mc_aux_2f->gam2_fwd, mc_aux_2f->bond_pay,
      mc_aux_2f->gam1_pay, mc_aux_2f->gam2_pay, mc_aux->covar, mc_aux->void_prm,
      /* Do PECS adjustment */
      tarn->pricing.do_pecs, 0, NULL, NULL, resetTargetNote,
      /*	Payoff function */
      TargetNotePayoff2F, /*	Result */
      prod_val);

  /* Calculate the initial value of the numeraire and discount the payoffs */
  df = swp_f_zr(tarn->market.lToday, mc_aux_2f->pay_date,
                tarn->market.szYieldCurve);
  df = exp(-df * mc_aux_2f->pay_time);
  for (i = 0; i < num_col; i++) {
    prod_val[i][0] *= df;
    prod_val[i][1] *= df;
  }

  /* Add PV of the first funding coupon
          iNum_Fund_DF = aux->ivFundEndIndex[aux->i1stCpn] -
     aux->ivFundStartIndex[aux->i1stCpn] + 1; dFundPV = 0.0; for ( i=0;
     i<iNum_Fund_DF; i++ )
          {
                  dFund_Start_DF = swp_f_df(tarn->market.lToday  ,
     tarn->deal.lvFundStartDates[i]  , tarn->market.szYieldCurve); dFund_End_DF
     = swp_f_df(tarn->market.lToday  , tarn->deal.lvFundEndDates[i]  ,
     tarn->market.szYieldCurve); dFund_cvgMargin = 1.0 - coverage(
     tarn->deal.lvFundStartDates[i]  , tarn->deal.lvFundEndDates[i]  ,
     tarn->deal.basisFund )
                                                                          * (
     tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i] ); dFundPV +=
     dFund_Start_DF - dFund_cvgMargin * dFund_End_DF;
          }
          prod_val[0][0] += dFundPV;
  */
  prod_val[0][0] += aux->dHistFundPV;
  prod_val[1][0] += aux->dHistCpnPV;

FREE_RETURN:
  free_TARN_MC_AUX(mc_aux);

  return err;
}

/* ------------------------------------------------------------------------------------------------------------------
 */
/* ------------------------------------------------------------------------------------------------------------------
 */
/* ------------------------------------------------------------------------------------------------------------------
 */
/* LGM2F MonteCarlo routine for calculating the TargetNote KO probabilities */

char *TargetNoteMC_KO_2F(char *szUnd, TARN_Struct *tarn, TARN_AUX *aux,
                         double **prod_val) {
  /* -----------------------------------------------------------------------------------------------------------------
   */
  /* Variable declaration */

  /* Target Note */
  TARN_MC_AUX MC_AUX, *mc_aux = &MC_AUX;
  TARN_MC_AUX_2F *mc_aux_2f;
  TARGET_NOTE_MODEL tarn_model = TN_2F;

  /* Book-keeping */
  int num_col = tarn->deal.nCouponDates;
  clock_t t1, t2;
  Err err = 0;

  /* End of variable declaration */
  /* -----------------------------------------------------------------------------------------------------------------
   */

  /* Start the clock */
  t1 = clock();

  /* set up the MC aux structure */
  init_TARN_MC_AUX(mc_aux, szUnd);
  if (err = alloc_TARN_MC_AUX(tarn, aux, mc_aux))
    goto FREE_RETURN;

  /* Set the target note values */
  initTargetNote(tarn);

  /* set up the event structures */
  if (err = set_TARN_EVENTS_2F(tarn, aux, mc_aux))
    goto FREE_RETURN;

  /* check the time */
  t2 = clock();
  smessage("Phase 1 -preprocessing  , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

  mc_aux_2f = (TARN_MC_AUX_2F *)mc_aux->ptrModel;
  err = mc_main_lgm2f(
      /*	Time data */
      tarn->calibration.num_prob_paths, num_col, mc_aux->time, mc_aux->date,
      mc_aux->nStps, 0, mc_aux_2f->fwd1, mc_aux_2f->fwd2, mc_aux_2f->exp1,
      mc_aux_2f->exp2, mc_aux_2f->phi1, mc_aux_2f->phi2, mc_aux_2f->phi12,
      mc_aux_2f->gam1_fwd, mc_aux_2f->gam2_fwd, mc_aux_2f->bond_pay,
      mc_aux_2f->gam1_pay, mc_aux_2f->gam2_pay, mc_aux->covar, mc_aux->void_prm,
      /* Do PECS adjustment */
      tarn->pricing.do_pecs, 0, NULL, NULL, resetTargetNote,
      /*	Payoff function */
      TargetNoteKO2F, /*	Result */
      prod_val);

FREE_RETURN:
  free_TARN_MC_AUX(mc_aux);

  return err;
}
