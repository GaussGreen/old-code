/* -------------------------------------------------------------------------------------------------------------------
 */
/*
CRIFUtil.c
Routines:
CRIF_calcHistory:  from the historical fixings  , calculate what has happened.
*/
#include "CRIFUtil.h"
#include "CRIFCaller.h"
#include "CRIFProdStruct.h"
#include "GenericMidatAutocal.h"
#include "GenericMidatCalib.h"
#include "LGMFwdVol.h"
#include "LGMSVMC.h"
#include "opfnctns.h"
#include "srt_h_allFx3F.h"
#include "srt_h_lgmtypes.h"
#include "srtaccess.h"

void crif_free_simularg(CRIF_SIMULARG sSimulArg) {
  if (sSimulArg->sModel) {
    free_LGMSV_model(sSimulArg->sModel);
    free(sSimulArg->sModel);
  }

  if (sSimulArg->sMCArg) {
    crif_free_mc_arg(sSimulArg->sMCArg);
    free(sSimulArg->sMCArg);
  }
}

Err crif_allocate_simularg(CRIF_SIMULARG sSimulArg) {
  sSimulArg->sModel = calloc(1, sizeof(LGMSV_model));
  sSimulArg->sMCArg = calloc(1, sizeof(crif_mcarg));

  if (!sSimulArg->sMCArg || !sSimulArg->sModel) {
    return ("Memory allocation faillure in crif_allocate_simularg");
  }

  return NULL;
}
Err crif_allocate_mc_arg(CRIF_MCARG sMCArg) {
  sMCArg->dDates = calloc(sMCArg->iNbTimes, sizeof(double));
  sMCArg->dSigma = calloc(sMCArg->iNbTimes, sizeof(double));
  sMCArg->dAlphaSV = calloc(sMCArg->iNbTimes, sizeof(double));
  sMCArg->dLambdaSV = calloc(sMCArg->iNbTimes, sizeof(double));
  sMCArg->dLevelSV = calloc(sMCArg->iNbTimes, sizeof(double));
  sMCArg->dRhoSV = calloc(sMCArg->iNbTimes, sizeof(double));
  sMCArg->dRho2SV = calloc(sMCArg->iNbTimes, sizeof(double));
  sMCArg->dLGMAlpha = calloc(sMCArg->iNbTimes, sizeof(double));
  sMCArg->dLGMRho = calloc(sMCArg->iNbTimes, sizeof(double));

  sMCArg->iEvalEvent = calloc(sMCArg->iNbTimes, sizeof(int));
  sMCArg->iIndexEvent = calloc(sMCArg->iNbTimes, sizeof(int));
  sMCArg->vPayoffParams = calloc(sMCArg->iNbTimes, sizeof(void *));

  sMCArg->sMCEBParams = calloc(1, sizeof(MCEBParams));
  sMCArg->iOptimise = calloc(sMCArg->iNbTimes, sizeof(double));
  sMCArg->dMarketValue = calloc(sMCArg->iNbTimes, sizeof(double));

  if (!sMCArg->dDates || !sMCArg->dSigma || !sMCArg->dAlphaSV ||
      !sMCArg->dLambdaSV || !sMCArg->dLevelSV || !sMCArg->dRhoSV ||
      !sMCArg->dRho2SV || !sMCArg->dLGMAlpha || !sMCArg->dLGMRho ||
      !sMCArg->iEvalEvent || !sMCArg->iIndexEvent || !sMCArg->vPayoffParams ||
      !sMCArg->sMCEBParams || !sMCArg->iOptimise || !sMCArg->dMarketValue) {
    return "Memory allocation faillure in crif_allocate_mc_arg";
  }

  return NULL;
}

static double crif_beta_func(long lStartDate, long lEndDate, double dLambda) {
  return ((1.0 - exp(-dLambda * (lEndDate - lStartDate) * YEARS_IN_DAY)) /
          dLambda);
}

Err crif_fill_eval_const_coupon(CMSCTS_MARKET sMarket, CRIF_DEAL sDeal,
                                LGMSV_MODEL sModel, int iIndexCoupon,
                                CRIF_COUPONARG sCouponArg) {
  Err err = NULL;
  CRIF_EXOTIC sExotic;
  CRIF_EXOTIC_AUX sExoticAux;
  SwapDP sFixedSwap;
  SrtDateList sFixedSwapDates;

  long lFixingDate, lStartDate, lTheoEndDate, lPayDate;
  double dLevel, dFloating, dCashSwap, dSwap;
  double *dAllDFDates = NULL, **AllSchedules = NULL;

  int iOldNbAllDF, iNbAllDF;
  int i, j;
  int iIndexDF;
  double dCoef;

  sExotic = sDeal->sExotic;
  sExoticAux = sExotic->sAux;

  sCouponArg->iIndexCoupon = iIndexCoupon;

  /* reconstruction of DF(t  ,Tpay) / DF(t  ,Tstar) */
  sCouponArg->tpay_tstar_alpha = swp_f_df(
      sModel->lTStarDate, sExotic->lPayDate[iIndexCoupon], sMarket->cYcName);
  sCouponArg->tpay_tstar_alpha = -log(sCouponArg->tpay_tstar_alpha);
  sCouponArg->tpay_tstar_beta = crif_beta_func(
      sModel->lTStarDate, sExotic->lPayDate[iIndexCoupon], sModel->dLambdaX);
  sCouponArg->tpay_tstar_gamma =
      0.5 * sCouponArg->tpay_tstar_beta * sCouponArg->tpay_tstar_beta;

  if (sModel->iOne2F == 2) {
    sCouponArg->tpay_tstar_beta_2 = crif_beta_func(
        sModel->lTStarDate, sExotic->lPayDate[iIndexCoupon], sModel->dLambdaX2);
    sCouponArg->tpay_tstar_gamma_2 =
        0.5 * sCouponArg->tpay_tstar_beta_2 * sCouponArg->tpay_tstar_beta_2;
    sCouponArg->tpay_tstar_gamma_12 =
        sCouponArg->tpay_tstar_beta * sCouponArg->tpay_tstar_beta_2;
  }

  /* Reconstruction of the DF for CMS calculation */
  iNbAllDF = 0;

  if (1.0 > 0) {
    /* Memory allocation */
    AllSchedules = calloc(1, sizeof(double *));

    if (!AllSchedules) {
      err = "Memory allocation faillure in crif_fill_eval_const_coupon";
      goto FREE_RETURN;
    }

    lFixingDate = sExotic->lPaidFixingDate[iIndexCoupon];
    lStartDate = add_unit(lFixingDate, sExotic->iFixingLag, SRT_BDAY,
                          MODIFIED_SUCCEEDING);

    if (fabs(sExotic->dPaidGearing[iIndexCoupon]) > 0) {
      /* Get the schedule */
      err = add_tenor(lStartDate, sExotic->cPaidTenor[iIndexCoupon],
                      NO_BUSDAY_CONVENTION, &lTheoEndDate);
      if (err)
        goto FREE_RETURN;

      err = swp_f_initSwapDP(lStartDate, lTheoEndDate, sExotic->cPaidFreq,
                             sExotic->cPaidBasis, &sFixedSwap);

      if (err)
        goto FREE_RETURN;

      sFixedSwapDates = SwapDP_to_DateList(&sFixedSwap, MODIFIED_SUCCEEDING);

      /* Calculate Level and Swap */
      dLevel = 0.0;

      for (j = 1; j < sFixedSwapDates.len; j++) {
        dLevel += coverage(sFixedSwapDates.date[j - 1], sFixedSwapDates.date[j],
                           sFixedSwap.basis_code) *
                  swp_f_df(sMarket->lToday, sFixedSwapDates.date[j],
                           sMarket->cYcName);
      }

      dFloating =
          swp_f_df(sMarket->lToday, sFixedSwapDates.date[0], sMarket->cYcName) -
          swp_f_df(sMarket->lToday,
                   sFixedSwapDates.date[sFixedSwapDates.len - 1],
                   sMarket->cYcName);

      dCashSwap = dFloating / dLevel;

      err = swp_f_ForwardRate(lStartDate, lTheoEndDate, sExotic->cPaidFreq,
                              sExotic->cPaidBasis, sMarket->cYcName,
                              sExotic->cPaidRef, &dSwap);

      if (err)
        goto FREE_RETURN;

      sCouponArg->dSpread = dSwap - dCashSwap;
      sCouponArg->iNbCMSDF = sFixedSwapDates.len;

      /* Save the dates */
      iOldNbAllDF = iNbAllDF;
      iNbAllDF += sFixedSwapDates.len;

      dAllDFDates = realloc(dAllDFDates, iNbAllDF * sizeof(double));
      AllSchedules[0] = calloc(sFixedSwapDates.len, sizeof(double));

      if (!dAllDFDates || !AllSchedules[0]) {
        err = "Memory allocation faillure in crif_fill_eval_const_fixing";
        goto FREE_RETURN;
      }

      for (j = 0; j < sFixedSwapDates.len; j++) {
        AllSchedules[0][j] = sFixedSwapDates.date[j];
        dAllDFDates[iOldNbAllDF + j] = sFixedSwapDates.date[j];

        if (j > 0) {
          sCouponArg->dCMSCoverage[j] =
              coverage(sFixedSwapDates.date[j - 1], sFixedSwapDates.date[j],
                       sFixedSwap.basis_code);
        }
      }

      swp_f_free_in_DateList(sFixedSwapDates);
    } else {
      AllSchedules[0] = NULL;
      sCouponArg->iNbCMSDF = 0;
    }

    /* Find the needed DF */
    num_f_sort_vector(iNbAllDF, dAllDFDates);
    num_f_unique_vector(&iNbAllDF, dAllDFDates);
    sCouponArg->iNbDF = iNbAllDF;

    /* Find the indexes for each CMS */
    for (i = 0; i < 1; i++) {
      iIndexDF = 0;

      for (j = 0; j < sCouponArg->iNbCMSDF; j++) {
        while (dAllDFDates[iIndexDF] < AllSchedules[i][j] - 0.5)
          iIndexDF++;
        sCouponArg->iIndexCMSDF[j] = iIndexDF;
      }
    }

    /* Reconstruction Formula */
    for (i = 0; i < sCouponArg->iNbDF; i++) {
      lPayDate = (long)(dAllDFDates[i] + 0.5);
      sCouponArg->tdf_tstar_alpha[i] =
          swp_f_df(sModel->lTStarDate, lPayDate, sMarket->cYcName);
      sCouponArg->tdf_tstar_alpha[i] = -log(sCouponArg->tdf_tstar_alpha[i]);
      sCouponArg->tdf_tstar_beta[i] =
          crif_beta_func(sModel->lTStarDate, lPayDate, sModel->dLambdaX);
      sCouponArg->tdf_tstar_gamma[i] =
          0.5 * sCouponArg->tdf_tstar_beta[i] * sCouponArg->tdf_tstar_beta[i];

      if (sModel->iOne2F == 2) {
        sCouponArg->tdf_tstar_beta_2[i] =
            crif_beta_func(sModel->lTStarDate, lPayDate, sModel->dLambdaX2);
        sCouponArg->tdf_tstar_gamma_2[i] = 0.5 *
                                           sCouponArg->tdf_tstar_beta_2[i] *
                                           sCouponArg->tdf_tstar_beta_2[i];
        sCouponArg->tdf_tstar_gamma_12[i] =
            sCouponArg->tdf_tstar_beta[i] * sCouponArg->tdf_tstar_beta_2[i];
      }
    }
  }

  /* Rescale the payoff coefficients if needed */
  if (sExotic->eCouponType[iIndexCoupon] == REGULAR_CRIF) {
    dCoef = sExotic->dCoverage[iIndexCoupon] * sExotic->dNot[iIndexCoupon];
    sCouponArg->dCoef = dCoef;
    sCouponArg->dAlpha = sExotic->dPrvCpngearing[iIndexCoupon];
    sCouponArg->dBeta = sExotic->dPaidMargin[iIndexCoupon];
    sCouponArg->dGamma = sExotic->dPaidGearing[iIndexCoupon];
    sCouponArg->dFloor = sExotic->dFloor[iIndexCoupon];
    sCouponArg->dCap = sExotic->dCap[iIndexCoupon];
  }

  /* Check for Floor and Cap */
  sCouponArg->iSavedCoupon = 0;

FREE_RETURN:

  if (dAllDFDates)
    free(dAllDFDates);
  if (AllSchedules) {
    for (i = 0; i < 1; i++) {
      if (AllSchedules[i])
        free(AllSchedules[i]);
    }

    free(AllSchedules);
  }

  return err;
}

Err crif_fill_eval_const_call(CMSCTS_MARKET sMarket, CRIF_DEAL sDeal,
                              LGMSV_MODEL sModel, int iIndexCall,
                              int iLastFundIdx, CRIF_CALLARG sCallArg) {
  Err err = NULL;
  CMSCTS_FUNDING sFunding;
  CMSCTS_FUND_AUX sFundingAux;
  CRIF_CALL sCall;
  CRIF_CALL_AUX sCallAux;
  int i;
  int iFundStart;
  long lPayDate;

  sFunding = sDeal->sFunding;
  sFundingAux = sFunding->sAux;
  sCall = sDeal->sCall;
  sCallAux = sCall->sAux;

  iFundStart = sCallAux->iStartFundIdx[iIndexCall];

  sCallArg->iIndexCall = iIndexCall;

  /* Reconstruction of DF(t  , Tfund) / DF(t  ,T*) */
  sCallArg->iNbFundCoupon = sFunding->iNbCpn - iFundStart;

  sCallArg->iNbFundCouponPartial = iLastFundIdx - iFundStart;

  for (i = 0; i <= sCallArg->iNbFundCoupon; i++) {
    if (i < sCallArg->iNbFundCoupon) {
      lPayDate = sFunding->lPayDate[i + iFundStart];
    } else {
      lPayDate = sFunding->lStartDate[iFundStart];
    }

    sCallArg->tpay_tstar_alpha[i] =
        swp_f_df(sModel->lTStarDate, lPayDate, sMarket->cYcName);
    sCallArg->tpay_tstar_alpha[i] = -log(sCallArg->tpay_tstar_alpha[i]);
    sCallArg->tpay_tstar_beta[i] =
        crif_beta_func(sModel->lTStarDate, lPayDate, sModel->dLambdaX);
    sCallArg->tpay_tstar_gamma[i] =
        0.5 * sCallArg->tpay_tstar_beta[i] * sCallArg->tpay_tstar_beta[i];

    if (sModel->iOne2F == 2) {
      sCallArg->tpay_tstar_beta_2[i] =
          crif_beta_func(sModel->lTStarDate, lPayDate, sModel->dLambdaX2);
      sCallArg->tpay_tstar_gamma_2[i] =
          0.5 * sCallArg->tpay_tstar_beta_2[i] * sCallArg->tpay_tstar_beta_2[i];
      sCallArg->tpay_tstar_gamma_12[i] =
          sCallArg->tpay_tstar_beta[i] * sCallArg->tpay_tstar_beta_2[i];
    }
  }

  /* reconstruction of DF(t  ,Tset) / DF(t  ,Tstar) */
  sCallArg->tset_tstar_alpha = swp_f_df(
      sModel->lTStarDate, sCall->lSettlDate[iIndexCall], sMarket->cYcName);
  sCallArg->tset_tstar_alpha = -log(sCallArg->tset_tstar_alpha);
  sCallArg->tset_tstar_beta = crif_beta_func(
      sModel->lTStarDate, sCall->lSettlDate[iIndexCall], sModel->dLambdaX);
  sCallArg->tset_tstar_gamma =
      0.5 * sCallArg->tset_tstar_beta * sCallArg->tset_tstar_beta;

  if (sModel->iOne2F == 2) {
    sCallArg->tset_tstar_beta_2 = crif_beta_func(
        sModel->lTStarDate, sCall->lSettlDate[iIndexCall], sModel->dLambdaX2);
    sCallArg->tset_tstar_gamma_2 =
        0.5 * sCallArg->tset_tstar_beta_2 * sCallArg->tset_tstar_beta_2;
    sCallArg->tset_tstar_gamma_12 =
        sCallArg->tset_tstar_beta * sCallArg->tset_tstar_beta_2;
  }

  sCallArg->mceb_sCouponArg = calloc(1, sizeof(crif_couponarg));

  err = crif_fill_eval_const_coupon(sMarket, sDeal, sModel,
                                    sCallAux->iStartCpnIdx[iIndexCall],
                                    sCallArg->mceb_sCouponArg);

  return err;
}

Err crif_fill_mc_arg(CMSCTS_MARKET sMarket, CRIF_DEAL sDeal,
                     CRIF_SIMULARG sSimulArg,
                     CRIF_PRICING_PARAMS sPricingParams) {
  Err err = NULL;

  CMSCTS_FUNDING sFunding;
  CMSCTS_FUND_AUX sFundingAux;
  CRIF_EXOTIC sExotic;
  CRIF_EXOTIC_AUX sExoticAux;
  CRIF_CALL sCall;
  CRIF_CALL_AUX sCallAux;
  CRIF_MCARG sMCArg;
  LGMSV_MODEL sModel;

  CRIF_PAYARG sPayArg;

  int i, j;
  double dAllTimes[10000], *dNewAllTimes = NULL;
  long lNbAllTimes;
  int iIndexFunding, iIndexCoupon, iIndexCouponForFixing, iIndexFixing,
      iIndexCall, iNbEvent;
  int iIndexTS, iLastFundIndex, iLastFundIndexCall;

  sFunding = sDeal->sFunding;
  sFundingAux = sFunding->sAux;
  sExotic = sDeal->sExotic;
  sExoticAux = sExotic->sAux;
  sCall = sDeal->sCall;
  sCallAux = sCall->sAux;
  sMCArg = sSimulArg->sMCArg;
  sModel = sSimulArg->sModel;

  /* First get all the event dates */
  dAllTimes[0] = 0.0;
  lNbAllTimes = 1;

  for (i = sExoticAux->iNumStartIndex; i <= sExoticAux->iEndIndex; i++) {
    if (sExotic->eCouponType[i] == REGULAR_CRIF &&
        sExotic->lPaidFixingDate[i] >= sMarket->lToday + sMarket->iEODFixFlag) {
      /* Add Fixing of Paid CMS */
      dAllTimes[lNbAllTimes] =
          (sExotic->lPaidFixingDate[i] - sMarket->lToday) * YEARS_IN_DAY;
      lNbAllTimes++;
    }
  }

  /* Add the Call Dates */
  if (!sCall->iIsExercised) {
    for (i = sCallAux->iStartIndex; i < sCall->iNbCall; i++) {
      dAllTimes[lNbAllTimes] =
          (sCall->lExeDate[i] - sMarket->lToday) * YEARS_IN_DAY;
      lNbAllTimes++;
    }
  }

  dNewAllTimes = calloc(lNbAllTimes, sizeof(double));
  if (!dNewAllTimes) {
    err = "Memory allocation faillure in crif_fill_mc_arg";
    goto FREE_RETURN;
  }

  memcpy(dNewAllTimes, dAllTimes, lNbAllTimes * sizeof(double));
  num_f_sort_vector(lNbAllTimes, dNewAllTimes);
  num_f_unique_vector(&lNbAllTimes, dNewAllTimes);

  /* Fill the vector with intermediary steps */
  err = fill_time_vector_max_time(lNbAllTimes, dNewAllTimes,
                                  sPricingParams->sSimulParams->dMCMaxTime,
                                  &sMCArg->iNbTimes, &sMCArg->dTimes);
  if (err)
    goto FREE_RETURN;

  /* All the memory allocation */
  err = crif_allocate_mc_arg(sMCArg);
  if (err)
    goto FREE_RETURN;

  /* Fill the informations */
  iIndexFunding = sFundingAux->iNumStartIndex;
  iIndexFixing = 0;

  if (!sCall->iIsExercised) {
    iIndexCall = sCallAux->iStartIndex;
  } else {
    iIndexCall = sCall->iNbCall;
  }

  iIndexCoupon = sExoticAux->iNumStartIndex;
  iIndexCouponForFixing = iIndexCoupon;

  /* Check for extra allocation of path dependent */
  for (i = iIndexCoupon; i < sExotic->iNbCpn; i++) {
    if (sExotic->eCouponType[i] == REGULAR_CRIF) {
      sMCArg->dPathInfos =
          calloc(sPricingParams->sSimulParams->lNbPaths + 1, sizeof(double));

      if (!sMCArg->dPathInfos) {
        err = "Memory allocation faillure in crif_fill_mc_arg";
        goto FREE_RETURN;
      }

      /* Initialisation */
      if (sExotic->lPaidFixingDate[i] <
          sMarket->lToday + sMarket->iEODFixFlag) {
        for (j = 0; j < sPricingParams->sSimulParams->lNbPaths + 1; j++) {
          sMCArg->dPathInfos[j] = sExotic->dPastPaidCpn[i];
        }
      }

      break;
    }
  }

  /* Check for Call Fees */
  sMCArg->iHasCallFee = 0;

  for (i = iIndexCall; i < sCall->iNbCall; i++) {
    if (fabs(sCall->dFee[i]) > 1.0E-08) {
      sMCArg->iHasCallFee = 1;
      break;
    }
  }

  /* Check for extra allocation for coupon */
  sMCArg->iHasFloorCap = 0;

  for (i = iIndexCoupon; i < sExotic->iNbCpn; i++) {
    if (sExotic->eCouponType[i] == REGULAR_CRIF) {
      sMCArg->dSavedCoupon =
          calloc(sPricingParams->sSimulParams->lNbPaths + 1, sizeof(double));
      sMCArg->dSavedDF =
          calloc(sPricingParams->sSimulParams->lNbPaths + 1, sizeof(double));
      sMCArg->iHasFloorCap = 1;

      if (!sMCArg->dSavedCoupon || !sMCArg->dSavedDF) {
        err = "Memory allocation faillure in crif_fill_mc_arg";
        goto FREE_RETURN;
      }

      break;
    }
  }

  iNbEvent = 0;

  for (i = 0; i < sMCArg->iNbTimes; i++) {
    /* Initialisation */
    sMCArg->dDates[i] =
        (long)(sMarket->lToday + sMCArg->dTimes[i] * DAYS_IN_YEAR + 0.5);
    sMCArg->vPayoffParams[i] = (void *)calloc(1, sizeof(crif_payarg));
    sMCArg->dMarketValue[i] = 0.0;
    sMCArg->iEvalEvent[i] = 0;
    sMCArg->iOptimise[i] = 0;

    sPayArg = (CRIF_PAYARG)(sMCArg->vPayoffParams[i]);

    sPayArg->iIsCall = 0;
    sPayArg->iIsCoupon = 0;
    sPayArg->iIsFixing = 0;
    sPayArg->sDeal = sDeal;
    sPayArg->sSimulArg = sSimulArg;
    sPayArg->iNumeraireType = sPricingParams->sSimulParams->iNumeraireType;
    sPayArg->dPathInfos = sMCArg->dPathInfos;
    sPayArg->dSavedCoupon = sMCArg->dSavedCoupon;
    sPayArg->dSavedDF = sMCArg->dSavedDF;

    /* Check for Call date */
    if (iIndexCall < sCall->iNbCall &&
        fabs(sMCArg->dDates[i] - sCall->lExeDate[iIndexCall]) < 0.5) {
      sPayArg->iIsCall = 1;
      sPayArg->sCallArg = calloc(1, sizeof(crif_callarg));

      if (!sPayArg->sCallArg) {
        err = "Memory allocation faillure in crif_fill_mc_arg";
        goto FREE_RETURN;
      }

      /* Look for next funding */
      if (iIndexCall < sCall->iNbCall - 1) {
        iLastFundIndexCall = sCallAux->iStartFundIdx[iIndexCall + 1];
      } else {
        iLastFundIndexCall = sFunding->iNbCpn;
      }

      iLastFundIndex = sCallAux->iStartFundIdx[iIndexCall] + 1;

      while (iLastFundIndex < iLastFundIndexCall &&
             sFundingAux->iFundingControll[iLastFundIndex] != 2) {
        iLastFundIndex++;
      }

      err = crif_fill_eval_const_call(sMarket, sDeal, sModel, iIndexCall,
                                      iLastFundIndex, sPayArg->sCallArg);

      if (err)
        goto FREE_RETURN;

      /* Update payoff with Funding */
      for (j = sCallAux->iStartFundIdx[iIndexCall]; j < iLastFundIndex; j++) {
        sMCArg->dMarketValue[i] +=
            -sDeal->dPayRec * sFundingAux->dMarketValue[j];
      }

      sMCArg->iEvalEvent[i] += 1;
      sMCArg->iIndexEvent[iNbEvent] = i;
      sMCArg->iOptimise[iNbEvent] = 1;
      iIndexCall++;
    }

    /* Check for Coupon date */
    if (iIndexCoupon < sExotic->iNbCpn &&
        (sExotic->eCouponType[iIndexCoupon] == REGULAR_CRIF) &&
        fabs(sMCArg->dDates[i] - sExotic->lPaidFixingDate[iIndexCoupon]) <
            0.5) {
      sPayArg->iIsCoupon = 1;

      sPayArg->sCouponArg = calloc(1, sizeof(crif_couponarg));

      if (!sPayArg->sCouponArg) {
        err = "Memory allocation faillure in crif_fill_mc_arg";
        goto FREE_RETURN;
      }

      err = crif_fill_eval_const_coupon(sMarket, sDeal, sModel, iIndexCoupon,
                                        sPayArg->sCouponArg);

      if (err)
        goto FREE_RETURN;

      sMCArg->iEvalEvent[i] = 1;
      sMCArg->iIndexEvent[iNbEvent] = i;
      sMCArg->iOptimise[iNbEvent] = 0;

      if (sExotic->eCouponType[iIndexCoupon] == REGULAR_CRIF) {
        sMCArg->dMarketValue[i] +=
            sDeal->dPayRec * sExoticAux->dMarketValue[iIndexCoupon];

        iIndexCouponForFixing++;
      }

      iIndexCoupon++;
    }

    if (sMCArg->iEvalEvent[i]) {
      iNbEvent++;
    }

    /* Fill diffusion params */
    iIndexTS = Get_Index(sMCArg->dTimes[i], sModel->dPWTime, sModel->iNbPWTime);

    sMCArg->dSigma[i] = sModel->dSigma[iIndexTS];
    sMCArg->dAlphaSV[i] = sModel->dAlpha[iIndexTS];
    sMCArg->dLambdaSV[i] = sModel->dLambdaEps[iIndexTS];
    sMCArg->dLevelSV[i] = sModel->dLvlEps[iIndexTS];
    sMCArg->dRhoSV[i] = sModel->dRho[iIndexTS];

    if (sModel->iOne2F == 2) {
      sMCArg->dRho2SV[i] = sModel->dRho2[iIndexTS];
      sMCArg->dLGMAlpha[i] = sModel->dLGMAlpha[iIndexTS];
      sMCArg->dLGMRho[i] = sModel->dLGMRho[iIndexTS];
    }
  }

  sMCArg->iNbEvent = iNbEvent;

FREE_RETURN:

  if (dNewAllTimes)
    free(dNewAllTimes);

  return err;
}

Err crif_payoff_1F(/* Event */
                   long path_index, double evt_date, double evt_time,
                   void *func_parm, double ft, double psi, double v, int nprod,
                   /* Vector of results to be updated */
                   double *prod_val, int *stop_path) {
  Err err = NULL;
  CRIF_DEAL sDeal;
  CMSCTS_FUNDING sFunding;
  CMSCTS_FUND_AUX sFundingAux;
  CRIF_EXOTIC sExotic;
  CRIF_EXOTIC_AUX sExoticAux;
  CRIF_CALL sCall;
  CRIF_CALL_AUX sCallAux;
  CRIF_SIMULARG sSimulArg;
  CRIF_MCARG sMCArg;
  LGMSV_MODEL sModel;

  CRIF_PAYARG sPayArg;
  CRIF_COUPONARG sCouponArg;
  CRIF_CALLARG sCallArg;

  double dFunding, dPartialFunding;
  double dFloat, dLevel, dCMS, dBasket, dCoupon;
  int i, j;
  int iStartFundIdx;
  double dFee;
  double DFR_pay, DFR_set, DFR_last;

  sPayArg = (CRIF_PAYARG)(func_parm);
  sDeal = sPayArg->sDeal;
  sSimulArg = sPayArg->sSimulArg;
  sCouponArg = sPayArg->sCouponArg;
  sCallArg = sPayArg->sCallArg;
  sMCArg = sSimulArg->sMCArg;
  sModel = sSimulArg->sModel;

  sFunding = sDeal->sFunding;
  sFundingAux = sFunding->sAux;
  sExotic = sDeal->sExotic;
  sExoticAux = sExotic->sAux;
  sCall = sDeal->sCall;
  sCallAux = sCall->sAux;

  memset(prod_val, 0, nprod * sizeof(double));

  if (sPayArg->iIsCall) {
    /* Value the Call */
    /* First evaluate funding and partial funding */

    dFunding = 0.0;
    dPartialFunding = 0.0;
    iStartFundIdx = sCallAux->iStartFundIdx[sCallArg->iIndexCall];

    for (i = 0; i < sCallArg->iNbFundCoupon; i++) {
      DFR_pay = exp(-sCallArg->tpay_tstar_alpha[i] -
                    sCallArg->tpay_tstar_beta[i] * ft -
                    sCallArg->tpay_tstar_gamma[i] * psi);

      dFunding += sFundingAux->dCouponPlusEx[iStartFundIdx + i] * DFR_pay;

      if (i < sCallArg->iNbFundCouponPartial) {
        if (i < sCallArg->iNbFundCouponPartial - 1) {
          dPartialFunding +=
              sFundingAux->dCouponPlusEx[iStartFundIdx + i] * DFR_pay;
        } else {
          dPartialFunding +=
              sFundingAux->dCouponPlusExPartial[iStartFundIdx + i] * DFR_pay;

          /* save last df */
          DFR_last = DFR_pay;
        }
      }
    }

    /* first coupon */
    DFR_pay =
        exp(-sCallArg->tpay_tstar_alpha[i] - sCallArg->tpay_tstar_beta[i] * ft -
            sCallArg->tpay_tstar_gamma[i] * psi);
    dFunding += sFunding->dNot[iStartFundIdx] * DFR_pay;
    dPartialFunding += sFunding->dNot[iStartFundIdx] * DFR_pay;

    prod_val[0] += sCall->dPayRec * dPartialFunding;
    prod_val[1] = prod_val[0];
    prod_val[2] = sCall->dPayRec * dFunding;
    prod_val[3] = v;

    // computes the relevant coupon datas for the next coupon //
    if (sCallArg->mceb_sCouponArg) {
      /* First evaluate all the DF */
      for (i = 0; i < sCallArg->mceb_sCouponArg->iNbDF; i++) {
        sCallArg->mceb_sCouponArg->dAllDF[i] =
            exp(-sCallArg->mceb_sCouponArg->tdf_tstar_alpha[i] -
                sCallArg->mceb_sCouponArg->tdf_tstar_beta[i] * ft -
                sCallArg->mceb_sCouponArg->tdf_tstar_gamma[i] * psi);
      }

      dBasket = 0.0;

      /* Reconstruct */
      for (i = 0; i < 1; i++) {
        if (sCallArg->mceb_sCouponArg->iNbCMSDF > 0) {
          dLevel = 0.0;

          for (j = 1; j < sCallArg->mceb_sCouponArg->iNbCMSDF; j++) {
            dLevel += sCallArg->mceb_sCouponArg->dCMSCoverage[j] *
                      sCallArg->mceb_sCouponArg
                          ->dAllDF[sCallArg->mceb_sCouponArg->iIndexCMSDF[j]];
          }

          dFloat = sCallArg->mceb_sCouponArg
                       ->dAllDF[sCallArg->mceb_sCouponArg->iIndexCMSDF[0]] -
                   sCallArg->mceb_sCouponArg
                       ->dAllDF[sCallArg->mceb_sCouponArg->iIndexCMSDF
                                    [sCallArg->mceb_sCouponArg->iNbCMSDF - 1]];
          dCMS = dFloat / dLevel + sCallArg->mceb_sCouponArg->dSpread;
          dBasket += dCMS;
        }
      }

      switch (sExotic->eCouponType[sCallArg->mceb_sCouponArg->iIndexCoupon]) {
      case REGULAR_CRIF: {
        dCoupon = sCallArg->mceb_sCouponArg->dAlpha *
                      sMCArg->dSavedCoupon[path_index] +
                  sCallArg->mceb_sCouponArg->dGamma * dBasket +
                  sCallArg->mceb_sCouponArg->dBeta;

        dCoupon = DMIN(DMAX(dCoupon, sCallArg->mceb_sCouponArg->dFloor),
                       sCallArg->mceb_sCouponArg->dCap);
        DFR_pay = exp(-sCallArg->mceb_sCouponArg->tpay_tstar_alpha -
                      sCallArg->mceb_sCouponArg->tpay_tstar_beta * ft -
                      sCallArg->mceb_sCouponArg->tpay_tstar_gamma * psi);
        break;
      }
      }

      // variables de regression
      if (fabs(sCallArg->mceb_sCouponArg->dGamma) > 1e-4) {
        prod_val[4] = dCoupon * DFR_pay *
                      sCallArg->mceb_sCouponArg->dCoef; // 1: le next coupon
        prod_val[5] = dCMS;                             // 2: le CMS
      } else
        prod_val[5] = dCMS; // 1: le CMS
    }

    /* then calculate the fee */
    if (sMCArg->iHasCallFee) {

      DFR_set =
          exp(-sCallArg->tset_tstar_alpha - sCallArg->tset_tstar_beta * ft -
              sCallArg->tset_tstar_gamma * psi);
      dFee = DFR_set * sCall->dFee[sCallArg->iIndexCall];
      prod_val[6] = dFee;
    }
  }

  if (sPayArg->iIsCoupon) {
    /* First evaluate all the DF */
    for (i = 0; i < sCouponArg->iNbDF; i++) {
      sCouponArg->dAllDF[i] = exp(-sCouponArg->tdf_tstar_alpha[i] -
                                  sCouponArg->tdf_tstar_beta[i] * ft -
                                  sCouponArg->tdf_tstar_gamma[i] * psi);
    }

    dBasket = 0.0;

    /* Reconstruct */
    for (i = 0; i < 1; i++) {
      if (sCouponArg->iNbCMSDF > 0) {
        dLevel = 0.0;

        for (j = 1; j < sCouponArg->iNbCMSDF; j++) {
          dLevel += sCouponArg->dCMSCoverage[j] *
                    sCouponArg->dAllDF[sCouponArg->iIndexCMSDF[j]];
        }

        dFloat =
            sCouponArg->dAllDF[sCouponArg->iIndexCMSDF[0]] -
            sCouponArg
                ->dAllDF[sCouponArg->iIndexCMSDF[sCouponArg->iNbCMSDF - 1]];
        dCMS = dFloat / dLevel + sCouponArg->dSpread;
        dBasket += dCMS;
      }
    }

    switch (sExotic->eCouponType[sCouponArg->iIndexCoupon]) {
    case REGULAR_CRIF: {
      dCoupon = sCouponArg->dAlpha * sMCArg->dSavedCoupon[path_index] +
                sCouponArg->dGamma * dBasket + sCouponArg->dBeta;

      dCoupon = DMIN(DMAX(dCoupon, sCouponArg->dFloor), sCouponArg->dCap);
      sMCArg->dSavedCoupon[path_index] = dCoupon;
      DFR_pay =
          exp(-sCouponArg->tpay_tstar_alpha - sCouponArg->tpay_tstar_beta * ft -
              sCouponArg->tpay_tstar_gamma * psi);

      prod_val[0] += -sCall->dPayRec * dCoupon * DFR_pay * sCouponArg->dCoef;

      break;
    }
    }
  }

  return NULL;
}

Err crif_payoff_2F(
    /* Event */
    long path_index, double evt_date, double evt_time, void *func_parm,
    double ft1, double ft2, double psi1, double psi2, double psi12, double v,
    int nprod,
    /* Vector of results to be updated */
    double *prod_val, int *stop_path) {
  Err err = NULL;
  CRIF_DEAL sDeal;
  CMSCTS_FUNDING sFunding;
  CMSCTS_FUND_AUX sFundingAux;
  CRIF_EXOTIC sExotic;
  CRIF_EXOTIC_AUX sExoticAux;
  CRIF_CALL sCall;
  CRIF_CALL_AUX sCallAux;
  CRIF_SIMULARG sSimulArg;
  CRIF_MCARG sMCArg;
  LGMSV_MODEL sModel;

  CRIF_PAYARG sPayArg;
  CRIF_COUPONARG sCouponArg;
  CRIF_CALLARG sCallArg;

  double dFunding, dPartialFunding;
  double dFloat, dLevel, dCMS, dBasket, dCoupon;
  int i, j;
  int iStartFundIdx;
  double dFee;
  double DFR_pay, DFR_set, DFR_last;

  sPayArg = (CRIF_PAYARG)(func_parm);
  sDeal = sPayArg->sDeal;
  sSimulArg = sPayArg->sSimulArg;
  sCouponArg = sPayArg->sCouponArg;
  sCallArg = sPayArg->sCallArg;
  sMCArg = sSimulArg->sMCArg;
  sModel = sSimulArg->sModel;

  sFunding = sDeal->sFunding;
  sFundingAux = sFunding->sAux;
  sExotic = sDeal->sExotic;
  sExoticAux = sExotic->sAux;
  sCall = sDeal->sCall;
  sCallAux = sCall->sAux;

  memset(prod_val, 0, nprod * sizeof(double));

  if (sPayArg->iIsCall) {
    /* Value the Call */
    /* First evaluate funding and partial funding */

    dFunding = 0.0;
    dPartialFunding = 0.0;
    iStartFundIdx = sCallAux->iStartFundIdx[sCallArg->iIndexCall];

    for (i = 0; i < sCallArg->iNbFundCoupon; i++) {
      DFR_pay = exp(-sCallArg->tpay_tstar_alpha[i] -
                    sCallArg->tpay_tstar_beta[i] * ft1 -
                    sCallArg->tpay_tstar_beta_2[i] * ft2 -
                    sCallArg->tpay_tstar_gamma[i] * psi1 -
                    sCallArg->tpay_tstar_gamma_2[i] * psi2 -
                    sCallArg->tpay_tstar_gamma_12[i] * psi12);

      dFunding += sFundingAux->dCouponPlusEx[iStartFundIdx + i] * DFR_pay;

      if (i < sCallArg->iNbFundCouponPartial) {
        if (i < sCallArg->iNbFundCouponPartial - 1) {
          dPartialFunding +=
              sFundingAux->dCouponPlusEx[iStartFundIdx + i] * DFR_pay;
        } else {
          dPartialFunding +=
              sFundingAux->dCouponPlusExPartial[iStartFundIdx + i] * DFR_pay;

          /* save last df */
          DFR_last = DFR_pay;
        }
      }
    }

    /* first coupon */
    DFR_pay = exp(-sCallArg->tpay_tstar_alpha[i] -
                  sCallArg->tpay_tstar_beta[i] * ft1 -
                  sCallArg->tpay_tstar_beta_2[i] * ft2 -
                  sCallArg->tpay_tstar_gamma[i] * psi1 -
                  sCallArg->tpay_tstar_gamma_2[i] * psi2 -
                  sCallArg->tpay_tstar_gamma_12[i] * psi12);

    dFunding += sFunding->dNot[iStartFundIdx] * DFR_pay;
    dPartialFunding += sFunding->dNot[iStartFundIdx] * DFR_pay;

    prod_val[0] += sCall->dPayRec * dPartialFunding;
    prod_val[1] = prod_val[0];
    prod_val[2] = sCall->dPayRec * dFunding; // le funding
    prod_val[3] = v;                         // la vol

    // computes the relevant coupon datas for the next coupon //
    if ((sCallArg->mceb_sCouponArg)) {
      /* First evaluate all the DF */
      for (i = 0; i < sCallArg->mceb_sCouponArg->iNbDF; i++) {
        sCallArg->mceb_sCouponArg->dAllDF[i] =
            exp(-sCallArg->mceb_sCouponArg->tdf_tstar_alpha[i] -
                sCallArg->mceb_sCouponArg->tdf_tstar_beta[i] * ft1 -
                sCallArg->mceb_sCouponArg->tdf_tstar_beta_2[i] * ft2 -
                sCallArg->mceb_sCouponArg->tdf_tstar_gamma[i] * psi1 -
                sCallArg->mceb_sCouponArg->tdf_tstar_gamma_2[i] * psi2 -
                sCallArg->mceb_sCouponArg->tdf_tstar_gamma_12[i] * psi12);
      }

      dBasket = 0.0;

      /* Reconstruct */
      for (i = 0; i < 1; i++) {
        if (sCallArg->mceb_sCouponArg->iNbCMSDF > 0) {
          dLevel = 0.0;

          for (j = 1; j < sCallArg->mceb_sCouponArg->iNbCMSDF; j++) {
            dLevel += sCallArg->mceb_sCouponArg->dCMSCoverage[j] *
                      sCallArg->mceb_sCouponArg
                          ->dAllDF[sCallArg->mceb_sCouponArg->iIndexCMSDF[j]];
          }

          dFloat = sCallArg->mceb_sCouponArg
                       ->dAllDF[sCallArg->mceb_sCouponArg->iIndexCMSDF[0]] -
                   sCallArg->mceb_sCouponArg
                       ->dAllDF[sCallArg->mceb_sCouponArg->iIndexCMSDF
                                    [sCallArg->mceb_sCouponArg->iNbCMSDF - 1]];
          dCMS = dFloat / dLevel + sCallArg->mceb_sCouponArg->dSpread;
          dBasket += dCMS;
        }
      }

      switch (sExotic->eCouponType[sCallArg->mceb_sCouponArg->iIndexCoupon]) {
      case REGULAR_CRIF: {
        dCoupon = sCallArg->mceb_sCouponArg->dAlpha *
                      sMCArg->dSavedCoupon[path_index] +
                  sCallArg->mceb_sCouponArg->dGamma * dBasket +
                  sCallArg->mceb_sCouponArg->dBeta;

        dCoupon = DMIN(DMAX(dCoupon, sCallArg->mceb_sCouponArg->dFloor),
                       sCallArg->mceb_sCouponArg->dCap);
        DFR_pay = exp(-sCallArg->mceb_sCouponArg->tpay_tstar_alpha -
                      sCallArg->mceb_sCouponArg->tpay_tstar_beta * ft1 -
                      sCallArg->mceb_sCouponArg->tpay_tstar_beta_2 * ft2 -
                      sCallArg->mceb_sCouponArg->tpay_tstar_gamma * psi1 -
                      sCallArg->mceb_sCouponArg->tpay_tstar_gamma_2 * psi2 -
                      sCallArg->mceb_sCouponArg->tpay_tstar_gamma_12 * psi12);

        break;
      }
      }

      // variables de regression
      prod_val[4] = dCoupon; // 1: le next coupon
      prod_val[5] = dCMS;    // 2: le CMS
      prod_val[6] =
          (sCall->dPayRec * dPartialFunding -
           dCoupon * DFR_pay * sCallArg->mceb_sCouponArg->dCoef) *
          (sCall->dPayRec * dPartialFunding -
           dCoupon * DFR_pay *
               sCallArg->mceb_sCouponArg->dCoef); // 3: le current coupon

      /* then calculate the fee */
      if (sMCArg->iHasCallFee) {
        DFR_set =
            exp(-sCallArg->tset_tstar_alpha - sCallArg->tset_tstar_beta * ft1 -
                sCallArg->tset_tstar_beta_2 * ft2 -
                sCallArg->tset_tstar_gamma * psi1 -
                sCallArg->tset_tstar_gamma_2 * psi2 -
                sCallArg->tset_tstar_gamma_12 * psi12);

        dFee = DFR_set * sCall->dFee[sCallArg->iIndexCall];
        prod_val[7] = dFee;
      }
    }

    /*		if (sCallArg->iNbFundCoupon > sCallArg->iNbFundCouponPartial)
                    {*/
    /*}
    else
    {
            prod_val[4] = 0.0;
    }*/
  }

  if (sPayArg->iIsCoupon) {
    /* First evaluate all the DF */
    for (i = 0; i < sCouponArg->iNbDF; i++) {
      sCouponArg->dAllDF[i] = exp(-sCouponArg->tdf_tstar_alpha[i] -
                                  sCouponArg->tdf_tstar_beta[i] * ft1 -
                                  sCouponArg->tdf_tstar_beta_2[i] * ft2 -
                                  sCouponArg->tdf_tstar_gamma[i] * psi1 -
                                  sCouponArg->tdf_tstar_gamma_2[i] * psi2 -
                                  sCouponArg->tdf_tstar_gamma_12[i] * psi12);
    }

    dBasket = 0.0;

    /* Reconstruct */
    for (i = 0; i < 1; i++) {
      if (sCouponArg->iNbCMSDF > 0) {
        dLevel = 0.0;

        for (j = 1; j < sCouponArg->iNbCMSDF; j++) {
          dLevel += sCouponArg->dCMSCoverage[j] *
                    sCouponArg->dAllDF[sCouponArg->iIndexCMSDF[j]];
        }

        dFloat =
            sCouponArg->dAllDF[sCouponArg->iIndexCMSDF[0]] -
            sCouponArg
                ->dAllDF[sCouponArg->iIndexCMSDF[sCouponArg->iNbCMSDF - 1]];
        dCMS = dFloat / dLevel + sCouponArg->dSpread;
        dBasket += dCMS;
      }
    }

    switch (sExotic->eCouponType[sCouponArg->iIndexCoupon]) {
    case REGULAR_CRIF: {
      dCoupon = sCouponArg->dAlpha * sMCArg->dSavedCoupon[path_index] +
                sCouponArg->dGamma * dBasket + sCouponArg->dBeta;

      dCoupon = DMIN(DMAX(dCoupon, sCouponArg->dFloor), sCouponArg->dCap);
      sMCArg->dSavedCoupon[path_index] = dCoupon;
      DFR_pay = exp(-sCouponArg->tpay_tstar_alpha -
                    sCouponArg->tpay_tstar_beta * ft1 -
                    sCouponArg->tpay_tstar_beta_2 * ft2 -
                    sCouponArg->tpay_tstar_gamma * psi1 -
                    sCouponArg->tpay_tstar_gamma_2 * psi2 -
                    sCouponArg->tpay_tstar_gamma_12 * psi12);

      prod_val[0] += -sCall->dPayRec * dCoupon * DFR_pay * sCouponArg->dCoef;

      break;
    }
    }
  }

  return NULL;
}

Err crif_adjust_payoff_2F(long lEvtIndex, long lNbPaths, long lNbProd,
                          double ***dSavedValues, void *vFuncParam,
                          MCEBPARAMS sMCEBParams) {
  return NULL;
}

Err crif_calculate_option(CMSCTS_MARKET sMarket, CRIF_DEAL sDeal,
                          CRIF_SIMULARG sSimulArg,
                          CRIF_PRICING_PARAMS sPricingParams,
                          CRIF_OUTPUTS sOutputs) {
  Err err = NULL;
  CMSCTS_FUNDING sFunding;
  CMSCTS_FUND_AUX sFundingAux;
  CRIF_EXOTIC sExotic;
  CRIF_EXOTIC_AUX sExoticAux;
  CRIF_CALL sCall;
  CRIF_CALL_AUX sCallAux;
  CRIF_MCARG sMCArg;
  LGMSV_MODEL sModel;

  double **dTempVal = NULL;
  double dDfStar;
  double dMarketIV;
  int i;

  sFunding = sDeal->sFunding;
  sFundingAux = sFunding->sAux;
  sExotic = sDeal->sExotic;
  sExoticAux = sExotic->sAux;
  sCall = sDeal->sCall;
  sCallAux = sCall->sAux;

  err = crif_fill_mc_arg(sMarket, sDeal, sSimulArg, sPricingParams);
  if (err)
    goto FREE_RETURN;

  sMCArg = sSimulArg->sMCArg;
  sModel = sSimulArg->sModel;

  /* Fill MCEB Params */
  mceb_set_default_params(sMCArg->sMCEBParams);

  sMCArg->sMCEBParams->iIsKO = 1;
  sMCArg->sMCEBParams->iCallCurrent = 1;
  sMCArg->sMCEBParams->iColPay = 0;
  sMCArg->sMCEBParams->iColBound = 2;
  sMCArg->sMCEBParams->iDoSmoothing = 0;
  sMCArg->sMCEBParams->dCallSpread =
      sPricingParams->sCouponPricingParams->dMCEBCallSpread;
  sMCArg->sMCEBParams->iMultiIndex = 1;
  sMCArg->sMCEBParams->iNbIndex = 3 + sModel->iOne2F;
  sMCArg->iAdjustIV = sPricingParams->sSimulParams->iAdjustIV;
  sMCArg->sMCEBParams->iCalcIV = 1;
  sMCArg->sMCEBParams->iAdjustIV = sPricingParams->sSimulParams->iAdjustIV;
  sMCArg->sMCEBParams->iHasNumeraire = 0;
  sMCArg->sMCEBParams->iColNumeraire = 0;
  sMCArg->sMCEBParams->iCalcOneTime =
      (sPricingParams->sReserveParams->iCalcSwitchAdjust > 0 ? 1 : 0);
  sMCArg->sMCEBParams->iCalcOneTimePartial =
      (sPricingParams->sReserveParams->iCalcSwitchAdjust > 0 ? 1 : 0);
  sMCArg->sMCEBParams->iCalcExeProba = sOutputs->iCalcExeProbas;
  sMCArg->sMCEBParams->iHasFees = sMCArg->iHasCallFee;
  sMCArg->sMCEBParams->iColFees = 5 + sModel->iOne2F;

  err = mceb_allocate_params(sMCArg->sMCEBParams, sMCArg->iNbEvent);
  if (err)
    goto FREE_RETURN;

  /* Save all Model Values */
  sMCArg->dModelValue = calloc(sMCArg->iNbEvent, sizeof(double));

  if (!sMCArg->dModelValue) {
    err = "Memory allocation faillure in crif_calculate_option";
    goto FREE_RETURN;
  }

  dDfStar = swp_f_df(sMarket->lToday, sModel->lTStarDate, sMarket->cYcName);

  /* Populates the Fwd IV's and calculates the Market IV */
  dMarketIV = 0.0;

  for (i = 0; i < sMCArg->sMCEBParams->lNbDates; i++) {
    sMCArg->sMCEBParams->dMarketFwdIV[i] =
        sMCArg->dMarketValue[sMCArg->iIndexEvent[i]] / dDfStar;
    dMarketIV += sMCArg->dMarketValue[sMCArg->iIndexEvent[i]] / dDfStar;
  }

  dTempVal = dmatrix(0, sMCArg->iNbEvent + 4 + 1, 0,
                     1 + 2 + sMCArg->sMCEBParams->iNbIndex +
                         sMCArg->sMCEBParams->iHasNumeraire);

  if (!dTempVal) {
    err = "Memory allocation faillure in crif_calculate_option";
    goto FREE_RETURN;
  }

  smessage("Launching MC  , Nb Steps: %d  , Nb Paths: %d", sMCArg->iNbTimes,
           sPricingParams->sSimulParams->lNbPaths);

  if (sModel->iOne2F == 1) {
    err = lgmSV_mc_balsam_rev(
        sMCArg->iNbTimes, sMCArg->iNbEvent, sMCArg->dTimes, sMCArg->dDates,
        sPricingParams->sSimulParams->lNbPaths, sModel->dLambdaX,
        sMCArg->dSigma, sMCArg->dAlphaSV, sMCArg->dLambdaSV, sMCArg->dLevelSV,
        sMCArg->dRhoSV, NULL, NULL, NULL,
        sPricingParams->sSimulParams->sLGMSVParams, sMCArg->vPayoffParams,
        sMCArg->iEvalEvent, sCall->iIsCallable, sMCArg->iOptimise,
        sMCArg->sMCEBParams, NULL, crif_payoff_1F, NULL, 7, dTempVal);
    if (err)
      goto FREE_RETURN;

    sOutputs->dExotic =
        (dTempVal[1][0] - dTempVal[0][0]) * dDfStar * sCall->dPayRec;
    if (sCall->iIsCallable)
      sOutputs->dCallInit =
          (dTempVal[7][0] - sCall->dPayRec * dTempVal[0][0]) * dDfStar;
    else
      sOutputs->dCallInit = 0;
  } else {
    err = lgmSV2F_mc_balsam_rev(
        sMCArg->iNbTimes, sMCArg->iNbEvent, sMCArg->dTimes, sMCArg->dDates,
        sPricingParams->sSimulParams->lNbPaths, sModel->dLambdaX,
        sModel->dLambdaX2, sMCArg->dSigma, sMCArg->dLGMAlpha, sMCArg->dLGMRho,
        sMCArg->dAlphaSV, sMCArg->dLambdaSV, sMCArg->dLevelSV, sMCArg->dRhoSV,
        sMCArg->dRho2SV, NULL, NULL, NULL, NULL, NULL, NULL,
        sPricingParams->sSimulParams->sLGMSVParams, sMCArg->vPayoffParams,
        sMCArg->iEvalEvent, sCall->iIsCallable, sMCArg->iOptimise,
        sMCArg->sMCEBParams, NULL, crif_payoff_2F, crif_adjust_payoff_2F, 8,
        dTempVal);
    if (err)
      goto FREE_RETURN;

    sOutputs->dExotic =
        (dTempVal[1][0] - dTempVal[0][0]) * dDfStar * sCall->dPayRec;
    if (sCall->iIsCallable)
      sOutputs->dCallInit =
          (dTempVal[8][0] - sCall->dPayRec * dTempVal[0][0]) * dDfStar;
    else
      sOutputs->dCallInit = 0;
  }

  sOutputs->dNewIV = 0.0;
  sOutputs->dExotic += sDeal->sExotic->sAux->dFixedPV;

  for (i = 0; i < sMCArg->sMCEBParams->lNbDates; i++) {
    sMCArg->dModelValue[i] *= dDfStar;
    sOutputs->dNewIV += sMCArg->dModelValue[i];

    if (sMCArg->sMCEBParams->dMarketFwdIV)
      sMCArg->sMCEBParams->dMarketFwdIV[i] *= dDfStar;
    if (sMCArg->sMCEBParams->dModelFwdIV)
      sMCArg->sMCEBParams->dModelFwdIV[i] *= dDfStar;
    if (sMCArg->sMCEBParams->iCalcOneTime)
      sMCArg->sMCEBParams->dOneTimeCall[i] *= dDfStar;
    if (sMCArg->sMCEBParams->iCalcOneTimePartial)
      sMCArg->sMCEBParams->dOneTimePartial[i] *= dDfStar;
  }

FREE_RETURN:

  if (dTempVal)
    free_dmatrix(dTempVal, 0, sMCArg->iNbEvent + 4, 0,
                 2 + sMCArg->sMCEBParams->iNbIndex +
                     sMCArg->sMCEBParams->iHasNumeraire);

  return err;
}

Err crif_calibrate_model(CMSCTS_MARKET sMarket, CRIF_DEAL sDeal,
                         CRIF_SIMULARG sSimulArg, CMSCTS_CALIB sCalibration,
                         CRIF_OUTPUTS sOutputs) {
  Err err = NULL;

  CRIF_EXOTIC sExotic;
  CRIF_CALL sCall;
  CRIF_CALL_AUX sCallAux;
  LGMSV_MODEL sModel;

  char *cPrimFreq, *cPrimBasis, *cPrimRef, *cSecFreq, *cSecBasis, *cSecRef;
  char **cPrimTenors = NULL, **cSecTenors = NULL;

  long lNbPrim, lNbSec;
  long *lPrimExeDates = NULL, *lSecExeDates = NULL;

  long lLongEndDate, lShortEndDate, temp;
  long lPrimEndDate, lSecEndDate;

  int *iPrimCalib = NULL, *iSecCalib = NULL;

  double *dPrimStrike = NULL, *dPrimStrikeS1 = NULL, *dPrimStrikeS2 = NULL,
         *dSecStrike = NULL, *dSecStrikeS1 = NULL, *dSecStrikeS2 = NULL,
         *dSecWeights = NULL;

  int iNbPWTime, iNbPWTimeNew;
  double dLambdaX;
  double *dPWTime = NULL, *dSigma = NULL, *dAlphaSV = NULL, *dLambdaSV = NULL,
         *dRhoSV = NULL, *dRho2SV = NULL;

  double dStepCalib;
  int iStartIndex, iStartCpnIndex;
  int iAddInArrears;

  char cDiagTenor[5] = "DIAG";

  int i, j;
  int iIndexCall, iIndexCpn;

  //	int		iNbSigTime;
  double *dSigTime = NULL;

  /*	double	*tsv[2]  ,
                          *tst[2];
          double	dTime0 = 1.0;

          double	dLambda;*/

  sExotic = sDeal->sExotic;
  sCall = sDeal->sCall;
  sCallAux = sCall->sAux;
  sModel = sSimulArg->sModel;

  /* Finds the number of coupon with gearing >0 */
  iStartCpnIndex = sExotic->sAux->iNumStartIndex;
  while (iStartCpnIndex < sExotic->iNbCpn &&
         fabs(sExotic->dPaidGearing[iStartCpnIndex]) < 1e-8) {
    iStartCpnIndex++;
  }

  if (iStartCpnIndex == sExotic->iNbCpn) // if degenerates into a midat
  {
    iStartCpnIndex = sExotic->sAux->iNumStartIndex;
  }

  switch (sCalibration->eCalibStrat) {
  case CMSCTS_DIAGCMS1: {
    iStartIndex = sCallAux->iStartIndex;
    lNbPrim = sCall->iNbCall - iStartIndex;
    lNbSec = sExotic->iNbCpn - iStartCpnIndex;
    break;
  }

  case CMSCTS_CMS1DIAG: {
    iStartIndex = sCallAux->iStartIndex;
    lNbSec = sCall->iNbCall - iStartIndex;
    lNbPrim = sExotic->iNbCpn - iStartCpnIndex;
    break;
  }
  }

  /* Memory allocation */
  cPrimTenors = calloc(lNbPrim + 1, sizeof(char *));
  lPrimExeDates = calloc(lNbPrim + 1, sizeof(long));
  iPrimCalib = calloc(lNbPrim + 1, sizeof(int));
  dPrimStrike = calloc(lNbPrim + 1, sizeof(double));
  dPrimStrikeS1 = calloc(lNbPrim + 1, sizeof(double));
  dPrimStrikeS2 = calloc(lNbPrim + 1, sizeof(double));

  cSecTenors = calloc(lNbSec + 1, sizeof(char *));
  lSecExeDates = calloc(lNbSec + 1, sizeof(long));
  iSecCalib = calloc(lNbSec + 1, sizeof(int));
  dSecStrike = calloc(lNbSec + 1, sizeof(double));
  dSecStrikeS1 = calloc(lNbSec + 1, sizeof(double));
  dSecStrikeS2 = calloc(lNbSec + 1, sizeof(double));
  dSecWeights = calloc(lNbSec + 1, sizeof(double));

  if (!cPrimTenors || !lPrimExeDates || !iPrimCalib || !dPrimStrike ||
      !dPrimStrikeS1 || !dPrimStrikeS2 || !cSecTenors || !lSecExeDates ||
      !iSecCalib || !dSecStrike || !dSecStrikeS1 || !dSecStrikeS2 ||
      !dSecWeights) {
    err = "Memory allocation faillure in crif_calibrate_model";
    goto FREE_RETURN;
  }

  /* Fill the calibration dates and common informations */
  if (sCallAux->iNbUsedCall > 0) {
    iIndexCall = sCallAux->iStartIndex;
  }

  if (sExotic->sAux->iNbUsedNumCpn > 0) {
    iIndexCpn = iStartCpnIndex;
  }

  iAddInArrears = 0;
  switch (sCalibration->eCalibStrat) {
  case CMSCTS_DIAGCMS1: {
    for (i = 0; i < lNbPrim + 1; i++) {
      if (i < lNbPrim) {
        if (sCallAux->iNbUsedCall > 0 && iIndexCall < sCall->iNbCall) {
          lPrimExeDates[i] = sCall->lExeDate[iIndexCall];
          iIndexCall++;
        } else {
          lPrimExeDates[i] =
              add_unit(sExotic->lStartDate[i + iStartIndex],
                       -sExotic->iFixingLag, SRT_BDAY, MODIFIED_SUCCEEDING);
        }
      } else {
        /* Check if we need an extra date */
        switch (sExotic->eCouponType[sExotic->iNbCpn - 1]) {
        case REGULAR_CRIF: {
          /* Add if In ARREARS fixing */
          if (sExotic->lPaidFixingDate[sExotic->iNbCpn - 1] >
              lPrimExeDates[i - 1] +
                  0.5 * (sExotic->lPayDate[sExotic->iNbCpn - 1] -
                         sExotic->lStartDate[sExotic->iNbCpn - 1])) {
            lPrimExeDates[i] = sExotic->lPaidFixingDate[sExotic->iNbCpn - 1];
            iAddInArrears = 1;
          }

          break;
        }
        }
      }

      iPrimCalib[i] = 1;
      dPrimStrikeS1[i] = sCalibration->dCalibSmileSTD;
      dPrimStrikeS2[i] = -sCalibration->dCalibSmileSTD;
    }

    for (j = 0; j < lNbSec + 1; j++) {
      if (j < lNbSec) {
        if (sExotic->sAux->iNbUsedNumCpn > 0 && iIndexCpn < sExotic->iNbCpn) {
          lSecExeDates[j] = sExotic->lPaidFixingDate[iIndexCpn];
          iIndexCpn++;
        } else {
          lSecExeDates[j] =
              add_unit(sExotic->lStartDate[j + iStartIndex],
                       -sExotic->iFixingLag, SRT_BDAY, MODIFIED_SUCCEEDING);
        }
      } else {
        /* Check if we need an extra date */
        switch (sExotic->eCouponType[sExotic->iNbCpn - 1]) {
        case REGULAR_CRIF: {
          /* Add if In ARREARS fixing */
          if (sExotic->lPaidFixingDate[sExotic->iNbCpn - 1] >
              lPrimExeDates[j - 1] +
                  0.5 * (sExotic->lPayDate[sExotic->iNbCpn - 1] -
                         sExotic->lStartDate[sExotic->iNbCpn - 1])) {
            lSecExeDates[j] = sExotic->lPaidFixingDate[sExotic->iNbCpn - 1];
            iAddInArrears = 1;
          }

          break;
        }
        }
      }

      iSecCalib[j] = 1;
      dSecStrikeS1[j] = sCalibration->dCalibSmileSTD;
      dSecStrikeS2[j] = -sCalibration->dCalibSmileSTD;
      dSecWeights[j] = 1.0;
    }

    break;
  }

  case CMSCTS_CMS1DIAG: {
    for (i = 0; i < lNbPrim + 1; i++) {
      if (i < lNbPrim) {
        if (sExotic->sAux->iNbUsedNumCpn > 0 && iIndexCpn < sExotic->iNbCpn) {
          lPrimExeDates[i] = sExotic->lPaidFixingDate[iIndexCpn];
          iIndexCpn++;
        } else {
          lPrimExeDates[i] =
              add_unit(sExotic->lStartDate[i + iIndexCpn], -sExotic->iFixingLag,
                       SRT_BDAY, MODIFIED_SUCCEEDING);
        }
      } else {
        /* Check if we need an extra date */
        switch (sExotic->eCouponType[sExotic->iNbCpn - 1]) {
        case REGULAR_CRIF: {
          /* Add if In ARREARS fixing */
          if (sExotic->lPaidFixingDate[sExotic->iNbCpn - 1] >
              lPrimExeDates[i - 1] +
                  0.5 * (sExotic->lPayDate[sExotic->iNbCpn - 1] -
                         sExotic->lStartDate[sExotic->iNbCpn - 1])) {
            lPrimExeDates[i] = sExotic->lPaidFixingDate[sExotic->iNbCpn - 1];
            iAddInArrears = 1;
          }

          break;
        }
        }
      }

      iPrimCalib[i] = 1;
      dPrimStrikeS1[i] = sCalibration->dCalibSmileSTD;
      dPrimStrikeS2[i] = -sCalibration->dCalibSmileSTD;
    }

    for (j = 0; j < lNbSec + 1; j++) {
      if (j < lNbSec) {
        if (sCallAux->iNbUsedCall > 0 && iIndexCall < sCall->iNbCall) {
          lSecExeDates[j] = sCall->lExeDate[iIndexCall];
          iIndexCall++;
        } else {
          lSecExeDates[j] =
              add_unit(sExotic->lStartDate[j + iStartIndex],
                       -sExotic->iFixingLag, SRT_BDAY, MODIFIED_SUCCEEDING);
        }
      } else {
        /* Check if we need an extra date */
        switch (sExotic->eCouponType[sExotic->iNbCpn - 1]) {
        case REGULAR_CRIF: {
          /* Add if In ARREARS fixing */
          if (sExotic->lPaidFixingDate[sExotic->iNbCpn - 1] >
              lPrimExeDates[j - 1] +
                  0.5 * (sExotic->lPayDate[sExotic->iNbCpn - 1] -
                         sExotic->lStartDate[sExotic->iNbCpn - 1])) {
            lSecExeDates[j] = sExotic->lPaidFixingDate[sExotic->iNbCpn - 1];
            iAddInArrears = 1;
          }

          break;
        }
        }
      }
      iSecCalib[j] = 1;
      dSecStrikeS1[j] = sCalibration->dCalibSmileSTD;
      dSecStrikeS2[j] = -sCalibration->dCalibSmileSTD;
      dSecWeights[j] = 1.0;
    }

    break;
  }
  }

  /* Check if we are allowed to use in arrears */
  if (!sCalibration->iUseInArrearsFixing)
    iAddInArrears = 0;

  /* Find the long and short end date */
  lLongEndDate = sMarket->lToday;
  lShortEndDate = (long)(sMarket->lToday + 100 * DAYS_IN_YEAR);

  for (i = iStartIndex; i < sExotic->iNbCpn; i++) {
    switch (sExotic->eCouponType[i]) {
    case REGULAR_CRIF: {
      if (fabs(sExotic->dPaidGearing[i]) > 1.0E-08) {
        err = add_tenor(sDeal->lTheoEndDate, sExotic->cPaidTenor[i],
                        NO_BUSDAY_CONVENTION, &temp);
        lLongEndDate = max(lLongEndDate, temp);
        lShortEndDate = min(lShortEndDate, temp);
      }

      break;
    }
    }
  }

  /* fill the strat calibration parameters */
  switch (sCalibration->eCalibStrat) {
  case CMSCTS_DIAGCMS1: {
    /* Primary Instruments */
    cPrimFreq = sMarket->cFixedFreq;
    cPrimBasis = sMarket->cFixedBasis;
    cPrimRef = sMarket->cRefRateName;
    lPrimEndDate = sDeal->lTheoEndDate;
    for (i = 0; i < lNbPrim; i++) {
      cPrimTenors[i] = &cDiagTenor[0];
    }

    /* Secondary Instruments */
    cSecFreq = sExotic->cPaidFreq;
    cSecBasis = sExotic->cPaidBasis;
    cSecRef = sExotic->cPaidRef;
    lSecEndDate = lShortEndDate;
    for (i = 0; i < lNbSec; i++) {
      cSecTenors[i] = sExotic->cPaidTenor[0];
    }

    break;
  }

  case CMSCTS_CMS1DIAG: {
    /* Primary Instruments */
    cPrimFreq = sExotic->cPaidFreq;
    cPrimBasis = sExotic->cPaidBasis;
    cPrimRef = sExotic->cPaidRef;
    lPrimEndDate = lShortEndDate;

    for (i = 0; i < lNbPrim; i++) {
      cPrimTenors[i] = sExotic->cPaidTenor[0];
    }

    /* Secondary Instruments */
    cSecFreq = sMarket->cFixedFreq;
    cSecBasis = sMarket->cFixedBasis;
    cSecRef = sMarket->cRefRateName;
    lSecEndDate = lLongEndDate;

    for (i = 0; i < lNbSec; i++) {
      cSecTenors[i] = &cDiagTenor[0];
    }

    break;
  }
  }

  dLambdaX = sCalibration->dLambda;

  /* Recopy Lambda calibration parameter */
  sCalibration->sLGMSVCalibParams->fix_lambda = sCalibration->iFixLambda;

  err = cpd_calib_diagonal_LGMSV_new_dlm(
      sMarket->cYcName, sMarket->cVcName, sMarket->GetCashVol,
      sMarket->cRefRateName,

      cPrimFreq, cPrimBasis, cPrimRef, lNbPrim, lPrimExeDates, iPrimCalib,
      cPrimTenors, lPrimEndDate, dPrimStrike, dPrimStrikeS1, dPrimStrikeS2,
      sCalibration->sPrimCalibParams,

      cSecFreq, cSecBasis, cSecRef, lNbSec, lSecExeDates, iSecCalib, cSecTenors,
      lSecEndDate, dSecStrike, dSecStrikeS1, dSecStrikeS2, dSecWeights,
      sCalibration->sSecCalibParams,

      sCalibration->sLGMSVCalibParams, sCalibration->iNbFactor, &dLambdaX,
      sCalibration->iNbTimesSV, sCalibration->dTimesSV, sCalibration->dAlphaSV,
      sCalibration->dLamSV, sCalibration->dRhoSV, sCalibration->dTStar,
      sCalibration->dLGMAlpha, sCalibration->dLGMGamma, sCalibration->dLGMRho,
      sCalibration->dRho2SV, sCalibration->sLGMSVNumerParams,

      &iNbPWTime, &dPWTime, &dSigma, &dAlphaSV, &dLambdaSV, &dRhoSV, &dRho2SV,

      sCalibration->sLMParams, sOutputs->sInstDatas);
  if (err)
    goto FREE_RETURN;

  /* Initialise the model with constant vol at the end */
  if (lNbPrim > 1) {
    dStepCalib = (lPrimExeDates[1] - lPrimExeDates[0]) * YEARS_IN_DAY;
  } else {
    dStepCalib = 1.0;
  }

  iNbPWTimeNew = iNbPWTime +
                 (int)(((sDeal->lTheoEndDate - sMarket->lToday) * YEARS_IN_DAY -
                        dPWTime[iNbPWTime - 1]) /
                       dStepCalib) +
                 1;

  dPWTime = (double *)realloc(dPWTime, iNbPWTimeNew * sizeof(double));
  dSigma = (double *)realloc(dSigma, iNbPWTimeNew * sizeof(double));
  dAlphaSV = (double *)realloc(dAlphaSV, iNbPWTimeNew * sizeof(double));
  dLambdaSV = (double *)realloc(dLambdaSV, iNbPWTimeNew * sizeof(double));
  dRhoSV = (double *)realloc(dRhoSV, iNbPWTimeNew * sizeof(double));
  if (dRho2SV)
    dRho2SV = (double *)realloc(dRho2SV, iNbPWTimeNew * sizeof(double));

  if (!dPWTime || !dSigma || !dAlphaSV || !dLambdaSV || !dRhoSV ||
      (!dRho2SV && sCalibration->iNbFactor == 2)) {
    err = "Memory allocation faillure (2) in cmscts_calibrate_model";
    goto FREE_RETURN;
  }

  for (i = 0; i < iNbPWTime; i++) {
    dAlphaSV[i] *= 2.0;
    dLambdaSV[i] *= 2.0;
  }

  for (i = iNbPWTime; i < iNbPWTimeNew; i++) {
    dPWTime[i] = dPWTime[i - 1] + dStepCalib;
    dSigma[i] = dSigma[iNbPWTime - 1];
    dAlphaSV[i] = dAlphaSV[iNbPWTime - 1];
    dLambdaSV[i] = dLambdaSV[iNbPWTime - 1];
    dRhoSV[i] = dRhoSV[iNbPWTime - 1];
    if (dRho2SV)
      dRho2SV[i] = dRho2SV[iNbPWTime - 1];
  }

  err = init_LGMSV_model(sModel, sMarket->lToday, sCalibration->iNbFactor,
                         iNbPWTimeNew, dLambdaX, dPWTime, dSigma, dAlphaSV,
                         dLambdaSV, dLambdaSV, dRhoSV, sCalibration->dTStar,
                         sCalibration->dLGMAlpha, sCalibration->dLGMGamma,
                         sCalibration->dLGMRho, dRho2SV);

  if (err)
    goto FREE_RETURN;

  ConvertTS_LGM_to_LGMSV(
      sModel->iNbPWTime, sModel->dPWTime, sModel->dSigma, sModel->dLambdaX,
      sModel->dTStar, sModel->iOne2F, sModel->dInitLGMAlpha, sModel->dLGMGamma,
      sModel->dInitLGMRho, sModel->dLGMAlpha, sModel->dLGMRho);
/*
        // calib du 2 Facteurs
        err = cpd_calib_diagonal_dlm(	sMarket->cYcName  ,
                                                                        sMarket->cVcName
   , sMarket->GetCashVol  , sMarket->cRefRateName  ,

                                                                        cPrimFreq
   , cPrimBasis  , cPrimRef  , lNbPrim  , lPrimExeDates  , iPrimCalib  ,
                                                                        cPrimTenors
   , sDeal->lTheoEndDate  , dPrimStrike  , sCalibration->sPrimCalibParams  ,

                                                                        cSecFreq
   , cSecBasis  , cSecRef  , lNbSec  , lSecExeDates  , iSecCalib  , cSecTenors ,
                                                                        sDeal->lTheoEndDate
   , dSecStrike  , dSecWeights  , sCalibration->sSecCalibParams  ,

                                                                        sCalibration->sLGMSVCalibParams
   ->fix_lambda  , sCalibration->iNbFactor  , 0  , &dTime0  , &dLambda  , 0  ,
                                                                        2  ,
                                                                        sCalibration->dLGMAlpha
   , sCalibration->dLGMGamma  , sCalibration->dLGMRho  ,


                                                                        //
   Shift Parameters 0  , NULL  , NULL  ,

                                                                        0  ,
                                                                        NULL  ,
                                                                        NULL  ,

                                                                        &iNbSigTime
   , &dSigTime  , &dSigma  , sCalibration->sLMParams  ,
                                                                        sOutputs->sLGM2F_FixedTau_InstDatas);
        if (err) goto FREE_RETURN;

        // Initialise the Und

        tsv[0] = dSigTime;
        tsv[1] = dSigma;

        for (i=0; i<iNbSigTime; i++)
        {
                dSigTime[i] = sMarket->lToday + DAYS_IN_YEAR * dSigTime[i];
        }

        tst[0] = &dTime0;
        tst[1] = &dLambdaX;

        dLambda = 1.0 / dLambda;
        dTime0 = sMarket->lToday + DAYS_IN_YEAR * dTime0;

        err = SrtInitIRUnd(	sSimulArg->cLGM_fixedTau_ModelName  ,
                                                sMarket->cYcName  ,
                                                (sCalibration->iNbFactor == 1?
   "LGM": "LGM2F")  , iNbSigTime  , 2  , tsv  , 1  , 2  , tst  , 0.0  ,
                                                sCalibration->dLGMAlpha  ,
                                                sCalibration->dLGMGamma  ,
                                                sCalibration->dLGMRho  ,
                                                0  ,  0  , 0  , 0  , 0  , 0);

        if (err) goto FREE_RETURN;
        */
FREE_RETURN:

  if (cPrimTenors)
    free(cPrimTenors);
  if (lPrimExeDates)
    free(lPrimExeDates);
  if (iPrimCalib)
    free(iPrimCalib);
  if (dPrimStrike)
    free(dPrimStrike);
  if (dPrimStrikeS1)
    free(dPrimStrikeS1);
  if (dPrimStrikeS2)
    free(dPrimStrikeS2);

  if (cSecTenors)
    free(cSecTenors);
  if (lSecExeDates)
    free(lSecExeDates);
  if (iSecCalib)
    free(iSecCalib);
  if (dSecStrike)
    free(dSecStrike);
  if (dSecStrikeS1)
    free(dSecStrikeS1);
  if (dSecStrikeS2)
    free(dSecStrikeS2);
  if (dSecWeights)
    free(dSecWeights);

  if (dPWTime)
    free(dPWTime);
  if (dSigma)
    free(dSigma);
  if (dAlphaSV)
    free(dAlphaSV);
  if (dLambdaSV)
    free(dLambdaSV);
  if (dRhoSV)
    free(dRhoSV);
  if (dRho2SV)
    free(dRho2SV);

  return err;
}

void crif_free_pay_arg(CRIF_PAYARG sPayArg) {
  if (sPayArg->sCallArg)
    free(sPayArg->sCallArg);
  if (sPayArg->sCouponArg)
    free(sPayArg->sCouponArg);
}

void crif_free_mc_arg(CRIF_MCARG sMCArg) {
  int i;
  CRIF_PAYARG sPayArg;

  if (sMCArg->dTimes)
    free(sMCArg->dTimes);
  if (sMCArg->dDates)
    free(sMCArg->dDates);
  if (sMCArg->dSigma)
    free(sMCArg->dSigma);
  if (sMCArg->dAlphaSV)
    free(sMCArg->dAlphaSV);
  if (sMCArg->dLambdaSV)
    free(sMCArg->dLambdaSV);
  if (sMCArg->dLevelSV)
    free(sMCArg->dLevelSV);
  if (sMCArg->dRhoSV)
    free(sMCArg->dRhoSV);
  if (sMCArg->dRho2SV)
    free(sMCArg->dRho2SV);
  if (sMCArg->dLGMAlpha)
    free(sMCArg->dLGMAlpha);
  if (sMCArg->dLGMRho)
    free(sMCArg->dLGMRho);

  if (sMCArg->iEvalEvent)
    free(sMCArg->iEvalEvent);
  if (sMCArg->iIndexEvent)
    free(sMCArg->iIndexEvent);

  if (sMCArg->vPayoffParams) {
    for (i = 0; i < sMCArg->iNbTimes; i++) {
      sPayArg = (CRIF_PAYARG)(sMCArg->vPayoffParams[i]);

      if (sPayArg) {
        crif_free_pay_arg(sPayArg);
        free(sPayArg);
      }
    }

    free(sMCArg->vPayoffParams);
  }

  if (sMCArg->dPathInfos)
    free(sMCArg->dPathInfos);
  if (sMCArg->dSavedCoupon)
    free(sMCArg->dSavedCoupon);
  if (sMCArg->dSavedDF)
    free(sMCArg->dSavedDF);

  if (sMCArg->iOptimise)
    free(sMCArg->iOptimise);
  if (sMCArg->dMarketValue)
    free(sMCArg->dMarketValue);

  if (sMCArg->sMCEBParams) {
    mceb_free_params(sMCArg->sMCEBParams);
    free(sMCArg->sMCEBParams);
  }
}
