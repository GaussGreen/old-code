
#include "CRIFProdStruct.h"
#include "CMSCTSPricing.h"
#include "srt_h_allFx3F.h"

#define BIG_BIG_NUMBER 1E40
#define ONE_MONTH 0.083333333

Err crif_allocate_deal(CRIF_DEAL sDeal) {
  sDeal->sFunding = calloc(1, sizeof(cmscts_funding));
  sDeal->sExotic = calloc(1, sizeof(crif_exotic));
  sDeal->sCall = calloc(1, sizeof(crif_call));

  if (!sDeal->sFunding || !sDeal->sExotic || !sDeal->sCall) {
    return "Memory allocation faillure in crif_allocate_deal";
  }

  return NULL;
}

void crif_free_deal(CRIF_DEAL sDeal) {
  if (sDeal->sFunding) {
    cmscts_free_funding(sDeal->sFunding);
    free(sDeal->sFunding);
  }

  if (sDeal->sExotic) {
    crif_free_exotic(sDeal->sExotic);
    free(sDeal->sExotic);
  }

  if (sDeal->sCall) {
    crif_free_call(sDeal->sCall);
    free(sDeal->sCall);
  }
}

Err crif_fill_funding(CRIF_DEAL sDeal, CMSCTS_MARKET sMarket) {
  Err err = NULL;
  CMSCTS_FUNDING sFunding;
  CRIF_EXOTIC sExotic;
  CMSCTS_FUND_AUX sAux;

  int iStartIndex, i;
  double dDomNotional, dFundNotional;

  /* Initialisation */
  sFunding = sDeal->sFunding;
  sExotic = sDeal->sExotic;

  /* Allocation */
  err = cmscts_funding_allocate_aux(sFunding);
  if (err)
    goto FREE_RETURN;

  sAux = sFunding->sAux;

  /* Recopy funding infos into Aux */
  memcpy(sAux->dNot, sFunding->dNot, sFunding->iNbCpn * sizeof(double));
  memcpy(sAux->dMargin, sFunding->dMargin, sFunding->iNbCpn * sizeof(double));
  memcpy(sAux->dFixCoupon, sFunding->dFixCoupon,
         sFunding->iNbCpn * sizeof(double));
  strcpy(sAux->cFundYcName, sMarket->cFundYcName);

  /*	Skip coupons fixed before today */
  iStartIndex = 0;
  while (iStartIndex < sFunding->iNbCpn &&
         sFunding->lPayDate[iStartIndex] <
             sMarket->lToday + sMarket->iEODPayFlag) {
    iStartIndex++;
  }

  sAux->iStartIndex = iStartIndex;

  /* Fill the Spreads */
  for (i = iStartIndex; i < sFunding->iNbCpn; i++) {
    if (sFunding->lStartDate[i] >= sMarket->lToday) {
      sAux->dInitSpread[i] = swp_f_spread(
          sFunding->lStartDate[i], sFunding->lPayDate[i], sFunding->cRefRate);
      sAux->dNewSpread[i] = swp_f_spread(
          sFunding->lStartDate[i], sFunding->lPayDate[i], sFunding->cRefRate);
    }
  }

  /* Conversion of foreign funding */
  if (sFunding->iCCY) {
    /* Check that we have a fixed notional */
    dDomNotional = sExotic->dNot[0];

    for (i = 1; i < sExotic->iNbCpn; i++) {
      if (fabs(dDomNotional - sExotic->dNot[i]) > 1.0E-08) {
        err = "No Term Structure of Notional allowed with foreign funding";
        goto FREE_RETURN;
      }
    }

    dFundNotional = sFunding->dNot[0];

    for (i = 1; i < sFunding->iNbCpn; i++) {
      if (fabs(dFundNotional - sFunding->dNot[i]) > 1.0E-08) {
        err = "No Term Structure of Notional allowed with foreign funding";
        goto FREE_RETURN;
      }
    }

    err = convert_funding_to_domestic(
        sMarket->lToday, sDeal->lStartDate, sMarket->iEODFixFlag,
        sMarket->iEODPayFlag, sMarket->dFundSpotFx, sMarket->lFundFxSpotDate,
        dDomNotional, sMarket->cYcName, sFunding->iNbCpn, sFunding->lFixDate,
        sFunding->lStartDate, sFunding->lPayDate, NULL, /* not used */
        sAux->cFundYcName, &dFundNotional, sAux->dNewSpread, sAux->dMargin,
        sAux->dFixCoupon, &sAux->lStartDate, &sAux->dFinalEquiEx,
        &sAux->dInitEquiEx);

    if (err)
      goto FREE_RETURN;

    /* Fill the new notional TS */
    for (i = 0; i < sFunding->iNbCpn; i++) {
      sAux->dNot[i] = dFundNotional;
    }
  }

  /* Calculate the Cash Fx */
  sAux->dCashFx =
      sMarket->dFundSpotFx *
      swp_f_df(sMarket->lToday, sMarket->lFundFxSpotDate, sMarket->cYcName) /
      swp_f_df(sMarket->lToday, sMarket->lFundFxSpotDate, sMarket->cFundYcName);

  /* Fill the Coupon Exchanges */
  for (i = iStartIndex; i < sFunding->iNbCpn; i++) {
    sAux->dCouponPlusEx[i] =
        (sFunding->dCoverage[i] * (sAux->dMargin[i] + sAux->dNewSpread[i]) -
         1.0) *
        sAux->dNot[i];
    sAux->dCouponPlusExPartial[i] = sAux->dCouponPlusEx[i];

    if (i < sFunding->iNbCpn - 1) {
      sAux->dCouponPlusEx[i] += sAux->dNot[i + 1];
    }
  }

FREE_RETURN:

  if (err) {
    cmscst_funding_free_aux(sFunding);
  }

  return err;
}

Err crif_interp_coupon_type(char *cCouponType, CRIF_COUPON_TYPE *eCouponType) {
  char str[30 + 1];

  if (!cCouponType)
    return "Empty string in crif_interp_coupon_type";

  strncpy(str, cCouponType, 30);
  str[30] = '\0';

  strupper(str);
  strip_white_space(str);

  if (!strcmp(str, "REGULAR")) {
    *eCouponType = REGULAR_CRIF;
    return 0;
  }
  return serror("unknown crif coupon type. %s", str);
}

Err crif_exotic_allocate_aux(CRIF_EXOTIC sExotic) {
  CRIF_EXOTIC_AUX sAux;

  sExotic->sAux = calloc(1, sizeof(crif_exotic_aux));
  sAux = sExotic->sAux;
  if (!sAux)
    return "Memory allocation faillure in crif_exotic_allocate_aux";

  sAux->iFloatingType = calloc(sExotic->iNbCpn, sizeof(int));
  sAux->iHasFloorCap = calloc(sExotic->iNbCpn, sizeof(int));
  sAux->dMarketValue = calloc(sExotic->iNbCpn, sizeof(double));
  sAux->dCVStrike = calloc(sExotic->iNbCpn, sizeof(double));

  if (!sAux->dCVStrike || !sAux->iFloatingType || !sAux->iHasFloorCap ||
      !sAux->dMarketValue) {
    return "Memory allocation error in crif_exotic_allocate_aux";
  }

  return NULL;
}

void crif_exotic_free_aux(CRIF_EXOTIC sExotic) {
  CRIF_EXOTIC_AUX sAux;

  sAux = (sExotic->sAux);

  if (sAux->iFloatingType)
    free(sAux->iFloatingType);
  if (sAux->iHasFloorCap)
    free(sAux->iHasFloorCap);
  if (sAux->dMarketValue)
    free(sAux->dMarketValue);
}

void crif_free_exotic(CRIF_EXOTIC sExotic) {

  if (sExotic->sAux) {
    crif_exotic_free_aux(sExotic);
    free(sExotic->sAux);
  }

  if (sExotic->lStartDate)
    free(sExotic->lStartDate);
  if (sExotic->lEndDate)
    free(sExotic->lEndDate);
  if (sExotic->lPayDate)
    free(sExotic->lPayDate);
  if (sExotic->dCoverage)
    free(sExotic->dCoverage);
  if (sExotic->dNot)
    free(sExotic->dNot);
  if (sExotic->eCouponType)
    free(sExotic->eCouponType);

  if (sExotic->lPaidFixingDate)
    free(sExotic->lPaidFixingDate);
  if (sExotic->dPaidMargin)
    free(sExotic->dPaidMargin);

  if (sExotic->dFloor)
    free(sExotic->dFloor);
  if (sExotic->dCap)
    free(sExotic->dCap);

  if (sExotic->dPrvCpngearing)
    free(sExotic->dPrvCpngearing);
  if (sExotic->dPaidGearing)
    free(sExotic->dPaidGearing);
  if (sExotic->cPaidTenor)
    free_svector_size(sExotic->cPaidTenor, 0, 0, 5);
  /*	if (sExotic->cPaidBasis) free_svector_size(sExotic->cPaidBasis  , 0  ,
     1-1  , 5); if (sExotic->cPaidFreq) free_svector_size(sExotic->cPaidFreq  ,
     0  , 1-1  , 5);
          if (sExotic->cPaidRef) free_svector_size(sExotic->cPaidRef  , 0  , 1-1
     , 10);	*/

  if (sExotic->dPastPaidCpn)
    free(sExotic->dPastPaidCpn);
}

Err crif_check_exotic(CRIF_EXOTIC sExotic) {
  int i;

  for (i = 0; i < sExotic->iNbCpn; i++) {
    if (i > 0) {
      /*	Check that start and pay dates are increasing */
      if (sExotic->lStartDate[i] < sExotic->lStartDate[i - 1]) {
        return "Start dates should be increasing in Exotic leg";
      }

      if (sExotic->lPayDate[i] < sExotic->lPayDate[i - 1]) {
        return "Pay dates should be increasing in Exotic leg";
      }
    }

    if (sExotic->lPayDate[i] < sExotic->lStartDate[i]) {
      return "Pay dates should be after start dates in Exotic leg";
    }
  }

  /*	OK */
  return NULL;
}

Err crif_fill_exotic(CRIF_DEAL sDeal, CMSCTS_MARKET sMarket,
                     CRIF_CPN_NUMPARAMS sNumerParams) {
  Err err = NULL;
  CRIF_EXOTIC sExotic;
  CRIF_EXOTIC_AUX sAux;

  int iStartIndex, i;

  /* Initialisation */
  sExotic = sDeal->sExotic;

  /* Memory allocation */
  err = crif_exotic_allocate_aux(sExotic);
  if (err)
    goto FREE_RETURN;

  sAux = sExotic->sAux;

  /*	Skip coupons fixed before today */
  iStartIndex = 0;
  while (iStartIndex < sExotic->iNbCpn &&
         sExotic->lPayDate[iStartIndex] <
             sMarket->lToday + sMarket->iEODPayFlag) {
    iStartIndex++;
  }

  sAux->iStartIndex = iStartIndex;

  for (i = iStartIndex; i < sExotic->iNbCpn; i++) {
    /* Fill the payoff */
    switch (sExotic->eCouponType[i]) {
    case REGULAR_CRIF: {
      sAux->iFloatingType[i] = 0;

      if (fabs(sExotic->dPaidGearing[i]) > 1.0E-08) {
        sAux->iFloatingType[i] = 1;
        break;
      }

      break;
    }

    default: {
      err = "Coupon type not recognised in crif_fill_exotic";
      goto FREE_RETURN;
    }
    }
  }

FREE_RETURN:

  return err;
}

Err crif_call_allocate_aux(CRIF_CALL sCall) {
  CRIF_CALL_AUX sAux;

  sCall->sAux = calloc(1, sizeof(crif_call_aux));
  sAux = sCall->sAux;
  if (!sAux)
    return "Memory allocation faillure in crif_call";

  sAux->iStartCpnIdx = calloc(sCall->iNbCall, sizeof(int));
  sAux->iStartFundIdx = calloc(sCall->iNbCall, sizeof(int));

  if (!sAux->iStartCpnIdx || !sAux->iStartFundIdx) {
    return "Memory allocation faillure in crif_call_allocate_aux";
  }

  return NULL;
}

void crif_call_free_aux(CRIF_CALL sCall) {
  CRIF_CALL_AUX sAux;

  sAux = sCall->sAux;

  if (sAux->iStartFundIdx)
    free(sAux->iStartFundIdx);
  if (sAux->iStartCpnIdx)
    free(sAux->iStartCpnIdx);
}

void crif_free_call(CRIF_CALL sCall) {
  if (sCall->sAux) {
    crif_call_free_aux(sCall);
    free(sCall->sAux);
  }

  if (sCall->lExeDate)
    free(sCall->lExeDate);
  if (sCall->lSettlDate)
    free(sCall->lSettlDate);
  if (sCall->dFee)
    free(sCall->dFee);
}

Err crif_check_call(CRIF_CALL sCall, CRIF_EXOTIC sExotic) {
  int i;
  CRIF_CALL_AUX sAux;

  sAux = sCall->sAux;

  for (i = sAux->iStartIndex; i < sCall->iNbCall; i++) {
    if (i > sAux->iStartIndex) {
      if (sCall->lExeDate[i] <= sCall->lExeDate[i - 1]) {
        return "Exercise dates should be increasing";
      }

      if (sCall->lSettlDate[i] <= sCall->lSettlDate[i - 1]) {
        return "Settlement dates should be increasing";
      }

      if (sAux->iStartFundIdx[i] < sAux->iStartFundIdx[i - 1]) {
        return "Number of funding coupons controlled by calls should be "
               "decreasing";
      }

      if (sAux->iStartCpnIdx[i] < sAux->iStartCpnIdx[i]) {
        return "Number of exotic coupons controlled by calls should be "
               "decreasing";
      }

      if (sAux->iStartFundIdx[i] <= sAux->iStartFundIdx[i - 1] &&
          sAux->iStartCpnIdx[i] <= sAux->iStartCpnIdx[i - 1]) {
        return serror("Calls %d and %d -indexed after today- are redundant", i,
                      i + 1);
      }
    }

    if (i > sCall->sAux->iStartIndex) {
      if (sExotic->sAux->iFloatingType[sCall->sAux->iStartCpnIdx[i - 1]] > 0 &&
          sCall->lExeDate[i] <=
              sExotic->lPaidFixingDate[sCall->sAux->iStartCpnIdx[i - 1]]) {
        return "Exercise date cannot be smaller or equal to the floating "
               "fixing date";
      }
    }
  }

  /*	OK */
  return NULL;
}

Err crif_fill_call(CRIF_DEAL sDeal, CMSCTS_MARKET sMarket) {
  CMSCTS_FUNDING sFunding;
  CRIF_EXOTIC sExotic;
  CRIF_CALL sCall;
  CRIF_CALL_AUX sAux;

  int iStartIndex, i, iFundIdx, iCpnIdx;
  Err err = NULL;

  /* Initialisation */
  sFunding = sDeal->sFunding;
  sExotic = sDeal->sExotic;
  sCall = sDeal->sCall;

  /* Memory allocation */
  err = crif_call_allocate_aux(sCall);
  if (err)
    goto FREE_RETURN;

  sAux = sCall->sAux;

  /*	Skip coupons fixed before today */
  iStartIndex = 0;
  while (iStartIndex < sCall->iNbCall &&
         sCall->lExeDate[iStartIndex] <
             sMarket->lToday + sMarket->iEODExeFlag) {
    iStartIndex++;
  }

  sAux->iStartIndex = iStartIndex;
  sAux->iNbUsedCall = sCall->iNbCall - sAux->iStartIndex;

  iFundIdx = 0;
  iCpnIdx = 0;

  for (i = 0; i < sCall->iNbCall; i++) {
    /* Find the first funding coupon to be called */
    while (iFundIdx < sFunding->iNbCpn &&
           sFunding->lStartDate[iFundIdx] < sCall->lExeDate[i]) {
      iFundIdx++;
    }
    if (iFundIdx == sFunding->iNbCpn) {
      err = serror("Call number %d does not control any coupon in funding leg",
                   i + 1);
      goto FREE_RETURN;
    }
    sAux->iStartFundIdx[i] = iFundIdx;

    sFunding->sAux->iFundingControll[iFundIdx] = 1;

    /* Find the first exotic coupon to be called */
    while (iCpnIdx < sExotic->iNbCpn &&
           sExotic->lStartDate[iCpnIdx] < sCall->lExeDate[i]) {
      iCpnIdx++;
    }
    if (iCpnIdx == sExotic->iNbCpn) {
      err = serror("Call number %d does not control any coupon in exotic leg",
                   i + 1);
      goto FREE_RETURN;
    }
    sAux->iStartCpnIdx[i] = iCpnIdx;
  }

  /* Exercised */
  if (sCall->iIsExercised) {
    i = 0;

    while (i < sCall->iNbCall &&
           fabs(sCall->lExeDate[i] - sCall->lExercisedDate) > 0.5) {
      i++;
    }

    if (i == sCall->iNbCall) {
      err = "Exercised Date is not in the Exercise Schedule";
      goto FREE_RETURN;
    }

    sAux->iIndexExercised = i;
  }

FREE_RETURN:

  return err;
}

Err crif_fill_check_all_struct(CRIF_DEAL sDeal, CMSCTS_MARKET sMarket,
                               CRIF_CPN_NUMPARAMS sNumerParams) {
  Err err = NULL;

  CMSCTS_FUNDING sFunding;
  CRIF_EXOTIC sExotic;

  /* Funding */
  sFunding = sDeal->sFunding;

  err = crif_fill_funding(sDeal, sMarket);
  if (err)
    goto FREE_RETURN;

  err = cmscts_check_funding(sDeal->sFunding);
  if (err)
    goto FREE_RETURN;

  /* Exotic */
  sExotic = sDeal->sExotic;

  err = crif_fill_exotic(sDeal, sMarket, sNumerParams);
  if (err)
    goto FREE_RETURN;

  err = crif_check_exotic(sExotic);
  if (err)
    goto FREE_RETURN;

  /* Theo End Date */
  sDeal->lTheoEndDate = sExotic->lPayDate[sExotic->iNbCpn - 1];

  /* Call */
  err = crif_fill_call(sDeal, sMarket);
  if (err)
    goto FREE_RETURN;

  err = crif_check_call(sDeal->sCall, sExotic);
  if (err)
    goto FREE_RETURN;

  /* Find the End Index for Funding */
  sFunding->sAux->iEndIndex = sFunding->iNbCpn - 1;

  if (sDeal->sCall->iNbCall > 0 && sDeal->sCall->iIsExercised) {
    sFunding->sAux->iEndIndex =
        sDeal->sCall->sAux->iStartFundIdx[sDeal->sCall->sAux->iIndexExercised] -
        1;
  }

  sFunding->sAux->iNbUsedCpn =
      sFunding->sAux->iEndIndex - sFunding->sAux->iStartIndex + 1;

  /* Find the End Index for Exotic */
  sExotic->sAux->iEndIndex = sExotic->iNbCpn - 1;

  if (sDeal->sCall->iNbCall > 0 && sDeal->sCall->iIsExercised) {
    sExotic->sAux->iEndIndex =
        sDeal->sCall->sAux->iStartCpnIdx[sDeal->sCall->sAux->iIndexExercised] -
        1;
  }

  sExotic->sAux->iNbUsedCpn =
      sExotic->sAux->iEndIndex - sExotic->sAux->iStartIndex + 1;

  /* find the start index exotic for numeric */
  if (sDeal->sCall->iIsExercised == 0) {
    // no closed form for the underlying so need to evaluate numerically from
    // the the coupon fixing aftertoday
    sExotic->sAux->iNumStartIndex = sExotic->sAux->iStartIndex;
    while (sExotic->sAux->iNumStartIndex < sExotic->iNbCpn &&
           sExotic->lPaidFixingDate[sExotic->sAux->iNumStartIndex] <
               sMarket->lToday + sMarket->iEODPayFlag) {
      sExotic->sAux->iNumStartIndex++;
    };
  } else {
    sExotic->sAux->iNumStartIndex = sExotic->iNbCpn;
  }

  sExotic->sAux->iNbUsedNumCpn =
      sExotic->sAux->iEndIndex - sExotic->sAux->iNumStartIndex + 1;

  /* find the start index funding for numeric */
  if (sExotic->sAux->iNumStartIndex < sExotic->iNbCpn) {
    sDeal->sFunding->sAux->iNumStartIndex = sDeal->sFunding->sAux->iStartIndex;

    while (sDeal->sFunding->sAux->iNumStartIndex < sDeal->sFunding->iNbCpn &&
           sDeal->sFunding->lStartDate[sDeal->sFunding->sAux->iNumStartIndex] <
               sExotic->lStartDate[sExotic->sAux->iNumStartIndex] - 10) {
      sDeal->sFunding->sAux->iNumStartIndex++;
    }

    if (sDeal->sFunding->sAux->iNumStartIndex == sDeal->sFunding->iNbCpn) {
      err = serror("Exotic coupon number %d does not have a funding equivalent",
                   sExotic->sAux->iNumStartIndex);
      goto FREE_RETURN;
    }
  } else {
    sDeal->sFunding->sAux->iNumStartIndex = sDeal->sFunding->iNbCpn;
  }

FREE_RETURN:

  return err;
}

Err crif_interp_calib_strat(const char *constStr, CMSCTS_CALIB_STRAT *val) {
  char str[30 + 1];

  if (!constStr)
    return "Empty string in crif_interp_calib_strat";

  strncpy(str, constStr, 30);
  str[30] = '\0';

  strupper(str);
  strip_white_space(str);

  if (!strcmp(str, "DIAGCMS")) {
    *val = CMSCTS_DIAGCMS1;
    return 0;
  }
  if (!strcmp(str, "CMSDIAG")) {
    *val = CMSCTS_CMS1DIAG;
    return 0;
  }

  return serror("unknown crif_calib_strat. %s", str);
}

Err crif_allocate_pricing_params(CRIF_PRICING_PARAMS sPricingParams) {
  sPricingParams->sCouponPricingParams = calloc(1, sizeof(crif_cpn_numparams));
  sPricingParams->sSimulParams = calloc(1, sizeof(cmscts_simul_params));
  sPricingParams->sReserveParams = calloc(1, sizeof(cmscts_reserve_params));

  if (!sPricingParams->sCouponPricingParams || !sPricingParams->sSimulParams ||
      !sPricingParams->sReserveParams) {
    return "Memory allocation faillure in crif_allocate_pricing_params";
  }

  sPricingParams->sSimulParams->sLGMSVParams = calloc(1, sizeof(LGMSVParam));

  return NULL;
}

void crif_free_pricing_params(CRIF_PRICING_PARAMS sPricingParams) {
  if (sPricingParams->sCouponPricingParams) {
    free(sPricingParams->sCouponPricingParams);
  }

  if (sPricingParams->sSimulParams) {
    if (sPricingParams->sSimulParams->sLGMSVParams) {
      free(sPricingParams->sSimulParams->sLGMSVParams);
    }

    free(sPricingParams->sSimulParams);
  }

  if (sPricingParams->sReserveParams) {
    free(sPricingParams->sReserveParams);
  }
}
