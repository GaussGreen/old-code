// ------------------------------------------------------------------------------------------------------------------
// //
//
// CRIF.cxx
//
//

#include "CRIF.h"
#include "CMSCTSPricing.h"
#include "CRIFProdStruct.h"

Err crif_calc_fund_value(CRIF_DEAL sDeal, CMSCTS_MARKET sMarket) {
  CMSCTS_FUNDING sFunding;
  CRIF_EXOTIC sExotic;
  CMSCTS_FUND_AUX sAux;
  int i;
  double df_start, df_end;
  Err err = NULL;

  sFunding = sDeal->sFunding;
  sExotic = sDeal->sExotic;
  sAux = sFunding->sAux;

  for (i = sAux->iStartIndex; i <= sAux->iEndIndex; i++) {
    df_end =
        swp_f_df(sMarket->lToday, sFunding->lPayDate[i], sMarket->cFundYcName);

    if (sFunding->lFixDate[i] < sMarket->lToday + sMarket->iEODFixFlag) {
      sAux->dMarketValue[i] = (sFunding->dFixCoupon[i] + sFunding->dMargin[i]) *
                              sFunding->dCoverage[i] * df_end *
                              sFunding->dNot[i] * sAux->dCashFx;
    } else {
      /* Libor calculation */
      df_start = swp_f_df(sMarket->lToday, sFunding->lStartDate[i],
                          sMarket->cFundYcName);

      sAux->dMarketValue[i] = (df_start - df_end +
                               (sFunding->dMargin[i] + sAux->dInitSpread[i]) *
                                   sFunding->dCoverage[i] * df_end) *
                              sFunding->dNot[i] * sAux->dCashFx;
    }

    /* for foreign funding add initial and final notional in Funding and Exotic
     * leg */
    if (sFunding->iCCY) {
      if (sFunding->lStartDate[i] > sMarket->lToday + sMarket->iEODPayFlag) {
        sAux->dMarketValue[i] -=
            sFunding->dNot[0] *
            swp_f_df(sMarket->lToday, sFunding->lStartDate[i],
                     sMarket->cFundYcName) *
            sFunding->sAux->dCashFx;
        sAux->dMarketValue[i] +=
            sExotic->dNot[0] * swp_f_df(sMarket->lToday,
                                        sFunding->lStartDate[i],
                                        sMarket->cYcName);
      }

      if (sFunding->lPayDate[i] > sMarket->lToday + sMarket->iEODPayFlag) {
        sAux->dMarketValue[i] +=
            sFunding->dNot[0] *
            swp_f_df(sMarket->lToday, sFunding->lPayDate[i],
                     sMarket->cFundYcName) *
            sFunding->sAux->dCashFx;
        sAux->dMarketValue[i] -=
            sExotic->dNot[0] *
            swp_f_df(sMarket->lToday, sFunding->lPayDate[i], sMarket->cYcName);
      }
    }
  }

  return err;
}

// calc  historical values for ratchet leg
Err crif_calc_exotic_value(CRIF_DEAL sDeal, CMSCTS_MARKET sMarket,
                           CRIF_PRICING_PARAMS sPricingParams,
                           CMSCTS_CALIB sCalibration) {
  CRIF_EXOTIC sExotic;
  CRIF_EXOTIC_AUX sAux;
  CRIF_CPN_NUMPARAMS sNumParams;

  SwapDP sFixedSwap;
  SrtDateList sFixedSwapDates;

  int i, j, iStartIndex;
  double dCouponValue, dPreviousCouponValue;

  double *dValues = NULL;
  double dPaidDf, dLevel, dFloating, dCashSwap;

  char **cPairNames = NULL;
  double **dPairCorrels = NULL, *dOutput = NULL;

  long lTheoEndDate;

  Err err = NULL;

  sExotic = sDeal->sExotic;
  sAux = sExotic->sAux;
  sNumParams = sPricingParams->sCouponPricingParams;

  /* Start Index */
  iStartIndex = sAux->iStartIndex;

  /* Must Value Past Coupons for the Ratchet*/
  dPreviousCouponValue = 0.0;
  for (i = iStartIndex; i <= sExotic->sAux->iEndIndex; i++) {
    dCouponValue = 0.0;
    if (sExotic->lPayDate[i] > sMarket->lToday) {
      dPaidDf =
          swp_f_df(sMarket->lToday, sExotic->lPayDate[i], sMarket->cYcName);
    } else {
      dPaidDf = 1.0;
    }

    /*  do the pricing */
    switch (sExotic->eCouponType[i]) {
    case REGULAR_CRIF: {

      if (sExotic->lPaidFixingDate[i] <
              sMarket->lToday + sMarket->iEODFixFlag ||
          sAux->iFloatingType[i] == 0) {
        /* Past Fixings Case  tto be clarified with MAD*/

        dCouponValue = sExotic->dPaidMargin[i];
        if (!sAux->iFloatingType[i] == 0) {
          dCouponValue += sExotic->dPrvCpngearing[i] * dPreviousCouponValue +
                          sExotic->dPaidGearing[i] * sExotic->dPastPaidCpn[i];
        }

        dCouponValue = DMAX(dCouponValue, sExotic->dFloor[i]);
        dCouponValue = DMIN(dCouponValue, sExotic->dCap[i]);

        if (sExotic->lPaidFixingDate[i] <
            sMarket->lToday + sMarket->iEODFixFlag) {
          sAux->dFixedPV +=
              dCouponValue * sExotic->dCoverage[i] * dPaidDf * sExotic->dNot[i];
          ;
        }
      } else {
        /* Get the schedule */
        err = add_tenor(sExotic->lStartDate[i], sExotic->cPaidTenor[i],
                        NO_BUSDAY_CONVENTION, &lTheoEndDate);
        if (err)
          goto FREE_RETURN;

        err = swp_f_initSwapDP(sExotic->lStartDate[i], lTheoEndDate,
                               sExotic->cPaidFreq, sExotic->cPaidBasis,
                               &sFixedSwap);

        if (err)
          goto FREE_RETURN;

        sFixedSwapDates = SwapDP_to_DateList(&sFixedSwap, MODIFIED_SUCCEEDING);

        /* Calculate Level and Swap */
        dLevel = 0.0;

        for (j = 1; j < sFixedSwapDates.len; j++) {
          dLevel += coverage(sFixedSwapDates.date[j - 1],
                             sFixedSwapDates.date[j], sFixedSwap.basis_code) *
                    swp_f_df(sMarket->lToday, sFixedSwapDates.date[j],
                             sMarket->cYcName);
        }

        dFloating = swp_f_df(sMarket->lToday, sFixedSwapDates.date[0],
                             sMarket->cYcName) -
                    swp_f_df(sMarket->lToday,
                             sFixedSwapDates.date[sFixedSwapDates.len - 1],
                             sMarket->cYcName);

        dCashSwap = dFloating / dLevel;

        dCouponValue = sExotic->dPaidMargin[i];
        dCouponValue += sExotic->dPrvCpngearing[i] * dPreviousCouponValue +
                        sExotic->dPaidGearing[i] * dCashSwap;
        dCouponValue = DMAX(dCouponValue, sExotic->dFloor[i]);
        dCouponValue = DMIN(dCouponValue, sExotic->dCap[i]);

        // store the CVstrike equal to the forward strike
        sAux->dCVStrike[i] =
            DMAX(0.0001, sExotic->dPaidMargin[i] +
                             sExotic->dPrvCpngearing[i] * dPreviousCouponValue);
      }

      break;
    }
    }

    dPreviousCouponValue = dCouponValue;
    dCouponValue *= sExotic->dCoverage[i] * dPaidDf * sExotic->dNot[i];
    sAux->dMarketValue[i] = dCouponValue;
  }

FREE_RETURN:
  return err;
}

Err crif_calc_mkt_value(CRIF_DEAL sDeal, CMSCTS_MARKET sMarket,
                        CRIF_PRICING_PARAMS sPricingParams,
                        CMSCTS_CALIB sCalibration) {
  Err err = NULL;

  err = crif_calc_fund_value(sDeal, sMarket);
  if (err)
    goto FREE_RETURN;

  err = crif_calc_exotic_value(sDeal, sMarket, sPricingParams, sCalibration);
  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  return err;
}

Err crif_calc_cv_value(CMSCTS_MARKET sMarket, CRIF_DEAL sDeal,
                       CRIF_SIMULARG sSimulArg,
                       CRIF_PRICING_PARAMS sPricingParams,
                       CRIF_OUTPUTS sOutputs) {
  Err err = NULL;
  CRIF_EXOTIC sExotic;
  CRIF_EXOTIC_AUX sAux;
  LGMSV_MODEL sModel;

  double **dTempVal = NULL;
  int i, iStartIndex;
  double dPreviousCouponValue;

  sExotic = sDeal->sExotic;
  sAux = sExotic->sAux;
  sModel = sSimulArg->sModel;

  /* Start Index */
  iStartIndex = sAux->iStartIndex;

  /* Must Value Past Coupons for the Ratchet*/
  dPreviousCouponValue = 0.0;
  for (i = iStartIndex; i <= sExotic->sAux->iEndIndex; i++) {
    /*  do the pricing */
    switch (sExotic->eCouponType[i]) {
    case REGULAR_CRIF: {

      if (sExotic->lPaidFixingDate[i] >
              sMarket->lToday + sMarket->iEODFixFlag &&
          sAux->iFloatingType[i] != 0) {
      }

      break;
    }
    }
  }

  return err;
}