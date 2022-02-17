/* -------------------------------------------------------------------------------------------------------------------
 */
/*
        TargetNoteUtil.c
        Routines:
                TARN_calcHistory:  from the historical fixings  , calculate what
   has happened.
*/
#include "TargetNoteUtil.h"
#include "CPDCalib.h"
#include "TargetNoteCaller.h"
#include "TargetNoteStatic.h"

/* -------------------------------------------------------------------------------------------------------------------
 */
/* TARN Structure */
init_TARN_Struct(TARN_Struct *tarn) {
  init_TARN_Market_Struct(&tarn->market);
  init_TARN_Deal_Struct(&tarn->deal);
  init_TARN_Model_Struct(&tarn->model);
  init_TARN_Calibration_Struct(&tarn->calibration);
  init_TARN_Pricing_Struct(&tarn->pricing);
  init_TARN_Output_Struct(&tarn->output);
}

copy_TARN_Struct(TARN_Struct *in_tarn, TARN_Struct *out_tarn) {
  /* initialize the output struct */
  init_TARN_Struct(out_tarn);

  /* copy each of the member structs */
  copy_TARN_Market_Struct(&in_tarn->market, &out_tarn->market);
  copy_TARN_Deal_Struct(&in_tarn->deal, &out_tarn->deal);
  copy_TARN_Model_Struct(&in_tarn->model, &out_tarn->model);
  copy_TARN_Calibration_Struct(&in_tarn->calibration, &out_tarn->calibration);
  copy_TARN_Pricing_Struct(&in_tarn->pricing, &out_tarn->pricing);
  copy_TARN_Output_Struct(&in_tarn->output, &out_tarn->output);
}

free_TARN_Struct(TARN_Struct *tarn) {
  free_TARN_Market_Struct(&tarn->market);
  free_TARN_Deal_Struct(&tarn->deal);
  free_TARN_Model_Struct(&tarn->model);
  free_TARN_Calibration_Struct(&tarn->calibration);
  free_TARN_Pricing_Struct(&tarn->pricing);
  free_TARN_Output_Struct(&tarn->output);
}
/* TARN Structure */
/* -------------------------------------------------------------------------------------------------------------------
 */

/*
init_TARN_All_Struct(
                                                TARN_Market_Struct*
tarn_market_struct  , TARN_Deal_Struct* tarn_deal_struct  , TARN_Model_Struct*
tarn_model_struct  , TARN_Calibration_Struct* tarn_calibration_struct  ,
                                                TARN_Pricing_Struct*
tarn_pricing_struct  , TARN_Output_Struct* tarn_output_struct
                                        )
{
        init_TARN_Market_Struct( tarn_market_struct );
        init_TARN_Deal_Struct( tarn_deal_struct );
        init_TARN_Model_Struct( tarn_model_struct );
        init_TARN_Calibration_Struct( tarn_calibration_struct );
        init_TARN_Pricing_Struct( tarn_pricing_struct );
        init_TARN_Output_Struct( tarn_output_struct );
}

free_TARN_All_Struct(
                                                TARN_Market_Struct*
tarn_market_struct  , TARN_Deal_Struct* tarn_deal_struct  , TARN_Model_Struct*
tarn_model_struct  , TARN_Calibration_Struct* tarn_calibration_struct  ,
                                                TARN_Pricing_Struct*
tarn_pricing_struct  , TARN_Output_Struct* tarn_output_struct
                                        )
{
        free_TARN_Market_Struct( tarn_market_struct );
        free_TARN_Deal_Struct( tarn_deal_struct );
        free_TARN_Model_Struct( tarn_model_struct );
        free_TARN_Calibration_Struct( tarn_calibration_struct );
        free_TARN_Pricing_Struct( tarn_pricing_struct );
        free_TARN_Output_Struct( tarn_output_struct );
}
*/

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Market Structure */
init_TARN_Market_Struct(TARN_Market_Struct *tarn_market_struct) {
  tarn_market_struct->szYieldCurve = calloc(256, sizeof(char));
  tarn_market_struct->szFundYieldCurve = calloc(256, sizeof(char));
  tarn_market_struct->szVolCurve = calloc(256, sizeof(char));
  tarn_market_struct->szVolCurveRefRate = calloc(256, sizeof(char));
  tarn_market_struct->szCcy = calloc(256, sizeof(char));
  tarn_market_struct->szLGM2FUnd = calloc(256, sizeof(char));
  tarn_market_struct->szLGMSVUnd = calloc(256, sizeof(char));
  tarn_market_struct->szLGM1FUnd = calloc(256, sizeof(char));
  tarn_market_struct->szLGM1FSVUnd = calloc(256, sizeof(char));
  tarn_market_struct->getCashVol = 0;
  tarn_market_struct->getDF = 0;
}

copy_TARN_Market_Struct(TARN_Market_Struct *in, TARN_Market_Struct *out) {
  out->lToday = in->lToday;
  strcpy(out->szYieldCurve, in->szYieldCurve);
  strcpy(out->szFundYieldCurve, in->szFundYieldCurve);
  strcpy(out->szVolCurve, in->szVolCurve);
  strcpy(out->szVolCurveRefRate, in->szVolCurveRefRate);
  strcpy(out->szCcy, in->szCcy);
  out->fx_fund_dom = in->fx_fund_dom;
  out->fx_fund_dom_spot_date = in->fx_fund_dom_spot_date;
  strcpy(out->szLGM2FUnd, in->szLGM2FUnd);
  strcpy(out->szLGMSVUnd, in->szLGMSVUnd);
  strcpy(out->szLGM1FUnd, in->szLGM1FUnd);
  strcpy(out->szLGM1FSVUnd, in->szLGM1FSVUnd);
  out->getCashVol = in->getCashVol;
  out->getDF = in->getDF;
  out->eodFixFlag = in->eodFixFlag;
  out->eodPayFlag = in->eodPayFlag;
}

free_TARN_Market_Struct(TARN_Market_Struct *tarn_market_struct) {
  free_and_zero(&tarn_market_struct->szYieldCurve);
  free_and_zero(&tarn_market_struct->szFundYieldCurve);
  free_and_zero(&tarn_market_struct->szVolCurve);
  free_and_zero(&tarn_market_struct->szVolCurveRefRate);
  free_and_zero(&tarn_market_struct->szCcy);
  free_and_zero(&tarn_market_struct->szLGM2FUnd);
  free_and_zero(&tarn_market_struct->szLGMSVUnd);
  free_and_zero(&tarn_market_struct->szLGM1FUnd);
  free_and_zero(&tarn_market_struct->szLGM1FSVUnd);
}
/* Market Structure */
/* -----------------------------------------------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Deal Structure */
init_TARN_Deal_Struct(TARN_Deal_Struct *tarn_deal_struct) {
  tarn_deal_struct->dvCoupon = 0;
  tarn_deal_struct->lvCpnStartDates = 0;
  tarn_deal_struct->lvCpnPayDates = 0;
  tarn_deal_struct->ivIsCapped = 0;
  tarn_deal_struct->ivIsFloored = 0;
  tarn_deal_struct->dvCap = 0;
  tarn_deal_struct->dvFloor = 0;
  tarn_deal_struct->dvGearing = 0;
  tarn_deal_struct->lvFltrFixingDates = 0;
  tarn_deal_struct->lvFltrStartDates = 0;
  tarn_deal_struct->lvFltrEndDates = 0;
  tarn_deal_struct->dvFltrSpread = 0;
  tarn_deal_struct->dvFltrHistFixings = 0;
  tarn_deal_struct->lvFundFixingDates = 0;
  tarn_deal_struct->lvFundStartDates = 0;
  tarn_deal_struct->lvFundEndDates = 0;
  tarn_deal_struct->dvFundMargin = 0;
  tarn_deal_struct->dvFundSpread = 0;
  tarn_deal_struct->dvFundHistFixings = 0;
  tarn_deal_struct->szBasisFund = calloc(8, sizeof(char));
  tarn_deal_struct->szvBasisFund = 0;
}

copy_TARN_Deal_Struct(TARN_Deal_Struct *in, TARN_Deal_Struct *out) {
  /* Local variables */
  int i;

  /* The target note details */
  out->dTarget = in->dTarget;
  out->couponType = in->couponType;
  out->bRepayCoupon = in->bRepayCoupon;
  out->dCpnNotional = in->dCpnNotional;

  /* The coupon schedule */
  out->nCouponDates = in->nCouponDates;
  out->dvCoupon = calloc(in->nCouponDates, sizeof(double));
  out->lvCpnStartDates = calloc(in->nCouponDates, sizeof(long));
  out->lvCpnPayDates = calloc(in->nCouponDates, sizeof(long));
  out->ivIsFloored = calloc(in->nCouponDates, sizeof(int));
  out->dvFloor = calloc(in->nCouponDates, sizeof(double));
  out->ivIsCapped = calloc(in->nCouponDates, sizeof(int));
  out->dvCap = calloc(in->nCouponDates, sizeof(double));
  out->basisCpn = in->basisCpn;
  out->dFixedCvg = in->dFixedCvg;

  /* The floater */
  out->dvGearing = calloc(in->nCouponDates, sizeof(double));
  out->lvFltrFixingDates = calloc(in->nCouponDates, sizeof(long));
  out->lvFltrStartDates = calloc(in->nCouponDates, sizeof(long));
  out->lvFltrEndDates = calloc(in->nCouponDates, sizeof(long));
  out->dvFltrSpread = calloc(in->nCouponDates, sizeof(double));
  out->basisFltr = in->basisFltr;
  out->dvFltrHistFixings = calloc(in->nCouponDates, sizeof(double));

  /* The funding */
  out->fund_ccy = in->fund_ccy;
  out->dFundingNotional = in->dFundingNotional;
  out->nFundDates = in->nFundDates;
  out->lvFundFixingDates = calloc(in->nFundDates, sizeof(long));
  out->lvFundStartDates = calloc(in->nFundDates, sizeof(long));
  out->lvFundEndDates = calloc(in->nFundDates, sizeof(long));
  out->dvFundMargin = calloc(in->nFundDates, sizeof(double));
  out->dvFundSpread = calloc(in->nFundDates, sizeof(double));
  out->szBasisFund = calloc(16, sizeof(char));
  out->szvBasisFund = calloc(in->nFundDates, sizeof(char *));
  out->basisFund = in->basisFund;
  out->dvFundHistFixings = calloc(in->nFundDates, sizeof(double));

  /* Loop over the each coupon and copy all the values */
  for (i = 0; i < out->nCouponDates; i++) {
    out->dvCoupon[i] = in->dvCoupon[i];
    out->lvCpnStartDates[i] = in->lvCpnStartDates[i];
    out->lvCpnPayDates[i] = in->lvCpnPayDates[i];
    out->ivIsFloored[i] = in->ivIsFloored[i];
    out->dvFloor[i] = in->dvFloor[i];
    out->ivIsCapped[i] = in->ivIsCapped[i];
    out->dvCap[i] = in->dvCap[i];
    out->dvGearing[i] = in->dvGearing[i];
    out->lvFltrFixingDates[i] = in->lvFltrFixingDates[i];
    out->lvFltrStartDates[i] = in->lvFltrStartDates[i];
    out->lvFltrEndDates[i] = in->lvFltrEndDates[i];
    out->dvFltrSpread[i] = in->dvFltrSpread[i];
    out->dvFltrHistFixings[i] = in->dvFltrHistFixings[i];
  }

  /* Loop over each funding date and copy all the values */
  strcpy(out->szBasisFund, in->szBasisFund);
  for (i = 0; i < out->nFundDates; i++) {
    out->lvFundFixingDates[i] = in->lvFundFixingDates[i];
    out->lvFundStartDates[i] = in->lvFundStartDates[i];
    out->lvFundEndDates[i] = in->lvFundEndDates[i];
    out->dvFundMargin[i] = in->dvFundMargin[i];
    out->dvFundSpread[i] = in->dvFundSpread[i];
    out->szvBasisFund[i] = calloc(16, sizeof(char));
    strcpy(out->szvBasisFund[i], in->szvBasisFund[i]);
    out->dvFundHistFixings[i] = in->dvFundHistFixings[i];
  }
}

free_TARN_Deal_Struct(TARN_Deal_Struct *tarn_deal_struct) {
  int i;
  free_and_zero(&tarn_deal_struct->dvCoupon);
  free_and_zero(&tarn_deal_struct->lvCpnStartDates);
  free_and_zero(&tarn_deal_struct->lvCpnPayDates);
  free_and_zero(&tarn_deal_struct->ivIsCapped);
  free_and_zero(&tarn_deal_struct->ivIsFloored);
  free_and_zero(&tarn_deal_struct->dvCap);
  free_and_zero(&tarn_deal_struct->dvFloor);
  free_and_zero(&tarn_deal_struct->dvGearing);
  free_and_zero(&tarn_deal_struct->lvFltrFixingDates);
  free_and_zero(&tarn_deal_struct->lvFltrStartDates);
  free_and_zero(&tarn_deal_struct->lvFltrEndDates);
  free_and_zero(&tarn_deal_struct->dvFltrSpread);
  free_and_zero(&tarn_deal_struct->dvFltrHistFixings);
  free_and_zero(&tarn_deal_struct->lvFundFixingDates);
  free_and_zero(&tarn_deal_struct->lvFundStartDates);
  free_and_zero(&tarn_deal_struct->lvFundEndDates);
  free_and_zero(&tarn_deal_struct->dvFundMargin);
  free_and_zero(&tarn_deal_struct->dvFundSpread);
  free_and_zero(&tarn_deal_struct->dvFundHistFixings);
  free_and_zero(&tarn_deal_struct->szBasisFund);
  for (i = 0; i < tarn_deal_struct->nFundDates; i++)
    free_and_zero(&(tarn_deal_struct->szvBasisFund[i]));
  //	free_and_zero( &( (void*) tarn_deal_struct->szvBasisFund ) );
}
/* Deal Structure */
/* -----------------------------------------------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Model Structure */
init_TARN_Model_Struct(TARN_Model_Struct *tarn_model_struct) {
  tarn_model_struct->lam_time = (double *)calloc(1, sizeof(double));
  tarn_model_struct->lam_time[0] = 10.0;
  tarn_model_struct->lam = (double *)calloc(1, sizeof(double));
  tarn_model_struct->nlam = 1;
  tarn_model_struct->smilepartime = 0;
  tarn_model_struct->alphaepsts = 0;
  tarn_model_struct->ldaepsts = 0;
  tarn_model_struct->rhoepsts = 0;
  tarn_model_struct->rho2epsts = 0;
  tarn_model_struct->smilepartime1F = 0;
  tarn_model_struct->alphaepsts1F = 0;
  tarn_model_struct->ldaepsts1F = 0;
  tarn_model_struct->rhoepsts1F = 0;
}

copy_TARN_Model_Struct(TARN_Model_Struct *in, TARN_Model_Struct *out) {
  /* local variables */
  int i;

  /* LGM2F model parameters */
  out->alpha = in->alpha;
  out->gamma = in->gamma;
  out->rho = in->rho;
  out->nlam = in->nlam;
  out->lam_time = calloc(in->nlam, sizeof(double));
  out->lam = calloc(in->nlam, sizeof(double));
  for (i = 0; i < out->nlam; i++) {
    out->lam_time[i] = in->lam_time[i];
    out->lam[i] = in->lam[i];
  }

  /* SV Model Parameters */
  out->nsmilepar = in->nsmilepar;
  out->smilepartime = calloc(in->nsmilepar, sizeof(double));
  out->alphaepsts = calloc(in->nsmilepar, sizeof(double));
  out->ldaepsts = calloc(in->nsmilepar, sizeof(double));
  out->rhoepsts = calloc(in->nsmilepar, sizeof(double));
  out->rho2epsts = calloc(in->nsmilepar, sizeof(double));
  out->tstar = in->tstar;
  for (i = 0; i < out->nsmilepar; i++) {
    out->smilepartime[i] = in->smilepartime[i];
    out->alphaepsts[i] = in->alphaepsts[i];
    out->ldaepsts[i] = in->ldaepsts[i];
    out->rhoepsts[i] = in->rhoepsts[i];
    out->rho2epsts[i] = in->rho2epsts[i];
  }

  /* 1FSV Parameters */
  out->nsmilepar1F = in->nsmilepar1F;
  out->smilepartime1F = calloc(in->nsmilepar1F, sizeof(double));
  out->alphaepsts1F = calloc(in->nsmilepar1F, sizeof(double));
  out->ldaepsts1F = calloc(in->nsmilepar1F, sizeof(double));
  out->rhoepsts1F = calloc(in->nsmilepar1F, sizeof(double));
  for (i = 0; i < out->nsmilepar1F; i++) {
    out->smilepartime1F[i] = in->smilepartime1F[i];
    out->alphaepsts1F[i] = in->alphaepsts1F[i];
    out->ldaepsts1F[i] = in->ldaepsts1F[i];
    out->rhoepsts1F[i] = in->rhoepsts1F[i];
  }
}

free_TARN_Model_Struct(TARN_Model_Struct *tarn_model_struct) {
  free_and_zero(&tarn_model_struct->lam_time);
  free_and_zero(&tarn_model_struct->lam);
  free_and_zero(&tarn_model_struct->smilepartime);
  free_and_zero(&tarn_model_struct->alphaepsts);
  free_and_zero(&tarn_model_struct->ldaepsts);
  free_and_zero(&tarn_model_struct->rhoepsts);
  free_and_zero(&tarn_model_struct->rho2epsts);
  free_and_zero(&tarn_model_struct->alphaepsts1F);
  free_and_zero(&tarn_model_struct->ldaepsts1F);
  free_and_zero(&tarn_model_struct->rhoepsts1F);
}
/* Model Structure */
/* -----------------------------------------------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Calibration Structure */
init_TARN_Calibration_Struct(TARN_Calibration_Struct *tarn_calibration_struct) {
  tarn_calibration_struct->szPrimRefRate = calloc(256, sizeof(char));
  tarn_calibration_struct->szPrimFreq = calloc(256, sizeof(char));
  tarn_calibration_struct->szPrimBasis = calloc(256, sizeof(char));
  tarn_calibration_struct->szPrimTenor = calloc(256, sizeof(char));
  tarn_calibration_struct->szSecRefRate = calloc(256, sizeof(char));
  tarn_calibration_struct->szSecFreq = calloc(256, sizeof(char));
  tarn_calibration_struct->szSecBasis = calloc(256, sizeof(char));
  tarn_calibration_struct->szSecTenor = calloc(256, sizeof(char));
}

copy_TARN_Calibration_Struct(TARN_Calibration_Struct *in,
                             TARN_Calibration_Struct *out) {
  /* Instruments calibration options */
  out->skip_start = in->skip_start;
  out->num_prob_paths = in->num_prob_paths;

  memcpy(&out->prim_param, &in->prim_param, sizeof(cpd_diag_calib_param));
  memcpy(&out->sec_param, &in->sec_param, sizeof(cpd_diag_calib_param));

  /* Volatility calibration options */
  out->nb_iter_LM = in->nb_iter_LM;
  out->precision_LM = in->precision_LM;
  /* LGM2F parameter options */
  out->fix_lambda = in->fix_lambda;
  /* LGM2V calibration options */
  LGMSV_Copy_CalibParams(&in->lgmsv_calib_params, &out->lgmsv_calib_params);
  /* The calibration instrument parameters */
  strcpy(out->szPrimRefRate, in->szPrimRefRate);
  strcpy(out->szPrimFreq, in->szPrimFreq);
  strcpy(out->szPrimBasis, in->szPrimBasis);
  strcpy(out->szPrimTenor, in->szPrimTenor);
  strcpy(out->szSecRefRate, in->szSecRefRate);
  strcpy(out->szSecFreq, in->szSecFreq);
  strcpy(out->szSecBasis, in->szSecBasis);
  strcpy(out->szSecTenor, in->szSecTenor);
  /*	SV Numerical params */
  out->iNbX = in->iNbX;
  out->iNbSigmaXLeft = in->iNbSigmaXLeft;
  out->iNbSigmaXRight = in->iNbSigmaXRight;
  out->dIntegParam = in->dIntegParam;
  out->iIntegMethod = in->iIntegMethod;
  out->dVolLimit = in->dVolLimit;
  out->iCalibLGM = in->iCalibLGM;
  out->dMinStd = in->dMinStd;
  out->dMaxStd = in->dMaxStd;
  out->numer_tstar = in->numer_tstar;
}

free_TARN_Calibration_Struct(TARN_Calibration_Struct *tarn_calibration_struct) {
  free_and_zero(&tarn_calibration_struct->szPrimRefRate);
  free_and_zero(&tarn_calibration_struct->szPrimFreq);
  free_and_zero(&tarn_calibration_struct->szPrimBasis);
  free_and_zero(&tarn_calibration_struct->szPrimTenor);
  free_and_zero(&tarn_calibration_struct->szSecRefRate);
  free_and_zero(&tarn_calibration_struct->szSecFreq);
  free_and_zero(&tarn_calibration_struct->szSecBasis);
  free_and_zero(&tarn_calibration_struct->szSecTenor);
}
/* Calibration Structure */
/* -----------------------------------------------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Pricing Structure */
init_TARN_Pricing_Struct(TARN_Pricing_Struct *tarn_pricing_struct) {}

copy_TARN_Pricing_Struct(TARN_Pricing_Struct *in, TARN_Pricing_Struct *out) {
  /* pricings required */
  out->iPrice2FCV = in->iPrice2FCV;
  out->iPriceSV = in->iPriceSV;
  /* MC */
  out->dSpread = in->dSpread;
  out->iSpreadType = in->iSpreadType;
  out->nStepT = in->nStepT;
  out->num_paths = in->num_paths;
  out->do_pecs = in->do_pecs;
}

free_TARN_Pricing_Struct(TARN_Pricing_Struct *tarn_pricing_struct) {}
/* Pricing Structure */
/* -----------------------------------------------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Output Structure */
init_TARN_Output_Struct(TARN_Output_Struct *tarn_output_struct) {
  /* return matrices */
  int row = 2, col = 2, i;
  tarn_output_struct->dmLGM2F_PV = (double **)calloc(row, sizeof(double *));
  tarn_output_struct->dmLGMSV_PV = (double **)calloc(row, sizeof(double *));
  tarn_output_struct->dmLGM1F_PV = (double **)calloc(row, sizeof(double *));
  tarn_output_struct->dmLGM1FSV_PV = (double **)calloc(row, sizeof(double *));
  for (i = 0; i < row; i++) {
    tarn_output_struct->dmLGM2F_PV[i] = (double *)calloc(col, sizeof(double));
    tarn_output_struct->dmLGMSV_PV[i] = (double *)calloc(col, sizeof(double));
    tarn_output_struct->dmLGM1F_PV[i] = (double *)calloc(col, sizeof(double));
    tarn_output_struct->dmLGM1FSV_PV[i] = (double *)calloc(col, sizeof(double));
  }

  /* set the size of nProb */
  tarn_output_struct->nProb = 0;
  tarn_output_struct->dmKO_prob = 0;

  /* underlying names */
  tarn_output_struct->szTARN_LGM2F_UND = calloc(256, sizeof(char));
  strcpy(tarn_output_struct->szTARN_LGM2F_UND, "TARN_LGM2F_UND");

  tarn_output_struct->szTARN_LGMSV_UND = calloc(256, sizeof(char));
  strcpy(tarn_output_struct->szTARN_LGMSV_UND, "TARN_LGMSV_UND");

  tarn_output_struct->szTARN_LGM1F_UND = calloc(256, sizeof(char));
  strcpy(tarn_output_struct->szTARN_LGM1F_UND, "TARN_LGM1F_UND");

  tarn_output_struct->szTARN_LGM1FSV_UND = calloc(256, sizeof(char));
  strcpy(tarn_output_struct->szTARN_LGM1FSV_UND, "TARN_LGM1FSV_UND");

  /* output instruments */
  cpd_init_calib_inst_data(&tarn_output_struct->inst_data_lgm1F);
  cpd_init_calib_inst_data(&tarn_output_struct->inst_data_lgm2F);
  cpd_init_calib_inst_data(&tarn_output_struct->inst_data_lgm1FSV);
  cpd_init_calib_inst_data(&tarn_output_struct->inst_data_lgmSV);
}

copy_TARN_Output_Struct(TARN_Output_Struct *in, TARN_Output_Struct *out) {
  /* local variables */
  int row = 2, col = 2, i, j;

  /* Output:  PV information */
  out->iIsKnockedOut = in->iIsKnockedOut;
  out->iIsExpired = in->iIsExpired;
  out->dCumulCoupon = in->dCumulCoupon;
  for (i = 0; i < row; i++)
    for (j = 0; j < col; j++) {
      out->dmLGM2F_PV[i][j] = in->dmLGM2F_PV[i][j];
      out->dmLGM1F_PV[i][j] = in->dmLGM1F_PV[i][j];
      out->dmLGM1FSV_PV[i][j] = in->dmLGM1FSV_PV[i][j];
      out->dmLGMSV_PV[i][j] = in->dmLGMSV_PV[i][j];
    }

  /* KO info */
  out->nProb = in->nProb;
  if (in->nProb > 0) {
    out->dmKO_prob = dmatrix(0, in->nProb - 1, 0, 1);
    for (i = 0; i < in->nProb; i++) {
      out->dmKO_prob[i][0] = in->dmKO_prob[i][0];
      out->dmKO_prob[i][1] = in->dmKO_prob[i][1];
    }
  } else
    out->dmKO_prob = 0;

  /* Output:  calibrated underlyings */
  strcpy(out->szTARN_LGM2F_UND, in->szTARN_LGM2F_UND);
  strcpy(out->szTARN_LGM1F_UND, in->szTARN_LGM1F_UND);
  strcpy(out->szTARN_LGM1FSV_UND, in->szTARN_LGM1FSV_UND);
  strcpy(out->szTARN_LGMSV_UND, in->szTARN_LGMSV_UND);

  /* Output:  calbiration instrument data */
  cpd_copy_calib_inst_data(&out->inst_data_lgm2F, &in->inst_data_lgm2F);
  cpd_copy_calib_inst_data(&out->inst_data_lgm1F, &in->inst_data_lgm1F);
  cpd_copy_calib_inst_data(&out->inst_data_lgm1FSV, &in->inst_data_lgm1FSV);
  cpd_copy_calib_inst_data(&out->inst_data_lgmSV, &in->inst_data_lgmSV);
}

free_TARN_Output_Struct(TARN_Output_Struct *tarn_output_struct) {
  int row = 2, col = 2, i;

  if (tarn_output_struct->dmLGM2F_PV) {
    for (i = 0; i < row; i++)
      free_and_zero(&tarn_output_struct->dmLGM2F_PV[i]);
    free_and_zero((void *)&tarn_output_struct->dmLGM2F_PV);
  }
  free_and_zero(&tarn_output_struct->szTARN_LGM2F_UND);

  if (tarn_output_struct->dmLGM1F_PV) {
    for (i = 0; i < row; i++)
      free_and_zero(&tarn_output_struct->dmLGM1F_PV[i]);
    free_and_zero((void *)&tarn_output_struct->dmLGM1F_PV);
  }
  free_and_zero(&tarn_output_struct->szTARN_LGM1F_UND);

  if (tarn_output_struct->dmLGM1FSV_PV) {
    for (i = 0; i < row; i++)
      free_and_zero(&tarn_output_struct->dmLGM1FSV_PV[i]);
    free_and_zero((void *)&tarn_output_struct->dmLGM1FSV_PV);
  }
  free_and_zero(&tarn_output_struct->szTARN_LGM1FSV_UND);

  if (tarn_output_struct->dmLGMSV_PV) {
    for (i = 0; i < row; i++)
      free_and_zero(&tarn_output_struct->dmLGMSV_PV[i]);
    free_and_zero((void *)&tarn_output_struct->dmLGMSV_PV);
  }
  free_and_zero(&tarn_output_struct->szTARN_LGMSV_UND);

  if (tarn_output_struct->dmKO_prob)
    free_dmatrix(tarn_output_struct->dmKO_prob, 0, tarn_output_struct->nProb, 0,
                 1);
  tarn_output_struct->dmKO_prob = 0;
  tarn_output_struct->nProb = 0;

  cpd_free_calib_inst_data(&tarn_output_struct->inst_data_lgm2F);
  cpd_free_calib_inst_data(&tarn_output_struct->inst_data_lgm1F);
  cpd_free_calib_inst_data(&tarn_output_struct->inst_data_lgm1FSV);
  cpd_free_calib_inst_data(&tarn_output_struct->inst_data_lgmSV);
}
/* Output Structure */
/* -----------------------------------------------------------------------------------------------------------------
 */

/* Return values
                                                int* out_i1stCpn  ,
                                                int* out_iIsKnockedOut  ,
                                                double* out_dCouponPV  ,
                                                double* out_dFundPV  ,
                                                double* out_dCumulCoupon )
*/

char *TARN_calcHistoryOld(TARN_Struct *tarn, TARN_AUX *aux) {
  /* Variable declaration */
  long lFirstFixing = tarn->market.lToday + !tarn->market.eodFixFlag;
  long lFirstPayment = tarn->market.lToday + !tarn->market.eodPayFlag;
  char *err = 0;
  int i = 0;
  double dCoupon = 0.0;
  int iTemp;

  /* initialize the return variables */
  tarn->output.iIsKnockedOut = 0;
  tarn->output.iIsExpired = 0;
  tarn->output.dCumulCoupon = 0.0;

  /* Find the value of the cumulative coupon and see if it has knocked out.  If
     it hasn't knocked out  , find the first event date */
  while ((!tarn->output.iIsKnockedOut) &&
         (aux->i1stCpn < tarn->deal.nCouponDates) &&
         (tarn->deal.lvFltrFixingDates[aux->i1stCpn] < lFirstFixing)) {
    dCoupon = tarn->deal.dvCoupon[aux->i1stCpn] +
              tarn->deal.dvGearing[aux->i1stCpn] *
                  tarn->deal.dvFltrHistFixings[aux->i1stCpn];
    if (tarn->deal.ivIsFloored[aux->i1stCpn])
      dCoupon = dCoupon < tarn->deal.dvFloor[aux->i1stCpn]
                    ? tarn->deal.dvFloor[aux->i1stCpn]
                    : dCoupon;
    if (tarn->deal.ivIsCapped[aux->i1stCpn])
      dCoupon = dCoupon > tarn->deal.dvCap[aux->i1stCpn]
                    ? tarn->deal.dvCap[aux->i1stCpn]
                    : dCoupon;
    if (tarn->output.dCumulCoupon + dCoupon * aux->dvCumulCvg[i] >=
        tarn->deal.dTarget) {
      tarn->output.iIsKnockedOut = aux->i1stCpn + 1;
      dCoupon = (tarn->deal.dTarget - tarn->output.dCumulCoupon) /
                aux->dvCouponCvg[aux->i1stCpn];
      tarn->output.dCumulCoupon = tarn->deal.dTarget;

    } else {
      tarn->output.dCumulCoupon += dCoupon * aux->dvCumulCvg[i];
      ++aux->i1stCpn;
    }
  }

  /* if there is only 1 coupon remaining and we have to repay the notional  , */
  /* then the cash flows are deterministic  , so it is as if we have knocked out
   */
  if (aux->i1stCpn == (tarn->deal.nCouponDates - 1) &&
      tarn->deal.bRepayCoupon) {
    tarn->output.iIsKnockedOut = tarn->deal.nCouponDates;
    dCoupon = (tarn->deal.dTarget - tarn->output.dCumulCoupon) /
              aux->dvCouponCvg[aux->i1stCpn];
  }

  /* If the swap has knocked out  , see if there are any remaining cash flows
   * and calculate their PV */
  if (tarn->output.iIsKnockedOut) {
    /* Calculate the coupon PV */
    if (tarn->deal.couponType == TN_SWAP) {
      if (tarn->deal.lvCpnPayDates[aux->i1stCpn] >= lFirstPayment)
        aux->dHistCpnPV = dCoupon * aux->dvCouponCvg[aux->i1stCpn] *
                          aux->dvCouponPayDF[aux->i1stCpn];
      //			aux->dHistCpnPV += ( tarn->deal.dTarget -
      //tarn->output.dCumulCoupon ) * aux->dvCouponPayDF[aux->i1stCpn];
    }
    if (tarn->deal.couponType == TN_SWAP) {
      if (tarn->deal.lvCpnPayDates[aux->i1stCpn] >= lFirstPayment)
        aux->dHistCpnPV = dCoupon * aux->dvCouponCvg[aux->i1stCpn] *
                          aux->dvCouponPayDF[aux->i1stCpn];
      //			aux->dHistCpnPV += ( tarn->deal.dTarget -
      //tarn->output.dCumulCoupon ) * aux->dvCouponPayDF[aux->i1stCpn];
    } else if (tarn->deal.couponType == TN_ZC) {
      aux->dHistCpnPV = tarn->deal.dTarget * aux->dvCouponPayDF[aux->i1stCpn];
    } else if (tarn->deal.couponType == TN_FINAL) {
      if (tarn->deal.lvCpnPayDates[aux->i1stCpn] >= lFirstPayment)
        aux->dHistCpnPV = dCoupon * aux->dvCouponCvg[aux->i1stCpn] *
                          aux->dvCouponPayDF[aux->i1stCpn];
      aux->dHistCpnPV += (tarn->deal.dTarget - tarn->output.dCumulCoupon) *
                         aux->dvCouponPayDF[aux->i1stCpn];
    }

    /* Calculate the funding PV */
    for (i = aux->ivFundStartIndex[aux->i1stCpn];
         i <= aux->ivFundEndIndex[aux->i1stCpn]; i++) {
      if (tarn->deal.lvFundFixingDates[i] < lFirstFixing &&
          tarn->deal.lvFundEndDates[i] >= lFirstPayment)
        aux->dHistFundPV +=
            (tarn->deal.dvFundHistFixings[i] + tarn->deal.dvFundMargin[i]) *
            aux->dvFundCvg[i] * aux->dvFundDomPayDF[i];

      if (tarn->deal.lvFundFixingDates[i] >= lFirstFixing &&
          tarn->deal.lvFundEndDates[i] >= lFirstPayment) {
        double dFund_Start_DF =
            swp_f_df(tarn->market.lToday, tarn->deal.lvFundStartDates[i],
                     tarn->market.szYieldCurve);
        double dFund_End_DF =
            swp_f_df(tarn->market.lToday, tarn->deal.lvFundEndDates[i],
                     tarn->market.szYieldCurve);
        double dFund_cvgMargin =
            1.0 - coverage(tarn->deal.lvFundStartDates[i],
                           tarn->deal.lvFundEndDates[i], tarn->deal.basisFund) *
                      (tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i]);
        aux->dHistFundPV += dFund_Start_DF - dFund_cvgMargin * dFund_End_DF;
      }
    }
    return 0;
  }

  /* if the swap has not knocked out  , calculate any residual coupon cash flows
   */
  if (tarn->deal.couponType == TN_SWAP) {
    if ((aux->i1stCpn > 0) &&
        (tarn->deal.lvCpnPayDates[aux->i1stCpn - 1] >= lFirstPayment))
      aux->dHistCpnPV = (tarn->deal.dvCoupon[aux->i1stCpn - 1] +
                         tarn->deal.dvGearing[aux->i1stCpn - 1] *
                             tarn->deal.dvFltrHistFixings[aux->i1stCpn - 1]) *
                        aux->dvCouponCvg[aux->i1stCpn - 1] *
                        aux->dvCouponPayDF[aux->i1stCpn - 1];
  }

  /* if the swap has not knocked out  , calculate any residual funding PV up to
   * the first fixing date */
  iTemp = aux->i1stCpn == 0 ? 0 : aux->i1stCpn - 1;
  for (i = aux->ivFundStartIndex[iTemp]; i <= aux->ivFundEndIndex[aux->i1stCpn];
       i++) {
    if (tarn->deal.lvFundFixingDates[i] < lFirstFixing &&
        tarn->deal.lvFundEndDates[i] >= lFirstPayment)
      aux->dHistFundPV +=
          (tarn->deal.dvFundHistFixings[i] + tarn->deal.dvFundMargin[i]) *
          aux->dvFundCvg[i] * aux->dvFundDomPayDF[i];

    if (tarn->deal.lvFundFixingDates[i] >= lFirstFixing &&
        tarn->deal.lvFundEndDates[i] >= lFirstPayment) {
      double dFund_Start_DF =
          swp_f_df(tarn->market.lToday, tarn->deal.lvFundStartDates[i],
                   tarn->market.szYieldCurve);
      double dFund_End_DF =
          swp_f_df(tarn->market.lToday, tarn->deal.lvFundEndDates[i],
                   tarn->market.szYieldCurve);
      double dFund_cvgMargin =
          1.0 - coverage(tarn->deal.lvFundStartDates[i],
                         tarn->deal.lvFundEndDates[i], tarn->deal.basisFund) *
                    (tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i]);
      aux->dHistFundPV += dFund_Start_DF - dFund_cvgMargin * dFund_End_DF;
    }
  }

  /* if there are no more fixings left  , then mark it as finished */
  if (aux->i1stCpn == tarn->deal.nCouponDates)
    tarn->output.iIsExpired = 1;
  else
    tarn->output.iIsExpired = 0;

  /* return */
  return 0;
}

/*
static double TARN_calcCoupon( TARN_Struct* tarn  , TARN_AUX* aux )
{
        double dCoupon = tarn->deal.dvCoupon[aux->i1stCpn] +
tarn->deal.dvGearing[aux->i1stCpn] * tarn->deal.dvFltrHistFixings[aux->i1stCpn];
        if ( tarn->deal.ivIsFloored[aux->i1stCpn] )
                dCoupon = dCoupon < tarn->deal.dvFloor[aux->i1stCpn] ?
tarn->deal.dvFloor[aux->i1stCpn] : dCoupon; if (
tarn->deal.ivIsCapped[aux->i1stCpn] ) dCoupon = dCoupon >
tarn->deal.dvCap[aux->i1stCpn] ? tarn->deal.dvCap[aux->i1stCpn] : dCoupon;
}
*/

char *TARN_calcHistory_ZC(TARN_Struct *tarn, TARN_AUX *aux) {
  /* Variable declaration */
  long lFirstFixing = tarn->market.lToday + !tarn->market.eodFixFlag;
  long lFirstPayment = tarn->market.lToday + !tarn->market.eodPayFlag;
  char *err = 0;
  int i = 0;
  double dCoupon = 0.0;
  int iTemp;

  /* Find the value of the cumulative coupon and see if it has knocked out.  If
     it hasn't knocked out  , find the first event date */
  while ((!tarn->output.iIsKnockedOut) &&
         (aux->i1stCpn < tarn->deal.nCouponDates) &&
         (tarn->deal.lvFltrFixingDates[aux->i1stCpn] < lFirstFixing)) {
    dCoupon = tarn->deal.dvCoupon[aux->i1stCpn] +
              tarn->deal.dvGearing[aux->i1stCpn] *
                  tarn->deal.dvFltrHistFixings[aux->i1stCpn];
    if (tarn->deal.ivIsFloored[aux->i1stCpn])
      dCoupon = dCoupon < tarn->deal.dvFloor[aux->i1stCpn]
                    ? tarn->deal.dvFloor[aux->i1stCpn]
                    : dCoupon;
    if (tarn->deal.ivIsCapped[aux->i1stCpn])
      dCoupon = dCoupon > tarn->deal.dvCap[aux->i1stCpn]
                    ? tarn->deal.dvCap[aux->i1stCpn]
                    : dCoupon;
    if (tarn->output.dCumulCoupon + dCoupon * aux->dvCumulCvg[i] >=
        tarn->deal.dTarget) {
      tarn->output.iIsKnockedOut = aux->i1stCpn + 1;
      dCoupon = (tarn->deal.dTarget - tarn->output.dCumulCoupon) /
                aux->dvCouponCvg[aux->i1stCpn];
      tarn->output.dCumulCoupon = tarn->deal.dTarget;
    } else {
      tarn->output.dCumulCoupon += dCoupon * aux->dvCumulCvg[i];
      ++aux->i1stCpn;
    }
  }

  /* if there is only 1 coupon remaining  , then the cash flows are
   * deterministic  , so it is as if we have knocked out */
  if (aux->i1stCpn == (tarn->deal.nCouponDates - 1) &&
      tarn->deal.bRepayCoupon) {
    tarn->output.iIsKnockedOut = tarn->deal.nCouponDates;
  }

  /* If the swap has knocked out  , see if there are any remaining cash flows
   * and calculate their PV */
  if (tarn->output.iIsKnockedOut) {
    /* Calculate the coupon PV */
    if (tarn->deal.lvCpnPayDates[aux->i1stCpn] >= lFirstPayment)
      aux->dHistCpnPV = tarn->deal.dTarget * aux->dvCouponPayDF[aux->i1stCpn];

    /* Calculate the funding PV */
    for (i = aux->ivFundStartIndex[aux->i1stCpn];
         i <= aux->ivFundEndIndex[aux->i1stCpn]; i++) {
      if (tarn->deal.lvFundFixingDates[i] < lFirstFixing &&
          tarn->deal.lvFundEndDates[i] >= lFirstPayment)
        aux->dHistFundPV +=
            (tarn->deal.dvFundHistFixings[i] + tarn->deal.dvFundMargin[i]) *
            aux->dvFundCvg[i] * aux->dvFundDomPayDF[i];

      if (tarn->deal.lvFundFixingDates[i] >= lFirstFixing &&
          tarn->deal.lvFundEndDates[i] >= lFirstPayment) {
        double dFund_Start_DF =
            swp_f_df(tarn->market.lToday, tarn->deal.lvFundStartDates[i],
                     tarn->market.szYieldCurve);
        double dFund_End_DF =
            swp_f_df(tarn->market.lToday, tarn->deal.lvFundEndDates[i],
                     tarn->market.szYieldCurve);
        double dFund_cvgMargin =
            1.0 - coverage(tarn->deal.lvFundStartDates[i],
                           tarn->deal.lvFundEndDates[i], tarn->deal.basisFund) *
                      (tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i]);
        aux->dHistFundPV += dFund_Start_DF - dFund_cvgMargin * dFund_End_DF;
      }
    }
    return 0;
  }

  /* if the swap has not knocked out  , calculate any residual funding PV up to
   * the first fixing date */
  iTemp = aux->i1stCpn == 0 ? 0 : aux->i1stCpn - 1;
  for (i = aux->ivFundStartIndex[iTemp]; i <= aux->ivFundEndIndex[aux->i1stCpn];
       i++) {
    if (tarn->deal.lvFundFixingDates[i] < lFirstFixing &&
        tarn->deal.lvFundEndDates[i] >= lFirstPayment)
      aux->dHistFundPV +=
          (tarn->deal.dvFundHistFixings[i] + tarn->deal.dvFundMargin[i]) *
          aux->dvFundCvg[i] * aux->dvFundDomPayDF[i];

    if (tarn->deal.lvFundFixingDates[i] >= lFirstFixing &&
        tarn->deal.lvFundEndDates[i] >= lFirstPayment) {
      double dFund_Start_DF =
          swp_f_df(tarn->market.lToday, tarn->deal.lvFundStartDates[i],
                   tarn->market.szYieldCurve);
      double dFund_End_DF =
          swp_f_df(tarn->market.lToday, tarn->deal.lvFundEndDates[i],
                   tarn->market.szYieldCurve);
      double dFund_cvgMargin =
          1.0 - coverage(tarn->deal.lvFundStartDates[i],
                         tarn->deal.lvFundEndDates[i], tarn->deal.basisFund) *
                    (tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i]);
      aux->dHistFundPV += dFund_Start_DF - dFund_cvgMargin * dFund_End_DF;
    }
  }

  /* return */
  return 0;
}

char *TARN_calcHistory_SWAP(TARN_Struct *tarn, TARN_AUX *aux) {
  /* Variable declaration */
  long lFirstFixing = tarn->market.lToday + !tarn->market.eodFixFlag;
  long lFirstPayment = tarn->market.lToday + !tarn->market.eodPayFlag;
  char *err = 0;
  int i = 0;
  double dCoupon = 0.0;
  int iTemp;

  /* Find the value of the cumulative coupon and see if it has knocked out.  If
     it hasn't knocked out  , find the first event date */
  while ((!tarn->output.iIsKnockedOut) &&
         (aux->i1stCpn < tarn->deal.nCouponDates) &&
         (tarn->deal.lvFltrFixingDates[aux->i1stCpn] < lFirstFixing)) {
    dCoupon = tarn->deal.dvCoupon[aux->i1stCpn] +
              tarn->deal.dvGearing[aux->i1stCpn] *
                  tarn->deal.dvFltrHistFixings[aux->i1stCpn];
    if (tarn->deal.ivIsFloored[aux->i1stCpn])
      dCoupon = dCoupon < tarn->deal.dvFloor[aux->i1stCpn]
                    ? tarn->deal.dvFloor[aux->i1stCpn]
                    : dCoupon;
    if (tarn->deal.ivIsCapped[aux->i1stCpn])
      dCoupon = dCoupon > tarn->deal.dvCap[aux->i1stCpn]
                    ? tarn->deal.dvCap[aux->i1stCpn]
                    : dCoupon;
    if (tarn->output.dCumulCoupon + dCoupon * aux->dvCumulCvg[i] >=
        tarn->deal.dTarget) {
      tarn->output.iIsKnockedOut = aux->i1stCpn + 1;
      dCoupon = (tarn->deal.dTarget - tarn->output.dCumulCoupon) /
                aux->dvCouponCvg[aux->i1stCpn];
      tarn->output.dCumulCoupon = tarn->deal.dTarget;
    } else {
      tarn->output.dCumulCoupon += dCoupon * aux->dvCumulCvg[i];
      ++aux->i1stCpn;
    }
  }

  /* if there is only 1 coupon remaining and the coupon is guaranteed  , then
   * the cash flows are deterministic  , */
  /* so it is as if we have knocked out */
  if (aux->i1stCpn == (tarn->deal.nCouponDates - 1) &&
      tarn->deal.bRepayCoupon) {
    tarn->output.iIsKnockedOut = tarn->deal.nCouponDates;
    dCoupon = (tarn->deal.dTarget - tarn->output.dCumulCoupon) /
              aux->dvCouponCvg[aux->i1stCpn];
  }

  /* If the swap has knocked out  , see if there are any remaining cash flows
   * and calculate their PV */
  if (tarn->output.iIsKnockedOut) {
    /* Calculate the coupon PV */
    if (tarn->deal.lvCpnPayDates[aux->i1stCpn] >= lFirstPayment)
      aux->dHistCpnPV = dCoupon * aux->dvCouponCvg[aux->i1stCpn] *
                        aux->dvCouponPayDF[aux->i1stCpn];

    /* Calculate the funding PV */
    for (i = aux->ivFundStartIndex[aux->i1stCpn];
         i <= aux->ivFundEndIndex[aux->i1stCpn]; i++) {
      if (tarn->deal.lvFundFixingDates[i] < lFirstFixing &&
          tarn->deal.lvFundEndDates[i] >= lFirstPayment)
        aux->dHistFundPV +=
            (tarn->deal.dvFundHistFixings[i] + tarn->deal.dvFundMargin[i]) *
            aux->dvFundCvg[i] * aux->dvFundDomPayDF[i];

      if (tarn->deal.lvFundFixingDates[i] >= lFirstFixing &&
          tarn->deal.lvFundEndDates[i] >= lFirstPayment) {
        double dFund_Start_DF =
            swp_f_df(tarn->market.lToday, tarn->deal.lvFundStartDates[i],
                     tarn->market.szYieldCurve);
        double dFund_End_DF =
            swp_f_df(tarn->market.lToday, tarn->deal.lvFundEndDates[i],
                     tarn->market.szYieldCurve);
        double dFund_cvgMargin =
            1.0 - coverage(tarn->deal.lvFundStartDates[i],
                           tarn->deal.lvFundEndDates[i], tarn->deal.basisFund) *
                      (tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i]);
        aux->dHistFundPV += dFund_Start_DF - dFund_cvgMargin * dFund_End_DF;
      }
    }
    return 0;
  }

  /* if the swap has not knocked out  , calculate any residual coupon cash flows
   */
  if ((aux->i1stCpn > 0) &&
      (tarn->deal.lvCpnPayDates[aux->i1stCpn - 1] >= lFirstPayment))
    aux->dHistCpnPV = (tarn->deal.dvCoupon[aux->i1stCpn - 1] +
                       tarn->deal.dvGearing[aux->i1stCpn - 1] *
                           tarn->deal.dvFltrHistFixings[aux->i1stCpn - 1]) *
                      aux->dvCouponCvg[aux->i1stCpn - 1] *
                      aux->dvCouponPayDF[aux->i1stCpn - 1];

  /* if the swap has not knocked out  , calculate any residual funding PV up to
   * the first fixing date */
  iTemp = aux->i1stCpn == 0 ? 0 : aux->i1stCpn - 1;
  for (i = aux->ivFundStartIndex[iTemp]; i <= aux->ivFundEndIndex[aux->i1stCpn];
       i++) {
    if (tarn->deal.lvFundFixingDates[i] < lFirstFixing &&
        tarn->deal.lvFundEndDates[i] >= lFirstPayment)
      aux->dHistFundPV +=
          (tarn->deal.dvFundHistFixings[i] + tarn->deal.dvFundMargin[i]) *
          aux->dvFundCvg[i] * aux->dvFundDomPayDF[i];

    if (tarn->deal.lvFundFixingDates[i] >= lFirstFixing &&
        tarn->deal.lvFundEndDates[i] >= lFirstPayment) {
      double dFund_Start_DF =
          swp_f_df(tarn->market.lToday, tarn->deal.lvFundStartDates[i],
                   tarn->market.szYieldCurve);
      double dFund_End_DF =
          swp_f_df(tarn->market.lToday, tarn->deal.lvFundEndDates[i],
                   tarn->market.szYieldCurve);
      double dFund_cvgMargin =
          1.0 - coverage(tarn->deal.lvFundStartDates[i],
                         tarn->deal.lvFundEndDates[i], tarn->deal.basisFund) *
                    (tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i]);
      aux->dHistFundPV += dFund_Start_DF - dFund_cvgMargin * dFund_End_DF;
    }
  }

  /* return */
  return 0;
}

char *TARN_calcHistory_FINAL(TARN_Struct *tarn, TARN_AUX *aux) {
  /* Variable declaration */
  long lFirstFixing = tarn->market.lToday + !tarn->market.eodFixFlag;
  long lFirstPayment = tarn->market.lToday + !tarn->market.eodPayFlag;
  char *err = 0;
  int i = 0;
  double dCoupon = 0.0;
  int iTemp;

  /* Find the value of the cumulative coupon and see if it has knocked out.  If
     it hasn't knocked out  , find the first event date */
  tarn->output.dAccretedCoupon = 0.0;
  while ((aux->i1stCpn < tarn->deal.nCouponDates) &&
         (tarn->deal.lvFltrFixingDates[aux->i1stCpn] < lFirstFixing)) {
    /* Calculate the coupon */
    dCoupon = tarn->deal.dvCoupon[aux->i1stCpn] +
              tarn->deal.dvGearing[aux->i1stCpn] *
                  tarn->deal.dvFltrHistFixings[aux->i1stCpn];
    if (tarn->deal.ivIsFloored[aux->i1stCpn])
      dCoupon = dCoupon < tarn->deal.dvFloor[aux->i1stCpn]
                    ? tarn->deal.dvFloor[aux->i1stCpn]
                    : dCoupon;
    if (tarn->deal.ivIsCapped[aux->i1stCpn])
      dCoupon = dCoupon > tarn->deal.dvCap[aux->i1stCpn]
                    ? tarn->deal.dvCap[aux->i1stCpn]
                    : dCoupon;

    /* Check if we have previously knocked out */
    if (tarn->output.iIsKnockedOut) {
      tarn->output.dAccretedCoupon *=
          1.0 + aux->dvCouponCvg[aux->i1stCpn] *
                    tarn->deal.dvFltrHistFixings[aux->i1stCpn - 1];
    }
    /* or if we knocked out this time */
    else if (tarn->output.dCumulCoupon + dCoupon * aux->dvCumulCvg[i] >=
             tarn->deal.dTarget) {
      tarn->output.iIsKnockedOut = aux->i1stCpn + 1;
      tarn->output.dAccretedCoupon = 1.0 + tarn->deal.dTarget;
      tarn->output.dCumulCoupon = tarn->deal.dTarget;
    }
    /* or add in the coupon */
    else {
      tarn->output.dCumulCoupon += dCoupon * aux->dvCumulCvg[i];
    }
    ++aux->i1stCpn;
  }

  /* if there is only 1 coupon remaining  , only then is there a cash flow */
  if (aux->i1stCpn == (tarn->deal.nCouponDates - 1)) {
    if (tarn->output.iIsKnockedOut) {
      tarn->output.dAccretedCoupon *=
          1.0 + aux->dvCouponCvg[aux->i1stCpn] *
                    tarn->deal.dvFltrHistFixings[aux->i1stCpn - 1];
    } else {
      tarn->output.iIsKnockedOut = tarn->deal.nCouponDates;
      tarn->output.dAccretedCoupon = 1.0 + tarn->deal.dTarget;
    }
    aux->dHistCpnPV =
        (tarn->output.dAccretedCoupon - 1.0) * aux->dvCouponPayDF[aux->i1stCpn];
  }
  /* or acrete to the next (i.e.  , first) coupon date */
  else {
    tarn->output.dAccretedCoupon *=
        1.0 + aux->dvCouponCvg[aux->i1stCpn] *
                  tarn->deal.dvFltrHistFixings[aux->i1stCpn - 1];
  }

  /* Calculate any residual funding PV up to the first fixing date */
  iTemp = aux->i1stCpn == 0 ? 0 : aux->i1stCpn - 1;
  for (i = aux->ivFundStartIndex[iTemp]; i <= aux->ivFundEndIndex[aux->i1stCpn];
       i++) {
    if (tarn->deal.lvFundFixingDates[i] < lFirstFixing &&
        tarn->deal.lvFundEndDates[i] >= lFirstPayment)
      aux->dHistFundPV +=
          (tarn->deal.dvFundHistFixings[i] + tarn->deal.dvFundMargin[i]) *
          aux->dvFundCvg[i] * aux->dvFundDomPayDF[i];

    if (tarn->deal.lvFundFixingDates[i] >= lFirstFixing &&
        tarn->deal.lvFundEndDates[i] >= lFirstPayment) {
      double dFund_Start_DF =
          swp_f_df(tarn->market.lToday, tarn->deal.lvFundStartDates[i],
                   tarn->market.szYieldCurve);
      double dFund_End_DF =
          swp_f_df(tarn->market.lToday, tarn->deal.lvFundEndDates[i],
                   tarn->market.szYieldCurve);
      double dFund_cvgMargin =
          1.0 - coverage(tarn->deal.lvFundStartDates[i],
                         tarn->deal.lvFundEndDates[i], tarn->deal.basisFund) *
                    (tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i]);
      aux->dHistFundPV += dFund_Start_DF - dFund_cvgMargin * dFund_End_DF;
    }
  }

  /* return */
  return 0;
}

char *TARN_calcHistory_SWAPKO(TARN_Struct *tarn, TARN_AUX *aux) {
  /* Variable declaration */
  long lFirstFixing = tarn->market.lToday + !tarn->market.eodFixFlag;
  long lFirstPayment = tarn->market.lToday + !tarn->market.eodPayFlag;
  char *err = 0;
  int i = 0;
  double dCoupon = 0.0;
  int iTemp;

  /* Find the value of the cumulative coupon and see if it has knocked out.  If
     it hasn't knocked out  , find the first event date */
  while ((!tarn->output.iIsKnockedOut) &&
         (aux->i1stCpn < tarn->deal.nCouponDates) &&
         (tarn->deal.lvFltrFixingDates[aux->i1stCpn] < lFirstFixing)) {
    dCoupon = tarn->deal.dvCoupon[aux->i1stCpn] +
              tarn->deal.dvGearing[aux->i1stCpn] *
                  tarn->deal.dvFltrHistFixings[aux->i1stCpn];
    if (tarn->deal.ivIsFloored[aux->i1stCpn])
      dCoupon = dCoupon < tarn->deal.dvFloor[aux->i1stCpn]
                    ? tarn->deal.dvFloor[aux->i1stCpn]
                    : dCoupon;
    if (tarn->deal.ivIsCapped[aux->i1stCpn])
      dCoupon = dCoupon > tarn->deal.dvCap[aux->i1stCpn]
                    ? tarn->deal.dvCap[aux->i1stCpn]
                    : dCoupon;
    if (tarn->output.dCumulCoupon + dCoupon * aux->dvCumulCvg[i] >=
        tarn->deal.dTarget) {
      tarn->output.iIsKnockedOut = aux->i1stCpn + 1;
      tarn->output.dCumulCoupon += dCoupon * aux->dvCumulCvg[i];
    } else {
      tarn->output.dCumulCoupon += dCoupon * aux->dvCumulCvg[i];
      ++aux->i1stCpn;
    }
  }

  /* If the swap has knocked out  , see if there are any remaining cash flows
   * and calculate their PV */
  if (tarn->output.iIsKnockedOut) {
    /* Calculate the coupon PV */
    if (tarn->deal.lvCpnPayDates[aux->i1stCpn] >= lFirstPayment)
      aux->dHistCpnPV = dCoupon * aux->dvCouponCvg[aux->i1stCpn] *
                        aux->dvCouponPayDF[aux->i1stCpn];

    /* Calculate the funding PV */
    for (i = aux->ivFundStartIndex[aux->i1stCpn];
         i <= aux->ivFundEndIndex[aux->i1stCpn]; i++) {
      if (tarn->deal.lvFundFixingDates[i] < lFirstFixing &&
          tarn->deal.lvFundEndDates[i] >= lFirstPayment)
        aux->dHistFundPV +=
            (tarn->deal.dvFundHistFixings[i] + tarn->deal.dvFundMargin[i]) *
            aux->dvFundCvg[i] * aux->dvFundDomPayDF[i];

      if (tarn->deal.lvFundFixingDates[i] >= lFirstFixing &&
          tarn->deal.lvFundEndDates[i] >= lFirstPayment) {
        double dFund_Start_DF =
            swp_f_df(tarn->market.lToday, tarn->deal.lvFundStartDates[i],
                     tarn->market.szYieldCurve);
        double dFund_End_DF =
            swp_f_df(tarn->market.lToday, tarn->deal.lvFundEndDates[i],
                     tarn->market.szYieldCurve);
        double dFund_cvgMargin =
            1.0 - coverage(tarn->deal.lvFundStartDates[i],
                           tarn->deal.lvFundEndDates[i], tarn->deal.basisFund) *
                      (tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i]);
        aux->dHistFundPV += dFund_Start_DF - dFund_cvgMargin * dFund_End_DF;
      }
    }
    return 0;
  }

  /* if the swap has not knocked out  , calculate any residual coupon cash flows
   */
  if ((aux->i1stCpn > 0) &&
      (tarn->deal.lvCpnPayDates[aux->i1stCpn - 1] >= lFirstPayment))
    aux->dHistCpnPV = (tarn->deal.dvCoupon[aux->i1stCpn - 1] +
                       tarn->deal.dvGearing[aux->i1stCpn - 1] *
                           tarn->deal.dvFltrHistFixings[aux->i1stCpn - 1]) *
                      aux->dvCouponCvg[aux->i1stCpn - 1] *
                      aux->dvCouponPayDF[aux->i1stCpn - 1];

  /* if the swap has not knocked out  , calculate any residual funding PV up to
   * the first fixing date */
  iTemp = aux->i1stCpn == 0 ? 0 : aux->i1stCpn - 1;
  for (i = aux->ivFundStartIndex[iTemp]; i <= aux->ivFundEndIndex[aux->i1stCpn];
       i++) {
    if (tarn->deal.lvFundFixingDates[i] < lFirstFixing &&
        tarn->deal.lvFundEndDates[i] >= lFirstPayment)
      aux->dHistFundPV +=
          (tarn->deal.dvFundHistFixings[i] + tarn->deal.dvFundMargin[i]) *
          aux->dvFundCvg[i] * aux->dvFundDomPayDF[i];

    if (tarn->deal.lvFundFixingDates[i] >= lFirstFixing &&
        tarn->deal.lvFundEndDates[i] >= lFirstPayment) {
      double dFund_Start_DF =
          swp_f_df(tarn->market.lToday, tarn->deal.lvFundStartDates[i],
                   tarn->market.szYieldCurve);
      double dFund_End_DF =
          swp_f_df(tarn->market.lToday, tarn->deal.lvFundEndDates[i],
                   tarn->market.szYieldCurve);
      double dFund_cvgMargin =
          1.0 - coverage(tarn->deal.lvFundStartDates[i],
                         tarn->deal.lvFundEndDates[i], tarn->deal.basisFund) *
                    (tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i]);
      aux->dHistFundPV += dFund_Start_DF - dFund_cvgMargin * dFund_End_DF;
    }
  }

  /* return */
  return 0;
}

char *TARN_calcHistory(TARN_Struct *tarn, TARN_AUX *aux) {
  /* Variable decalarations */
  char *err = 0;

  /* initialize the return variables */
  tarn->output.iIsKnockedOut = 0;
  tarn->output.iIsExpired = 0;
  tarn->output.dCumulCoupon = 0.0;

  /* Calculate the history */
  switch (tarn->deal.couponType) {
  case TN_ZC:
    err = TARN_calcHistory_ZC(tarn, aux);
    break;
  case TN_SWAP:
    err = TARN_calcHistory_SWAP(tarn, aux);
    break;
  case TN_SWAPKO:
    err = TARN_calcHistory_SWAPKO(tarn, aux);
    break;
  case TN_FINAL:
    err = TARN_calcHistory_FINAL(tarn, aux);
    break;
  default:
    err = "Unknown error";
  }
  if (err)
    return err;

  /* if there are no more fixings left  , then mark it as finished */
  if (aux->i1stCpn == tarn->deal.nCouponDates)
    tarn->output.iIsExpired = 1;
  else
    tarn->output.iIsExpired = 0;

  /* end */
  return err;
}
