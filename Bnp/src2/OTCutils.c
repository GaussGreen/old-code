/* ==========================================================
        FILENAME :			OTCutils.c

        PURPOSE:			Calculates One Time Callables using a
   smile on the Fx by a copula based method

        AUTHOR:				J. Dinh

        DATE:				25-FEB-2004
   ========================================================== */

#include "Fx3FBetaDLMCalculations.h"
#include "Fx3FBetaDLMCalibration.h"
#include "Fx3FBetaDLMUtil.h"

#include "CopulaGaussian.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "OTCutils.h"
#include "SABRCalibRRBT.h"
#include "math.h"
#include "num_h_interp.h"
#include "opfnctns.h"
#include "opsabrgeneric.h"
#include "opsabrgenericinterp.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srtaccess.h"
#include "swp_h_vol.h"

#include "nag.h"
#include "nag_stdlib.h"
#include "nage04.h"
#include "nagx02.h"

/* =========================================================================
 Macro
========================================================================= */
#define POS_VAL(X) ((X) > 0 ? (X) : 0)

#define CALL_VAL(FWD, STRIKE, D, S) ((FWD)*norm((D) + (S)) - (STRIKE)*norm((D)))

#define PUT_VAL(FWD, STRIKE, D, S)                                             \
  (-(FWD)*norm(-(D) - (S)) + (STRIKE)*norm(-(D)))

#define OPT_VAL_MACRO(TYPE, FWD, STRIKE, STD, HALF_STD)                        \
  ((TYPE) == 0                                                                 \
       ? 0.0                                                                   \
       : ((TYPE) == 1                                                          \
              ? POS_VAL((FWD) - (STRIKE))                                      \
              : ((TYPE) == 2                                                   \
                     ? POS_VAL((STRIKE) - (FWD))                               \
                     : ((TYPE) == 3 ? CALL_VAL((FWD), (STRIKE),                \
                                               log((FWD) / (STRIKE)) / (STD) - \
                                                   (HALF_STD),                 \
                                               (STD))                          \
                                    : PUT_VAL((FWD), (STRIKE),                 \
                                              log((FWD) / (STRIKE)) / (STD) -  \
                                                  (HALF_STD),                  \
                                              (STD))))))

/* =========================================================================
OTC STRUCTURE
==========================================================================*/
void cpd_init_otc_params(
    otcpd_params *otc_params, CPD_STR cpd, int ncall, int nSimul, int Do_pecs,
    int nPoints, int nStd, double CummulPrecision, int CummulLinear, double std,
    int nIter, int PayoffFunction, double fwdSmileVisuNstd,
    /* otc */
    int fwdVolMethod,        /*	0=Sliding 3F vol; 1=Converging 3F vol; 2=Sliding
                                Cvg Sbeta; 3=Cvg Cvg Sbeta */
    int smileOtc,            /*	use smile parameters for the OTC */
    int smileFee,            /*	use smile parameters for the OTC Fees */
    int smileModel, int otc, /*	which call to keep */
    double BMpi,
    /* Correl */
    int firstIndex, int secondIndex, int firstLong, int secondLong,
    long correlTstar, char *CPDsigma, char *CPDalpha, char *CPDbeta,
    char *CPDrho,
    /* Change the Funding for speed */
    int FundingSpeedUp,
    /* Do not calculate all */
    int keepFirstN, int oneOutOfN,
    /*Fast MC*/
    int FMC_do, double FMC_precision, int FMC_min_paths,
    /*Which GMA*/
    int use_GMA,
    /* Change the calibration strategy */
    int LongShort) {
  int i;

  otc_params->Seed = -123456789;

  otc_params->Options_number = cpd->num_calls;

  otc_params->BMnStd = nStd;
  otc_params->BMalpha = 0.0;
  otc_params->BMbeta = 1.0;
  otc_params->BMlambda = 0.0;
  otc_params->BMpi = BMpi;
  otc_params->BMsigma = 0.1;
  otc_params->BMfwd1 = 100.0;
  otc_params->BMfwd2 = 100.0;
  otc_params->BMsig1 = 0.1;
  otc_params->BMsig2 = 0.1;

  otc_params->COPULAgauss = NULL;
  otc_params->COPULAdo_pecs = Do_pecs;
  otc_params->COPULAmatrix = NULL;
  otc_params->COPULAnSimul = nSimul;
  otc_params->CUMULnPoints = nPoints;
  otc_params->CUMULnStd = nStd;
  otc_params->CUMUL = NULL;
  otc_params->CUMULPrecision = CummulPrecision;
  otc_params->CUMULlinear = CummulLinear;

  otc_params->FwdVolMethod = fwdVolMethod;
  otc_params->FSMsigma = CPDsigma;
  otc_params->FSMalpha = CPDalpha;
  otc_params->FSMbeta = CPDbeta;
  otc_params->FSMrho = CPDrho;

  otc_params->OTC = otc;
  otc_params->OTCsmile = smileOtc;
  otc_params->OTCsmileFee = smileFee;
  otc_params->OTCfwdSmileVisuStd = fwdSmileVisuNstd;
  otc_params->OTCpayoff = PayoffFunction;

  if (otc_params->OTCpayoff == 6)
    otc_params->OTCsmileAndFlat = 1;
  else
    otc_params->OTCsmileAndFlat = 0;

  otc_params->smile_sigma = 0.1;
  otc_params->smile_sigmaBeta = 0.1;
  otc_params->smile_alpha = 0.0;
  otc_params->smile_beta = 1.0;
  otc_params->smile_rho = 1.0;

  otc_params->smileModel = smileModel;

  otc_params->SSLnIter = nIter;
  otc_params->SSLstd = std;
  otc_params->SSLshift = 0.0;
  otc_params->SSLvolDown = 0.1;
  otc_params->SSLvolUp = 0.1;

  otc_params->CORRELfirstIndex = firstIndex;
  otc_params->CORRELsecondIndex = secondIndex;
  otc_params->CORRELfirstLong = firstLong;
  otc_params->CORRELsecondLong = secondLong;
  otc_params->CORRELtStar = correlTstar;

  if (cpd->type == -1)
    otc_params->KO_do_optim = 1;
  else
    otc_params->KO_do_optim = 0;

  otc_params->KO_FxGauss = NULL;
  otc_params->KO_FxMatrix = NULL;
  otc_params->KO_RatesGauss = NULL;
  otc_params->KO_RatesMatrix = NULL;
  otc_params->KO_FMC = FMC_do;
  otc_params->KO_FMC_precision = FMC_precision;
  otc_params->KO_FMC_min_paths = FMC_min_paths;

  otc_params->SpeedFunding = FundingSpeedUp;

  otc_params->use_GMA = use_GMA;
  otc_params->LongShort = LongShort;

  otc_params->OTCcalculate = NULL;
  otc_params->OTCcalculate =
      (int *)calloc(otc_params->Options_number, sizeof(int));

  if (oneOutOfN == 1) {
    for (i = 0; i < otc_params->Options_number; i++)
      otc_params->OTCcalculate[i] = 1;
  } else {
    otc_params->OTCcalculate[otc_params->Options_number - 1] = 1;

    for (i = 0; i < otc_params->Options_number - 1; i++) {
      if (i < keepFirstN ||
          (i >= keepFirstN && (double)(i - keepFirstN + 1) / oneOutOfN ==
                                  (i - keepFirstN + 1) / oneOutOfN))
        otc_params->OTCcalculate[i] = 1;
    }
  }
}

Err cpd_free_otc_params(otcpd_params *otc_params) {
  Err err = NULL;

  if (otc_params->COPULAgauss)
    free_dmatrix(otc_params->COPULAgauss, 0, otc_params->COPULAnSimul - 1, 0,
                 3 - 1);

  if (otc_params->CUMUL)
    free_dmatrix(otc_params->CUMUL, 0, 2 - 1, 0, otc_params->CUMULnPoints - 1);

  if (otc_params->COPULAmatrix)
    free_dmatrix(otc_params->COPULAmatrix, 0, otc_params->COPULAnSimul - 1, 0,
                 3 - 1 + otc_params->OTCsmileAndFlat);

  if (otc_params->KO_FxGauss)
    free_dmatrix(otc_params->KO_FxGauss, 0, otc_params->COPULAnSimul - 1, 0,
                 otc_params->Options_number - 1);

  if (otc_params->KO_FxMatrix)
    free_dmatrix(otc_params->KO_FxMatrix, 0, otc_params->COPULAnSimul - 1, 0,
                 otc_params->Options_number - 1);

  if (otc_params->KO_RatesGauss)
    free_dmatrix(otc_params->KO_RatesGauss, 0, otc_params->COPULAnSimul - 1, 0,
                 2 - 1);

  if (otc_params->KO_RatesMatrix)
    free_dmatrix(otc_params->KO_RatesMatrix, 0, otc_params->COPULAnSimul - 1, 0,
                 2 - 1);

  if (otc_params->OTCcalculate)
    free(otc_params->OTCcalculate);

  otc_params->COPULAgauss = NULL;
  otc_params->CUMUL = NULL;
  otc_params->COPULAmatrix = NULL;
  otc_params->KO_FxGauss = NULL;
  otc_params->KO_FxMatrix = NULL;
  otc_params->KO_RatesGauss = NULL;
  otc_params->KO_RatesMatrix = NULL;
  otc_params->OTCcalculate = NULL;

  return err;
}

void cpd_init_otc_precalc(otcpd_precalc *otc_precalc) {
  otc_precalc->df_const_dom_fx_val = NULL;
  otc_precalc->df_lin_dom_fx_val = NULL;
  otc_precalc->df_const_for_fx_val = NULL;
  otc_precalc->df_lin_for_fx_val = NULL;

  otc_precalc->df_const_pd_pay = NULL;
  otc_precalc->df_lin_pd_pay = NULL;

  otc_precalc->df_const_fund_pay = NULL;
  otc_precalc->df_lin_fund_pay = NULL;

  otc_precalc->df_const_fund_start = NULL;
  otc_precalc->df_lin_fund_start = NULL;

  otc_precalc->pd_fwd3F_vol = NULL;
  otc_precalc->pd_fwdsmile_vol = NULL;
  otc_precalc->pd_fwdsmile_alpha = NULL;
  otc_precalc->pd_fwdsmile_beta = NULL;
  otc_precalc->pd_fwdsmile_rho = NULL;
  otc_precalc->pd_fwdsmile_pi = NULL;

  otc_precalc->FxAdj = NULL;
  otc_precalc->FeeFxAdj = NULL;

  otc_precalc->NotFxAdj = NULL;
  otc_precalc->NotFeeFxAdj = NULL;

  otc_precalc->KO_barrier = NULL;
}

void cpd_free_otc_precalc(otcpd_precalc *otc_precalc, int nCall,
                          int nCoupons_pd, int nCoupons_fund) {
  if (otc_precalc->df_const_dom_fx_val)
    free_dmatrix(otc_precalc->df_const_dom_fx_val, 0, nCall - 1, 0,
                 nCoupons_pd + 2);
  if (otc_precalc->df_lin_dom_fx_val)
    free_dmatrix(otc_precalc->df_lin_dom_fx_val, 0, nCall - 1, 0,
                 nCoupons_pd + 2);
  if (otc_precalc->df_const_for_fx_val)
    free_dmatrix(otc_precalc->df_const_for_fx_val, 0, nCall - 1, 0,
                 nCoupons_pd + 2);
  if (otc_precalc->df_lin_for_fx_val)
    free_dmatrix(otc_precalc->df_lin_for_fx_val, 0, nCall - 1, 0,
                 nCoupons_pd + 2);
  if (otc_precalc->df_const_pd_pay)
    free_dmatrix(otc_precalc->df_const_pd_pay, 0, nCall - 1, 0,
                 nCoupons_pd + 2);
  if (otc_precalc->df_lin_pd_pay)
    free_dmatrix(otc_precalc->df_lin_pd_pay, 0, nCall - 1, 0, nCoupons_pd + 2);
  if (otc_precalc->df_const_fund_pay)
    free_dmatrix(otc_precalc->df_const_fund_pay, 0, nCall - 1, 0,
                 nCoupons_fund + 2);
  if (otc_precalc->df_lin_fund_pay)
    free_dmatrix(otc_precalc->df_lin_fund_pay, 0, nCall - 1, 0,
                 nCoupons_fund + 2);
  if (otc_precalc->df_const_fund_start)
    free(otc_precalc->df_const_fund_start);
  if (otc_precalc->df_lin_fund_start)
    free(otc_precalc->df_lin_fund_start);

  if (otc_precalc->pd_fwd3F_vol)
    free_dmatrix(otc_precalc->pd_fwd3F_vol, 0, nCall, 0, nCoupons_pd + 2);
  if (otc_precalc->pd_fwdsmile_alpha)
    free_dmatrix(otc_precalc->pd_fwdsmile_alpha, 0, nCall, 0, nCoupons_pd + 2);
  if (otc_precalc->pd_fwdsmile_beta)
    free_dmatrix(otc_precalc->pd_fwdsmile_beta, 0, nCall, 0, nCoupons_pd + 2);
  if (otc_precalc->pd_fwdsmile_rho)
    free_dmatrix(otc_precalc->pd_fwdsmile_rho, 0, nCall, 0, nCoupons_pd + 2);
  if (otc_precalc->pd_fwdsmile_pi)
    free_dmatrix(otc_precalc->pd_fwdsmile_pi, 0, nCall, 0, nCoupons_pd + 2);
  if (otc_precalc->pd_fwdsmile_vol)
    free_dmatrix(otc_precalc->pd_fwdsmile_vol, 0, nCall, 0, nCoupons_pd + 2);

  if (otc_precalc->FxAdj)
    free_dmatrix(otc_precalc->FxAdj, 0, nCall - 1, 0, nCoupons_pd - 1);
  if (otc_precalc->FeeFxAdj)
    free(otc_precalc->FeeFxAdj);

  if (otc_precalc->NotFxAdj)
    free_dmatrix(otc_precalc->NotFxAdj, 0, nCall - 1, 0, 2);
  if (otc_precalc->NotFeeFxAdj)
    free(otc_precalc->NotFeeFxAdj);

  if (otc_precalc->KO_barrier)
    free_dmatrix(otc_precalc->KO_barrier, 0, nCall - 1, 0, nCoupons_pd - 1);
}

Err cpd_alloc_otc_precalc(otcpd_precalc *otc_precalc, int nCall,
                          int nCoupons_pd, int nCoupons_fund, int nMktDates) {
  Err err = NULL;
  otc_precalc->df_const_dom_fx_val = dmatrix(0, nCall - 1, 0, nCoupons_pd + 2);
  otc_precalc->df_lin_dom_fx_val = dmatrix(0, nCall - 1, 0, nCoupons_pd + 2);
  otc_precalc->df_const_for_fx_val = dmatrix(0, nCall - 1, 0, nCoupons_pd + 2);
  otc_precalc->df_lin_for_fx_val = dmatrix(0, nCall - 1, 0, nCoupons_pd + 2);
  otc_precalc->df_const_pd_pay = dmatrix(0, nCall - 1, 0, nCoupons_pd + 2);
  otc_precalc->df_lin_pd_pay = dmatrix(0, nCall - 1, 0, nCoupons_pd + 2);
  otc_precalc->df_const_fund_pay = dmatrix(0, nCall - 1, 0, nCoupons_fund + 2);
  otc_precalc->df_lin_fund_pay = dmatrix(0, nCall - 1, 0, nCoupons_fund + 2);
  otc_precalc->df_const_fund_start = calloc(nCall, sizeof(double));
  otc_precalc->df_lin_fund_start = calloc(nCall, sizeof(double));

  otc_precalc->pd_fwd3F_vol = dmatrix(0, nCall, 0, nCoupons_pd + 2);
  otc_precalc->pd_fwdsmile_vol = dmatrix(0, nCall, 0, nCoupons_pd + 2);
  otc_precalc->pd_fwdsmile_alpha = dmatrix(0, nCall, 0, nCoupons_pd + 2);
  otc_precalc->pd_fwdsmile_beta = dmatrix(0, nCall, 0, nCoupons_pd + 2);
  otc_precalc->pd_fwdsmile_rho = dmatrix(0, nCall, 0, nCoupons_pd + 2);
  otc_precalc->pd_fwdsmile_pi = dmatrix(0, nCall, 0, nCoupons_pd + 2);

  otc_precalc->FxAdj = dmatrix(0, nCall - 1, 0, nCoupons_pd - 1);
  otc_precalc->FeeFxAdj = calloc(nCoupons_pd, sizeof(double));

  otc_precalc->NotFxAdj = dmatrix(0, nCall - 1, 0, 2);
  otc_precalc->NotFeeFxAdj =
      calloc(nCall + 1, sizeof(double)); /* +1 for the last notional */

  otc_precalc->KO_barrier = dmatrix(0, nCall - 1, 0, nCoupons_pd - 1);

  if (!otc_precalc->pd_fwd3F_vol || !otc_precalc->pd_fwdsmile_vol ||
      !otc_precalc->pd_fwdsmile_alpha || !otc_precalc->pd_fwdsmile_beta ||
      !otc_precalc->pd_fwdsmile_rho || !otc_precalc->pd_fwdsmile_pi ||
      !otc_precalc->FxAdj || !otc_precalc->FeeFxAdj ||
      !otc_precalc->KO_barrier) {
    err = "cpd_alloc_otc_precalc: memory allocation failure";
  }

  return err;
}

/* =========================================================================
 Precalculations
========================================================================= */
Err cpd_otc_precalc(otcpd_precalc *otc_precalc, otcpd_params *OTC_params,
                    CPD_STR cpd, CPD_UND und, SMILE_VOL_MARKET smile_mkt,
                    double *mergeTimes, int mergeNtimes, double *sigDom,
                    double *sigFor, double *sigFx, double *corDF,
                    double *corDFx, double *corFFx) {
  int i, j, k;
  double phiDom, phiFor, beta, fund_lda, fund_phi;
  double *forwardVec = NULL;
  Err err = NULL;

  err = cpd_alloc_otc_precalc(otc_precalc, cpd->num_calls, cpd->pd_leg->num_cpn,
                              cpd->fund_leg->num_cpn, smile_mkt->num_vols);
  if (err)
    goto FREE_RETURN;

  /*=====================================
  Is the funding optimizeable?
          only false in the case of a $ funded
          CPD with fixing in advance
  =====================================*/

  if (cpd->fund_leg->dom_for && OTC_params->SpeedFunding) {
    for (k = 0; k < cpd->num_calls; k++) {
      if (cpd->fund_leg->cpn[cpd->call[k].fund_idx].pay_date >
          cpd->pd_leg->cpn[cpd->call[k].pd_idx].fx_val_date)
        OTC_params->SpeedFunding = 0;
    }
  }

  /* ===================================
  discount factors precalculations
  =================================== */
  if (cpd->fund_leg->dom_for)
    fund_lda = und->lda_for;
  else
    fund_lda = und->lda_dom;

  for (i = 0; i < cpd->num_calls; i++) {
    err = OTCdfPhi(und->lda_dom, cpd->call[i].ex_time, und->sigma_time_rates,
                   und->sigma_dom, und->sigma_n_rates, &phiDom);

    if (err)
      goto FREE_RETURN;

    err = OTCdfPhi(und->lda_for, cpd->call[i].ex_time, und->sigma_time_rates,
                   und->sigma_for, und->sigma_n_rates, &phiFor);

    if (err)
      goto FREE_RETURN;

    if (cpd->fund_leg->dom_for)
      fund_phi = phiFor;
    else
      fund_phi = phiDom;

    beta = (1.0 - exp(-fund_lda *
                      (cpd->fund_leg->cpn[cpd->call[i].fund_idx].start_time -
                       cpd->call[i].ex_time))) /
           fund_lda;
    otc_precalc->df_lin_fund_start[i] = -beta;
    otc_precalc->df_const_fund_start[i] = -0.5 * beta * beta * fund_phi;

    for (j = cpd->call[i].pd_idx; j < cpd->pd_leg->num_cpn; j++) {
      beta = (1.0 - exp(-und->lda_dom * (cpd->pd_leg->cpn[j].fx_val_time -
                                         cpd->call[i].ex_time))) /
             und->lda_dom;
      otc_precalc->df_lin_dom_fx_val[i][j] = -beta;
      otc_precalc->df_const_dom_fx_val[i][j] = -0.5 * beta * beta * phiDom;

      if (cpd->pd_leg->cpn[j].fx_val_time != cpd->pd_leg->cpn[j].pay_time) {
        beta = (1.0 - exp(-und->lda_dom * (cpd->pd_leg->cpn[j].pay_time -
                                           cpd->call[i].ex_time))) /
               und->lda_dom;
        otc_precalc->df_lin_pd_pay[i][j] = -beta;
        otc_precalc->df_const_pd_pay[i][j] = -0.5 * beta * beta * phiDom;
      } else {
        otc_precalc->df_lin_pd_pay[i][j] = otc_precalc->df_lin_dom_fx_val[i][j];
        otc_precalc->df_const_pd_pay[i][j] =
            otc_precalc->df_const_dom_fx_val[i][j];
      }

      beta = (1.0 - exp(-und->lda_for * (cpd->pd_leg->cpn[j].fx_val_time -
                                         cpd->call[i].ex_time))) /
             und->lda_for;
      otc_precalc->df_lin_for_fx_val[i][j] = -beta;
      otc_precalc->df_const_for_fx_val[i][j] = -0.5 * beta * beta * phiFor;
    }

    for (j = cpd->call[i].fund_idx; j < cpd->fund_leg->num_cpn; j++) {
      beta = (1.0 - exp(-fund_lda * (cpd->fund_leg->cpn[j].pay_time -
                                     cpd->call[i].ex_time))) /
             fund_lda;
      otc_precalc->df_lin_fund_pay[i][j] = -beta;
      otc_precalc->df_const_fund_pay[i][j] = -0.5 * beta * beta * fund_phi;
    }

    beta = (1.0 - exp(-und->lda_dom * (cpd->pd_leg->not_ref.fx_val_time -
                                       cpd->call[i].ex_time))) /
           und->lda_dom;
    otc_precalc->df_lin_dom_fx_val[i][cpd->pd_leg->num_cpn] = -beta;
    otc_precalc->df_const_dom_fx_val[i][cpd->pd_leg->num_cpn] =
        -0.5 * beta * beta * phiDom;

    if (cpd->pd_leg->not_ref.fx_val_time != cpd->pd_leg->not_ref.pay_time) {
      beta = (1.0 - exp(-und->lda_dom * (cpd->pd_leg->not_ref.pay_time -
                                         cpd->call[i].ex_time))) /
             und->lda_dom;
      otc_precalc->df_lin_pd_pay[i][cpd->pd_leg->num_cpn] = -beta;
      otc_precalc->df_const_pd_pay[i][cpd->pd_leg->num_cpn] =
          -0.5 * beta * beta * phiDom;
    } else {
      otc_precalc->df_lin_pd_pay[i][cpd->pd_leg->num_cpn] =
          otc_precalc->df_lin_dom_fx_val[i][cpd->pd_leg->num_cpn];
      otc_precalc->df_const_pd_pay[i][cpd->pd_leg->num_cpn] =
          otc_precalc->df_const_dom_fx_val[i][cpd->pd_leg->num_cpn];
    }

    beta = (1.0 - exp(-und->lda_for * (cpd->pd_leg->not_ref.fx_val_time -
                                       cpd->call[i].ex_time))) /
           und->lda_for;
    otc_precalc->df_lin_for_fx_val[i][cpd->pd_leg->num_cpn] = -beta;
    otc_precalc->df_const_for_fx_val[i][cpd->pd_leg->num_cpn] =
        -0.5 * beta * beta * phiFor;

    /* Should be Only if KO or usenotstr */
    for (j = 0; j < 2; j++) {
      if (i + j < cpd->num_calls) {
        beta = (1.0 - exp(-und->lda_dom * (cpd->call[i + j].fx_val_time -
                                           cpd->call[i].ex_time))) /
               und->lda_dom;
        otc_precalc->df_lin_dom_fx_val[i][cpd->pd_leg->num_cpn + j + 1] = -beta;
        otc_precalc->df_const_dom_fx_val[i][cpd->pd_leg->num_cpn + j + 1] =
            -0.5 * beta * beta * phiDom;

        if (cpd->call[i + j].fx_val_time != cpd->call[i + j].set_time) {
          beta = (1.0 - exp(-und->lda_dom * (cpd->call[i + j].set_time -
                                             cpd->call[i].ex_time))) /
                 und->lda_dom;
          otc_precalc->df_lin_pd_pay[i][cpd->pd_leg->num_cpn + j + 1] = -beta;
          otc_precalc->df_const_pd_pay[i][cpd->pd_leg->num_cpn + j + 1] =
              -0.5 * beta * beta * phiDom;
        } else {
          otc_precalc->df_lin_pd_pay[i][cpd->pd_leg->num_cpn + j + 1] =
              otc_precalc->df_lin_dom_fx_val[i][cpd->pd_leg->num_cpn + j + 1];
          otc_precalc->df_const_pd_pay[i][cpd->pd_leg->num_cpn + j + 1] =
              otc_precalc->df_const_dom_fx_val[i][cpd->pd_leg->num_cpn + j + 1];
        }

        beta = (1.0 - exp(-und->lda_for * (cpd->call[i + j].fx_val_time -
                                           cpd->call[i].ex_time))) /
               und->lda_for;
        otc_precalc->df_lin_for_fx_val[i][cpd->pd_leg->num_cpn + j + 1] = -beta;
        otc_precalc->df_const_for_fx_val[i][cpd->pd_leg->num_cpn + j + 1] =
            -0.5 * beta * beta * phiFor;

        beta = (1.0 - exp(-fund_lda *
                          (cpd->call[i + j].set_time - cpd->call[i].ex_time))) /
               fund_lda;
        otc_precalc->df_lin_fund_pay[i][cpd->fund_leg->num_cpn + j + 1] = -beta;
        otc_precalc->df_const_fund_pay[i][cpd->fund_leg->num_cpn + j + 1] =
            -0.5 * beta * beta * fund_phi;
      }
    }
  }

  /* ===================================
  Fx Adjustments for coupons
  =================================== */
  for (i = 0; i < cpd->num_calls; i++) {
    for (j = 0; j < cpd->pd_leg->num_cpn; j++) {
      if (cpd->pd_leg->cpn[j].fx_fix_time > cpd->call[i].set_time) {
        err = Fx3DtsFwdPayAdjustment_corr(
            cpd->call[i].ex_time, cpd->pd_leg->cpn[j].fx_val_time,
            cpd->pd_leg->cpn[j].fx_val_time, cpd->pd_leg->cpn[j].pay_time,
            cpd->pd_leg->cpn[j].fx_fix_time, und->sigma_time_rates,
            und->sigma_n_rates, und->sigma_dom, und->lda_dom, und->sigma_for,
            und->lda_for, und->sigma_time_fx, und->sigma_fx, und->sigma_n_fx,
            und->corr_times, und->correl_dom_for, und->correl_dom_fx,
            und->correl_for_fx, und->corr_n_times, &(otc_precalc->FxAdj[i][j]));

        if (err)
          goto FREE_RETURN;

        otc_precalc->FxAdj[i][j] = exp(otc_precalc->FxAdj[i][j]);
      } else {
        otc_precalc->FxAdj[i][j] = 1.0;
      }
    }
  }

  for (j = 0; j < cpd->pd_leg->num_cpn; j++) {
    err = Fx3DtsFwdPayAdjustment_corr(
        0.0, cpd->pd_leg->cpn[j].fx_val_time, cpd->pd_leg->cpn[j].fx_val_time,
        cpd->pd_leg->cpn[j].pay_time, cpd->pd_leg->cpn[j].fx_fix_time,
        und->sigma_time_rates, und->sigma_n_rates, und->sigma_dom, und->lda_dom,
        und->sigma_for, und->lda_for, und->sigma_time_fx, und->sigma_fx,
        und->sigma_n_fx, und->corr_times, und->correl_dom_for,
        und->correl_dom_fx, und->correl_for_fx, und->corr_n_times,
        &(otc_precalc->FeeFxAdj[j]));

    if (err)
      goto FREE_RETURN;

    otc_precalc->FeeFxAdj[j] = exp(otc_precalc->FeeFxAdj[j]);
  }

  /* ===================================
  Fx Adjustments for Not
  =================================== */
  for (i = 0; i < cpd->num_calls; i++) {
    for (j = 0; j < 2; j++) {
      if (i + j < cpd->num_calls && cpd->call->use_opt_str) {
        if (cpd->call[i + j].fx_fix_time > cpd->call[i].set_time) {
          err = Fx3DtsFwdPayAdjustment_corr(
              cpd->call[i].ex_time, cpd->call[i + j].fx_val_time,
              cpd->call[i + j].fx_val_time, cpd->call[i + j].set_time,
              cpd->call[i + j].fx_fix_time, und->sigma_time_rates,
              und->sigma_n_rates, und->sigma_dom, und->lda_dom, und->sigma_for,
              und->lda_for, und->sigma_time_fx, und->sigma_fx, und->sigma_n_fx,
              und->corr_times, und->correl_dom_for, und->correl_dom_fx,
              und->correl_for_fx, und->corr_n_times,
              &(otc_precalc->NotFxAdj[i][j]));

          if (err)
            goto FREE_RETURN;

          otc_precalc->NotFxAdj[i][j] = exp(otc_precalc->NotFxAdj[i][j]);
        } else {
          otc_precalc->NotFxAdj[i][j] = 1.0;
        }
      }
    }

    err = Fx3DtsFwdPayAdjustment_corr(
        cpd->call[i].ex_time, cpd->pd_leg->not_ref.fx_val_time,
        cpd->pd_leg->not_ref.fx_val_time, cpd->pd_leg->not_ref.pay_time,
        cpd->pd_leg->not_ref.fx_fix_time, und->sigma_time_rates,
        und->sigma_n_rates, und->sigma_dom, und->lda_dom, und->sigma_for,
        und->lda_for, und->sigma_time_fx, und->sigma_fx, und->sigma_n_fx,
        und->corr_times, und->correl_dom_for, und->correl_dom_fx,
        und->correl_for_fx, und->corr_n_times, &(otc_precalc->NotFxAdj[i][2]));

    if (err)
      goto FREE_RETURN;

    otc_precalc->NotFxAdj[i][j] = exp(otc_precalc->NotFxAdj[i][2]);
  }

  if (cpd->call->use_opt_str) {
    for (j = 0; j < cpd->num_calls; j++) {
      err = Fx3DtsFwdPayAdjustment_corr(
          0.0, cpd->call[j].fx_val_time, cpd->call[j].fx_val_time,
          cpd->call[j].set_time, cpd->call[j].fx_fix_time,
          und->sigma_time_rates, und->sigma_n_rates, und->sigma_dom,
          und->lda_dom, und->sigma_for, und->lda_for, und->sigma_time_fx,
          und->sigma_fx, und->sigma_n_fx, und->corr_times, und->correl_dom_for,
          und->correl_dom_fx, und->correl_for_fx, und->corr_n_times,
          &(otc_precalc->NotFeeFxAdj[j]));

      if (err)
        goto FREE_RETURN;

      otc_precalc->NotFeeFxAdj[j] = exp(otc_precalc->NotFeeFxAdj[j]);
    }
  } else {
    for (j = 0; j < cpd->num_calls; j++)
      otc_precalc->NotFeeFxAdj[j] = 1.0;
  }
  err = Fx3DtsFwdPayAdjustment_corr(
      0.0, cpd->pd_leg->not_ref.fx_val_time, cpd->pd_leg->not_ref.fx_val_time,
      cpd->pd_leg->not_ref.pay_time, cpd->pd_leg->not_ref.fx_fix_time,
      und->sigma_time_rates, und->sigma_n_rates, und->sigma_dom, und->lda_dom,
      und->sigma_for, und->lda_for, und->sigma_time_fx, und->sigma_fx,
      und->sigma_n_fx, und->corr_times, und->correl_dom_for, und->correl_dom_fx,
      und->correl_for_fx, und->corr_n_times,
      &(otc_precalc->NotFeeFxAdj[cpd->num_calls]));

  if (err)
    goto FREE_RETURN;

  otc_precalc->NotFeeFxAdj[cpd->num_calls] =
      exp(otc_precalc->NotFeeFxAdj[cpd->num_calls]);

  /* ===================================
  Barrier Adjustments
  =================================== */
  if (cpd->type) {
    for (i = 0; i < cpd->num_calls; i++) {
      for (j = 0; j < cpd->num_calls; j++) {
        if (j < i) {
          err = Fx3DtsFwdPayAdjustment_corr(
              0.0, cpd->call[j].ex_time2bd, cpd->call[j].ex_time2bd,
              cpd->call[i].ex_time2bd, cpd->call[j].ex_time,
              und->sigma_time_rates, und->sigma_n_rates, und->sigma_dom,
              und->lda_dom, und->sigma_for, und->lda_for, und->sigma_time_fx,
              und->sigma_fx, und->sigma_n_fx, und->corr_times,
              und->correl_dom_for, und->correl_dom_fx, und->correl_for_fx,
              und->corr_n_times, &(otc_precalc->KO_barrier[i][j]));

          if (err)
            goto FREE_RETURN;

          otc_precalc->KO_barrier[i][j] = exp(otc_precalc->KO_barrier[i][j]);
        } else {
          otc_precalc->KO_barrier[i][j] = 1.0;
        }
      }
    }
  }

  /* ===================================
   Forward SABR
  =================================== */
  forwardVec = calloc(cpd->pd_leg->num_cpn + 1 + 2, sizeof(double));

  if (!forwardVec) {
    err = "otc_caller: memory allocation failure";
    goto FREE_RETURN;
  }

  for (i = 0; i < cpd->pd_leg->num_cpn; i++)
    forwardVec[i] =
        und->spot_fx *
        swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->for_yc) /
        swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->dom_yc);
  ;

  forwardVec[cpd->pd_leg->num_cpn] =
      und->spot_fx *
      swp_f_df(und->today, cpd->pd_leg->not_ref.fx_val_date, und->for_yc) /
      swp_f_df(und->today, cpd->pd_leg->not_ref.fx_val_date, und->dom_yc);

  for (j = 0; j < cpd->num_calls; j++) {
    for (i = 0; i < 2; i++) {
      if (j + i < cpd->num_calls) {
        forwardVec[cpd->pd_leg->num_cpn + 1 + i] =
            und->spot_fx *
            swp_f_df(und->today, cpd->call[j + i].fx_val_date, und->for_yc) /
            swp_f_df(und->today, cpd->call[j + i].fx_val_date, und->dom_yc);
      }
    }

    if (OTC_params->OTCpayoff == 6) {
      /* =================================
      First Populate the 3F fwd vols and then the fwd smiles
      ================================= */

      err = OTCgetFwdSABR(
          forwardVec, cpd, und, smile_mkt, j, otc_precalc->pd_fwd3F_vol[j],
          otc_precalc->pd_fwdsmile_alpha[j], otc_precalc->pd_fwdsmile_beta[j],
          otc_precalc->pd_fwdsmile_rho[j], otc_precalc->pd_fwdsmile_pi[j],
          cpd->pd_leg->num_cpn + 1, mergeTimes, mergeNtimes, sigDom,
          und->lda_dom, sigFor, und->lda_for, sigFx, corDF, corDFx, corFFx,
          OTC_params->FSMsigma, OTC_params->FSMalpha, OTC_params->FSMbeta,
          OTC_params->FSMrho,
          4); /* 0=Sliding 3F vol; 1=Converging 3F vol; 2=Sliding Cvg Sbeta;
                 3=Cvg Cvg Sbeta; 4=3F*/

      if (err)
        goto FREE_RETURN;
    }

    err = OTCgetFwdSABR(
        forwardVec, cpd, und, smile_mkt, j, otc_precalc->pd_fwdsmile_vol[j],
        otc_precalc->pd_fwdsmile_alpha[j], otc_precalc->pd_fwdsmile_beta[j],
        otc_precalc->pd_fwdsmile_rho[j], otc_precalc->pd_fwdsmile_pi[j],
        cpd->pd_leg->num_cpn + 1, mergeTimes, mergeNtimes, sigDom, und->lda_dom,
        sigFor, und->lda_for, sigFx, corDF, corDFx, corFFx,
        OTC_params->FSMsigma, OTC_params->FSMalpha, OTC_params->FSMbeta,
        OTC_params->FSMrho,
        OTC_params
            ->FwdVolMethod); /* 0=Sliding 3F vol; 1=Converging 3F vol; 2=Sliding
                                Cvg Sbeta; 3=Cvg Cvg Sbeta; 4=3F*/

    if (err)
      goto FREE_RETURN;
  }

FREE_RETURN:
  if (forwardVec)
    free(forwardVec);
  return err;
}

/* =========================================================================
 Calibration of the sum of shifted lognormal parameters to SABR to obtain a
copula based estimation of OTC
========================================================================= */
Err CalibSSLtoSABR(double t,   /* Maturity in years */
                   double fwd, /* Forward */
                   double std, /* # std to calibrate skew */
                   /* SABR Parameters */
                   double sigma, double alpha, double beta, double rho,
                   /* Results */
                   double *VolUp, double *VolDown, double *Fwd_Shift,
                   double *error,
                   /* Parameters */
                   int nIter) {
  double Vol_SL, Beta_SL, k1, v1, k2, v2;

  double *dStrikes = NULL, *dWeights = NULL, *dParams = NULL,
         *dMarketprices = NULL;
  long *lUseparams = NULL;

  int nProducts, nParams;

  int volIter, iterLimit = 10;
  double oldPrice, price, oldVol, newVol, vol, tolerance, expUp, expDown;

  Err err = NULL;

  nProducts = 3;

  /* Memory allocation */
  dStrikes = dvector(1, nProducts);
  dMarketprices = dvector(1, nProducts);
  dWeights = dvector(1, nProducts);

  if (!dStrikes || !dMarketprices || !dWeights) {
    err = "CalibSSLtoSABR: Memory allocation error";
    goto FREE_RETURN;
  }

  /* Find the strikes and vols to prices to */
  k1 = fwd * exp(std * sigma * sqrt(t));
  k2 = fwd * exp(-std * sigma * sqrt(t));

  srt_f_optsarbvol(fwd, k1, t, sigma, alpha, beta, rho, SRT_LOGNORMAL,
                   SRT_LOGNORMAL, &v1);

  srt_f_optsarbvol(fwd, k2, t, sigma, alpha, beta, rho, SRT_LOGNORMAL,
                   SRT_LOGNORMAL, &v2);
  dStrikes[1] = k1;
  dStrikes[2] = fwd;
  dStrikes[3] = k2;

  dMarketprices[1] = srt_f_optblksch(fwd, k1, v1, t, 1.0, SRT_CALL, PREMIUM);
  dMarketprices[2] =
      srt_f_optblksch(fwd, fwd, sigma, t, 1.0, SRT_CALL, PREMIUM);
  dMarketprices[3] = srt_f_optblksch(fwd, k2, v2, t, 1.0, SRT_CALL, PREMIUM);

  /* Finds a first guess */
  err = find_shifted_log_params(t, fwd, k1, v1, k2, v2, Fwd_Shift, &Vol_SL,
                                &Beta_SL);

  if (err)
    goto FREE_RETURN;

  nParams = 6;
  lUseparams = lngvector(1, nParams);
  dParams = dvector(1, nParams);

  if (!lUseparams || !dParams) {
    err = "CalibSSLtoSABR: Memory allocation error";
    goto FREE_RETURN;
  }

  dParams[1] = fwd;
  dParams[2] = t;
  dParams[3] = Vol_SL;
  dParams[4] = *Fwd_Shift;
  dParams[5] = alpha / 2.0 * sqrt(t);
  dParams[6] = -alpha / 2.0 * sqrt(t);

  /* Adjust the ATM vol to match the market */
  volIter = 0;
  tolerance = fwd / 100000;
  expUp = exp(dParams[5]);
  expDown = exp(dParams[6]);
  oldVol = Vol_SL;
  newVol = Vol_SL * 0.95;
  err = SSLprice(fwd, t, fwd, oldVol * expUp, oldVol * expDown, dParams[4], 1,
                 &oldPrice);
  if (err)
    goto FREE_RETURN;
  err = SSLprice(fwd, t, fwd, newVol * expUp, newVol * expDown, dParams[4], 1,
                 &price);
  if (err)
    goto FREE_RETURN;

  do {
    vol = newVol -
          (price - dMarketprices[2]) / ((price - oldPrice) / (newVol - oldVol));
    oldVol = newVol;
    newVol = vol;
    oldPrice = price;
    err = SSLprice(fwd, t, fwd, vol * expUp, vol * expDown, dParams[4], 1,
                   &price);
    if (err)
      goto FREE_RETURN;
    volIter++;
  } while (fabs(price - dMarketprices[2]) > tolerance && volIter < iterLimit);

  if (fabs(price - dMarketprices[2]) < tolerance) {
    dParams[3] = vol;
  }

  /* Minimisation */
  dWeights[1] = 0.3;
  dWeights[2] = 0.3; /* Put more weight ATM ?*/
  dWeights[3] = 0.3;

  /* Call Levenberg-Marcquardt */
  lUseparams[1] = 0;
  lUseparams[2] = 0;
  lUseparams[3] = 0;
  lUseparams[4] = 1;
  lUseparams[5] = 1;
  lUseparams[6] = 1;

  err = levenberg_marquardt_select(dStrikes, dMarketprices, dWeights, nProducts,
                                   dParams, lUseparams, nParams, nIter,
                                   SSLGradient, error);

  if (dParams[5] > dParams[6]) {
    *VolUp = dParams[3] * exp(dParams[5]);
    *VolDown = dParams[3] * exp(dParams[6]);
  } else {
    *VolUp = dParams[3] * exp(dParams[6]);
    *VolDown = dParams[3] * exp(dParams[5]);
  }
  *Fwd_Shift = dParams[4];

FREE_RETURN:

  if (dStrikes)
    free_dvector(dStrikes, 1, nProducts);
  if (dMarketprices)
    free_dvector(dMarketprices, 1, nProducts);
  if (dWeights)
    free_dvector(dWeights, 1, nProducts);
  if (dParams)
    free_dvector(dParams, 1, nParams);
  if (lUseparams)
    free_lngvector(lUseparams, 1, nParams);
  return err;
}

Err SSLGradient(double strike, double *ParamSSL, double *price,
                double *gradient, int nbParamsSSL) {
  Err err = NULL;
  double vol, volUp, volDown, volpert1, volpert2;
  double shift, forward, maturity;
  double shiftedpriceshift;
  double shiftedpricevolup;
  double shiftedpricevoldown;
  double bumpSL;
  double bumpVolup;
  double bumpVoldown;

  forward = ParamSSL[1];
  maturity = ParamSSL[2];
  vol = ParamSSL[3];
  shift = ParamSSL[4];
  volpert1 = ParamSSL[5];
  volpert2 = ParamSSL[6];

  /* Bump parameters */
  bumpSL = shift / 1000;
  bumpVolup = volpert1 / 1000;
  bumpVoldown = volpert2 / 1000;

  if (fabs(bumpSL) < 0.0001)
    bumpSL = 0.0001;

  if (fabs(bumpVolup) < 0.0001)
    bumpVolup = 0.0001;

  if (fabs(bumpVoldown) < 0.0001)
    bumpVoldown = 0.0001;

  /* Gradient calculations */
  volUp = vol * exp(volpert1);
  volDown = vol * exp(volpert2);

  err = SSLprice(forward, maturity, strike, volUp, volDown, shift, 1, price);
  if (err)
    goto FREE_RETURN;
  err = SSLprice(forward, maturity, strike, volUp, volDown, shift + bumpSL, 1,
                 &shiftedpriceshift);
  if (err)
    goto FREE_RETURN;
  err = SSLprice(forward, maturity, strike, volUp * exp(bumpVolup), volDown,
                 shift, 1, &shiftedpricevolup);
  if (err)
    goto FREE_RETURN;
  err = SSLprice(forward, maturity, strike, volUp, volDown * exp(bumpVoldown),
                 shift, 1, &shiftedpricevoldown);
  if (err)
    goto FREE_RETURN;

  gradient[1] = 0.0;
  gradient[2] = 0.0;
  gradient[3] = 0.0;
  gradient[4] = (shiftedpriceshift - *price) / bumpSL;
  gradient[5] = (shiftedpricevolup - *price) / bumpVolup;
  gradient[6] = (shiftedpricevoldown - *price) / bumpVoldown;

FREE_RETURN:
  return err;
}

/* =========================================================================
 Get the cumulative distribution
========================================================================= */
Err SSLCumulative(double Forward, double Maturity, double VolUp, double VolDown,
                  double FwdShift, int Npoints, int Nstd, double **Cumulative) {
  Err err = NULL;
  double volATM;
  double *firstGuessGrid = NULL, *firstGuessCumul = NULL;
  int i;
  double temp;

  firstGuessGrid = calloc(Npoints, sizeof(double));
  firstGuessCumul = calloc(Npoints, sizeof(double));

  if (!firstGuessGrid || !firstGuessCumul) {
    err = "SSLCumulative: memory allocation failure";
    goto FREE_RETURN;
  }

  /* ========================================
  FIRST GUESS ON THE GRID
  ======================================== */

  /*Takes the highest vol to construct the grid */
  volATM = fabs(VolUp);
  if (volATM < fabs(VolDown))
    volATM = fabs(VolDown);

  /* Calculates the grid */
  if (FwdShift < 0.0 && fabs(FwdShift) > Forward)
    firstGuessGrid[0] =
        -FwdShift + (Forward + FwdShift) * exp(Nstd * volATM * sqrt(Maturity));
  else
    firstGuessGrid[0] =
        (Forward + FwdShift) * exp(-Nstd * volATM * sqrt(Maturity)) - FwdShift;

  firstGuessCumul[0] = 0.0;

  for (i = 1; i < Npoints; i++) {
    if (FwdShift < 0.0 && fabs(FwdShift) > Forward)
      firstGuessGrid[i] =
          -FwdShift +
          (Forward + FwdShift) * exp(-Nstd * (2.0 * (double)i / Npoints - 1.0) *
                                     volATM * sqrt(Maturity));
    else
      firstGuessGrid[i] =
          (Forward + FwdShift) * exp(Nstd * (2.0 * (double)i / Npoints - 1.0) *
                                     volATM * sqrt(Maturity)) -
          FwdShift;

    err = SSLcumul(Forward, Maturity, firstGuessGrid[i], VolUp, VolDown,
                   FwdShift, &(firstGuessCumul[i]));
    if (err)
      goto FREE_RETURN;
  }

  /* ========================================
  GET THE GRID
  ======================================== */
  for (i = 0; i < Npoints; i++) {
    Cumulative[0][i] = firstGuessCumul[0] +
                       (double)i / (Npoints - 1) *
                           (firstGuessCumul[Npoints - 1] - firstGuessCumul[0]);
    Cumulative[0][i] = interp(firstGuessCumul, firstGuessGrid, Npoints,
                              Cumulative[0][i], 0, &temp);
  }

  for (i = 0; i < Npoints; i++) {
    err = SSLcumul(Forward, Maturity, Cumulative[0][i], VolUp, VolDown,
                   FwdShift, &(Cumulative[1][i]));
    if (err)
      goto FREE_RETURN;
  }

FREE_RETURN:
  if (firstGuessGrid)
    free(firstGuessGrid);
  if (firstGuessCumul)
    free(firstGuessCumul);

  return err;
}

Err BMCumulative(double Forward, double Maturity, double BMbeta, double BMfwd1,
                 double BMfwd2, double BMsig1, double BMsig2, double BMpi,
                 int nPoints, int nStd, double **Cumulative) {
  Err err = NULL;
  double *firstGuessGrid = NULL, *firstGuessCumul = NULL;
  double gridMin1, gridMax1, gridMin2, gridMax2, Nmin, Nmax, ATMvol;
  int i;
  double temp;
  double strikeBump, price;
  double expSigSqrtT;

  strikeBump = Forward / 10000;

  firstGuessGrid = calloc(nPoints, sizeof(double));
  firstGuessCumul = calloc(nPoints, sizeof(double));

  if (!firstGuessGrid || !firstGuessCumul) {
    err = "SSLCumulative: memory allocation failure";
    goto FREE_RETURN;
  }

  /* ========================================
  FIRST GUESS ON THE GRID
  ======================================== */

  /* Calculates the grid */
  gridMin1 = BMfwd1 *
             exp(-nStd * BMsig1 * pow(Forward, BMbeta - 1.0) * sqrt(Maturity));
  gridMin2 = BMfwd2 *
             exp(-nStd * BMsig2 * pow(Forward, BMbeta - 1.0) * sqrt(Maturity));
  gridMax1 = BMfwd1 *
             exp(+nStd * BMsig1 * pow(Forward, BMbeta - 1.0) * sqrt(Maturity));
  gridMax2 = BMfwd2 *
             exp(+nStd * BMsig2 * pow(Forward, BMbeta - 1.0) * sqrt(Maturity));

  if (gridMin1 > gridMin2)
    gridMin1 = gridMin2;
  if (gridMax1 < gridMax2)
    gridMax1 = gridMax2;

  err = BMM_Option_Price_From_States(Forward, Maturity, BMbeta, BMfwd1, BMfwd2,
                                     BMsig1, BMsig2, BMpi, SRT_PUT, &price);

  if (err)
    goto FREE_RETURN;

  err = srt_f_optimpvol(price, Forward, Forward, Maturity, 1.0, SRT_PUT,
                        SRT_LOGNORMAL, &ATMvol);

  if (err)
    goto FREE_RETURN;

  Nmin = 1.0 / (ATMvol * sqrt(Maturity)) * log(gridMin1 / Forward);
  Nmax = 1.0 / (ATMvol * sqrt(Maturity)) * log(gridMax1 / Forward);

  expSigSqrtT = exp((Nmax - Nmin) / (nPoints - 1) * ATMvol * sqrt(Maturity));

  for (i = 0; i < nPoints; i++) {
    if (i == 0)
      firstGuessGrid[i] = gridMin1;
    else
      firstGuessGrid[i] = firstGuessGrid[i - 1] * expSigSqrtT;

    err = BMcumul(Maturity, BMfwd1, BMfwd2, BMsig1, BMsig2, BMpi, BMbeta,
                  firstGuessGrid[i], strikeBump, &(firstGuessCumul[i]));

    if (err)
      goto FREE_RETURN;
  }

  /* ========================================
  GET THE GRID
  ======================================== */
  for (i = 0; i < nPoints; i++) {
    Cumulative[0][i] = firstGuessCumul[0] +
                       (double)i / (nPoints - 1) *
                           (firstGuessCumul[nPoints - 1] - firstGuessCumul[0]);
    Cumulative[0][i] = interp(firstGuessCumul, firstGuessGrid, nPoints,
                              Cumulative[0][i], 0, &temp);
  }

  for (i = 0; i < nPoints; i++) {
    err = BMcumul(Maturity, BMfwd1, BMfwd2, BMsig1, BMsig2, BMpi, BMbeta,
                  Cumulative[0][i], strikeBump, &(Cumulative[1][i]));

    if (err)
      goto FREE_RETURN;
  }

FREE_RETURN:
  if (firstGuessGrid)
    free(firstGuessGrid);
  if (firstGuessCumul)
    free(firstGuessCumul);

  return err;
}

Err SSLcumul(double Forward, double Maturity, double Strike, double VolUp,
             double VolDown, double Shift, double *Cumul) {
  Err err = NULL;
  double stdevUp, stdevDown;
  double d2Up, d2Down;
  int cp;

  if (Shift < 0.0 && fabs(Shift) > Forward) {
    cp = 1;
    if (Strike > fabs(Shift)) {
      err = "SSLcumul: Strike not admissible";
      return err;
    }
  } else
    cp = -1;

  stdevUp = sqrt(Maturity);
  stdevDown = fabs(VolDown) * stdevUp;
  stdevUp = fabs(VolUp) * stdevUp;

  d2Up = (log((Forward + Shift) / (Strike + Shift)) -
          0.5 * VolUp * VolUp * Maturity) /
         stdevUp;
  d2Down = (log((Forward + Shift) / (Strike + Shift)) -
            0.5 * VolDown * VolDown * Maturity) /
           stdevDown;

  *Cumul = 0.5 * (norm(cp * d2Up) + norm(cp * d2Down));

  return err;
}

Err BMcumul(double Maturity, double BMfwd1, double BMfwd2, double BMsig1,
            double BMsig2, double BMpi, double BMbeta, double strike,
            double strikeBump, double *cumul) {
  double priceDown, priceUp;
  double stdev1, stdev2;
  double d2_1, d2_2;

  Err err = NULL;

  if (fabs(BMbeta - 1.0) > 1.0e-10) // No closed form solution
  {
    err = BMM_Option_Price_From_States(strike - strikeBump, Maturity, BMbeta,
                                       BMfwd1, BMfwd2, BMsig1, BMsig2, BMpi,
                                       SRT_PUT, &priceDown);

    if (err)
      goto FREE_RETURN;

    err = BMM_Option_Price_From_States(strike + strikeBump, Maturity, BMbeta,
                                       BMfwd1, BMfwd2, BMsig1, BMsig2, BMpi,
                                       SRT_PUT, &priceUp);

    if (err)
      goto FREE_RETURN;

    *cumul = 0.5 * (priceUp - priceDown) / strikeBump;
  } else {
    stdev2 = sqrt(Maturity);
    stdev1 = fabs(BMsig1) * stdev2;
    stdev2 = fabs(BMsig2) * stdev2;

    d2_1 = (log(BMfwd1 / strike) - 0.5 * BMsig1 * BMsig1 * Maturity) / stdev1;
    d2_2 = (log(BMfwd2 / strike) - 0.5 * BMsig2 * BMsig2 * Maturity) / stdev2;

    *cumul = BMpi * norm(-d2_1) + (1.0 - BMpi) * norm(-d2_2);
  }

FREE_RETURN:
  return err;
}

Err SSLprice(double Forward, double Maturity, double Strike, double VolUp,
             double VolDown, double Shift, int IsCall, double *Price) {
  Err err = NULL;

  if (Shift < 0.0 && fabs(Shift) > Forward) {
    *Price =
        srt_f_optblksch(-Forward - Shift, -Strike - Shift, fabs(VolUp),
                        Maturity, 1.0, IsCall ? SRT_PUT : SRT_CALL, PREMIUM);
    *Price +=
        srt_f_optblksch(-Forward - Shift, -Strike - Shift, fabs(VolDown),
                        Maturity, 1.0, IsCall ? SRT_PUT : SRT_CALL, PREMIUM);
  } else {
    *Price =
        srt_f_optblksch(Forward + Shift, Strike + Shift, fabs(VolUp), Maturity,
                        1.0, IsCall ? SRT_CALL : SRT_PUT, PREMIUM);
    *Price +=
        srt_f_optblksch(Forward + Shift, Strike + Shift, fabs(VolDown),
                        Maturity, 1.0, IsCall ? SRT_CALL : SRT_PUT, PREMIUM);
  }
  *Price *= 0.5;
  return err;
}

Err GenerateRealisations(double Mean1, double Mean2, double **cumulative,
                         int Nsimul, int Npoints,
                         int GetGaussOrMarginal, // 0:Gauss / 1:Marginal
                         int FxIndex, int isLinear, double **matrix) {
  Err err = NULL;
  double test1, test2;
  int i, index1, index2, index3, method;

  for (i = 0; i < Nsimul; i++) {
    matrix[i][0] += Mean1;
    matrix[i][1] += Mean2;

    if (GetGaussOrMarginal) {
      matrix[i][FxIndex] = norm(matrix[i][FxIndex]);

      err = OTCutilsGetIndex(cumulative[1], Npoints, matrix[i][FxIndex],
                             isLinear, &index1, &index2, &index3, &method);

      if (err)
        goto FREE_RETURN;

      err = OTCutilsInterpolate(cumulative[0][index1], cumulative[0][index2],
                                cumulative[0][index3], cumulative[1][index1],
                                cumulative[1][index2], cumulative[1][index3],
                                matrix[i][FxIndex], method, &test1, &test2);

      if (err)
        goto FREE_RETURN;

      if (isLinear) {
        matrix[i][FxIndex] = test1; // max between the two
      } else {
        if (test1 > cumulative[0][index1] && test1 < cumulative[0][index3])
          matrix[i][FxIndex] = test1;
        else if (test2 > cumulative[0][index1] && test2 < cumulative[0][index3])
          matrix[i][FxIndex] = test2;
        else if (matrix[i][FxIndex] > cumulative[1][Npoints - 1])
          matrix[i][FxIndex] =
              cumulative[0][Npoints - 1]; // min between the two
        else if (matrix[i][FxIndex] < cumulative[1][0])
          matrix[i][FxIndex] = cumulative[0][0]; // max between the two
        else {
          err = "no solution in interpolation";
          goto FREE_RETURN;
        }
      }
    }
  }

FREE_RETURN:
  return err;
}

/*	Generate BalSam numbers */
Err OTC_balsam_generation(long nbPaths, /* Must be odd */
                          long nbSteps, double **matrix, long *seed) {
  double *init_gauss = NULL;
  double step, prob;
  long i, j;
  int rand;

  init_gauss = dvector(0, nbPaths - 1);
  nbPaths -= 1;
  nbPaths /= 2;

  if (!init_gauss)
    return "Memory allocation failure in OTC_balsam_generation";

  step = 0.5 / (nbPaths + 1);
  prob = step;

  /* Generation of the fractiles of the gaussian */

  for (i = 0; i < nbPaths; i++) {
    init_gauss[i] = inv_cumnorm_fast(prob);
    init_gauss[nbPaths + i + 1] = -init_gauss[i];
    prob += step;
  }
  init_gauss[nbPaths] = 0.0;

  /* shuffle it */
  nbPaths *= 2;
  nbPaths += 1;
  for (j = 0; j < nbSteps; j++) {
    for (i = 0; i < nbPaths - 1; i++) {
      /* rand = random_int(nbPaths-1-i  , &seed) + i; */
      rand = i + (int)((nbPaths - i) * uniform(seed));
      matrix[i][j] = init_gauss[rand];
      init_gauss[rand] = init_gauss[i];
      init_gauss[i] = matrix[i][j];
    }
    matrix[nbPaths - 1][j] = init_gauss[nbPaths - 1];
  }

  if (init_gauss) {
    free_dvector(init_gauss, 0, nbPaths - 1);
  }

  return NULL;
}

/* ===============================================
        Discount factor
=============================================== */
Err OTCdfBis(double realisation, double lambda, double Tstart, double Tend,
             double Tstar, /* measure Date always equals the OTC Date */
             double *volTimes, double *vol, int nVolTimes, double *df) {
  double Phi, beta, expLambda;
  int i, endIndex;
  Err err = NULL;

  Phi = 0.0;
  endIndex = Get_Index(Tstart, volTimes, nVolTimes);
  expLambda = exp(-lambda * (Tend - Tstart));

  if (lambda < 1e-10) {
    err = "OTCdfBis: lambda must be positive";
    goto FREE_RETURN;
  }
  beta = (1 - expLambda) / lambda;

  if (endIndex == 0) {
    Phi += vol[0] * vol[0] * Phi_Func(2 * lambda, Tstart, 0.0, Tstart);
  } else {
    Phi += vol[0] * vol[0] * Phi_Func(2 * lambda, Tstart, 0.0, volTimes[0]);

    for (i = 1; i < endIndex; i++)
      Phi += vol[i] * vol[i] *
             Phi_Func(2 * lambda, Tstart, volTimes[i - 1], volTimes[i]);

    Phi += vol[endIndex] * vol[endIndex] *
           Phi_Func(2 * lambda, Tstart, volTimes[i - 1], Tstart);
  }

  *df = -0.5 * beta * beta * Phi;
  *df -= (1 - expLambda) * realisation / lambda;
  *df = exp(*df);

FREE_RETURN:
  return err;
}

Err OTCdfPhi(double lambda, double Tstart, double *volTimes, double *vol,
             int nVolTimes, double *Phi) {
  int i, endIndex;
  Err err = NULL;

  *Phi = 0.0;
  endIndex = Get_Index(Tstart, volTimes, nVolTimes);

  if (lambda < 1e-10) {
    err = "OTCdfBis: lambda must be positive";
    goto FREE_RETURN;
  }

  if (endIndex == 0) {
    *Phi += vol[0] * vol[0] * Phi_Func(2 * lambda, Tstart, 0.0, Tstart);
  } else {
    *Phi += vol[0] * vol[0] * Phi_Func(2 * lambda, Tstart, 0.0, volTimes[0]);

    for (i = 1; i < endIndex; i++)
      *Phi += vol[i] * vol[i] *
              Phi_Func(2 * lambda, Tstart, volTimes[i - 1], volTimes[i]);

    *Phi += vol[endIndex] * vol[endIndex] *
            Phi_Func(2 * lambda, Tstart, volTimes[i - 1], Tstart);
  }

FREE_RETURN:
  return err;
}

Err OTCdf(double realisation, double lambda, double Tstart, double Tend,
          double Tstar, /* measure Date always equals the OTC Date */
          double Phi, double *df) {
  double beta, expLambda;
  Err err = NULL;

  expLambda = exp(-lambda * (Tend - Tstart));

  if (lambda < 1e-10) {
    err = "OTCdfBis: lambda must be positive";
    goto FREE_RETURN;
  }
  beta = (1 - expLambda) / lambda;

  *df = -0.5 * beta * beta * Phi;
  *df -= (1 - expLambda) * realisation / lambda;
  *df = exp(*df);

FREE_RETURN:
  return err;
}

Err OTCgetMoments(double endTime, double ldaDom, double ldaFor, double TStar,
                  double *mergeTimes, int mergeNtimes, double *SigmaDom,
                  double *SigmaFor, double *SigmaFx, double *CorrDomFor,
                  double *CorrDomFx, double *CorrForFx,

                  double *dom_std, double *for_fwd, double *for_std,

                  double *dom_for_cov, double *dom_fx_cov, double *for_fx_cov) {
  int j;
  double t1, t2;
  double var_rates_dom, var_rates_for, var_ffx;
  double covar_dom_for, covar_dom_ffx, covar_for_ffx;
  double sig_dom, sig_for, sig_fx, sig_dom2, sig_for2, sig_fx2, sig_domfor,
      sig_domfx, sig_forfx;

  int start_index, end_index;

  Err err = NULL;

  if (fabs(endTime) < 1.0e-10)
    endTime = 1.0e-10;

  start_index = Get_Index(0.0, mergeTimes, mergeNtimes);
  end_index = Get_Index(endTime, mergeTimes, mergeNtimes);

  var_ffx = 0.0;
  var_rates_dom = 0.0;
  var_rates_for = 0.0;
  covar_dom_for = 0.0;
  covar_dom_ffx = 0.0;
  covar_for_ffx = 0.0;

  for (j = start_index; j < end_index + 1; j++) {
    if (j > start_index)
      t1 = mergeTimes[j - 1];
    else
      t1 = 0.0; /* First part */

    if (j == end_index || start_index == end_index)
      t2 = endTime; /* Last part */
    else
      t2 = mergeTimes[j];

    sig_dom = SigmaDom[j];
    sig_for = SigmaFor[j];
    sig_fx = SigmaFx[j];
    sig_dom2 = sig_dom * sig_dom;
    sig_for2 = sig_for * sig_for;
    sig_fx2 = sig_fx * sig_fx;
    sig_domfor = CorrDomFor[j] * sig_dom * sig_for;
    sig_domfx = CorrDomFx[j] * sig_dom * sig_fx;
    sig_forfx = CorrForFx[j] * sig_for * sig_fx;

    /* Expectation and Variance */
    var_rates_dom += sig_dom2 * Phi_Func(2.0 * ldaDom, TStar, t1, t2);
    var_rates_for += sig_for2 * Phi_Func(2.0 * ldaFor, TStar, t1, t2);

    covar_dom_for += sig_domfor * Phi_Func(ldaDom + ldaFor, TStar, t1, t2);

    covar_dom_ffx += sig_domfx * Phi_Func(ldaDom, TStar, t1, t2) -
                     sig_domfor * Gamma_Func(ldaFor, ldaDom, TStar, t1, t2) +
                     sig_dom2 * Gamma_Func(ldaDom, ldaDom, TStar, t1, t2);

    covar_for_ffx += sig_forfx * Phi_Func(ldaFor, TStar, t1, t2) -
                     sig_for2 * Gamma_Func(ldaFor, ldaFor, TStar, t1, t2) +
                     sig_domfor * Gamma_Func(ldaDom, ldaFor, TStar, t1, t2);

    /* Fx component */
    var_ffx += sig_fx2 * (t2 - t1);

    /* Fx / Rates component */
    var_ffx += -2.0 * sig_forfx * Etha_Func(ldaFor, TStar, t1, t2) +
               2.0 * sig_domfx * Etha_Func(ldaDom, TStar, t1, t2);

    /* Rates / Rates component */
    var_ffx += sig_for2 * Psi_Func(ldaFor, ldaFor, TStar, t1, t2) +
               sig_dom2 * Psi_Func(ldaDom, ldaDom, TStar, t1, t2) -
               2.0 * sig_domfor * Psi_Func(ldaDom, ldaFor, TStar, t1, t2);
  }

  *for_fwd = -covar_for_ffx;

  *for_std = var_rates_for;
  *dom_std = var_rates_dom;

  *dom_std = sqrt(*dom_std);
  *for_std = sqrt(*for_std);

  *dom_for_cov = covar_dom_for / *dom_std / *for_std;
  *dom_fx_cov = covar_dom_ffx / *dom_std / sqrt(var_ffx);
  *for_fx_cov = covar_for_ffx / *for_std / sqrt(var_ffx);

  return err;
}

/* ====================================
Convexity adjustment for Fx(Tval) paid at Tpay
==================================== */
Err OTCcvxtAdj(double t, double Tval, double Tpay, double ldaDom, double ldaFor,
               double *mergeTimes, int mergeNtimes, double *SigmaDom,
               double *SigmaFor, double *SigmaFx, double *CorrDomFor,
               double *CorrDomFx, double *CorrForFx,
               /* Results */
               double *adjustment) {
  int j;
  double t1, t2;
  double covar_dom_ffx;
  double sig_dom, sig_for, sig_fx, sig_dom2, sig_domfor, sig_domfx;

  int start_index, end_index;

  Err err = NULL;

  covar_dom_ffx = 0.0;

  if (fabs(Tval - Tpay) > 1.0e-10) {
    start_index = Get_Index(t, mergeTimes, mergeNtimes);
    end_index = Get_Index(Tval, mergeTimes, mergeNtimes);

    for (j = start_index; j < end_index + 1; j++) {
      if (j > start_index)
        t1 = mergeTimes[j - 1];
      else
        t1 = t; /* First part */

      if (j == end_index || start_index == end_index)
        t2 = Tval; /* Last part */
      else
        t2 = mergeTimes[j];

      sig_dom = SigmaDom[j];
      sig_for = SigmaFor[j];
      sig_fx = SigmaFx[j];
      sig_dom2 = sig_dom * sig_dom;
      sig_domfor = CorrDomFor[j] * sig_dom * sig_for;
      sig_domfx = CorrDomFx[j] * sig_dom * sig_fx;

      /* Cumulative Covariance */
      covar_dom_ffx +=
          sig_domfx * Phi_Func(ldaDom, 0.0, t1, t2) -
          sig_domfor * Gamma2_Func(ldaFor, ldaDom, Tval, 0.0, t1, t2) +
          sig_dom2 * Gamma2_Func(ldaDom, ldaDom, Tval, 0.0, t1, t2);
    }
    covar_dom_ffx *= 1.0 / ldaDom * (exp(-ldaDom * Tval) - exp(-ldaDom * Tpay));
  }

  *adjustment = exp(covar_dom_ffx);

  return err;
}

Err OTCcorrelateVariables(double domStd, double forStd, double corDF,
                          double corDFx, double corFFx, long nbPaths,
                          double **inputMatrix, double **outputMatrix) {
  long i;
  Err err = NULL;

  double b, c, d, e, f, x, y, z;

  b = corDF;
  c = 1 - b * b;
  c = sqrt(c);

  d = corDFx;

  if (c != 0) {
    e = (corFFx - b * d) / c;
    f = 1 - d * d - e * e;
    if (f < 0) {
      return "Check your correlations";
    }
    f = sqrt(f);
  } else {
    if (fabs(corDFx - corFFx) < 1e-10) {
      e = 1 - d * d;
      f = 0;
      if (e < 0) {
        return "Check your correlations";
      }
      e = sqrt(e);
    } else {
      return "Check your correlations";
    }
  }

  for (i = 0; i < nbPaths; i++) {
    x = inputMatrix[i][0];
    y = inputMatrix[i][1];
    z = inputMatrix[i][2];

    outputMatrix[i][0] = domStd * x;
    outputMatrix[i][1] = forStd * (b * x + c * y);
    outputMatrix[i][2] = d * x + e * y + f * z;
  }

  return err;
}

/* Forward SABR parameters */
Err OTCgetFwdSABR(double *forward, CPD_STR cpd, CPD_UND und,
                  SMILE_VOL_MARKET smile_mkt, int otc, double *optionsSigma,
                  double *optionsAlpha, double *optionsBeta, double *optionsRho,
                  double *optionsPi, int noptions, double *mergeTimes,
                  int mergeNtimes, double *sigDom, double ldaDom,
                  double *sigFor, double ldaFor, double *sigFx, double *corDF,
                  double *corDFx, double *corFFx, char *CPDsigma,
                  char *CPDalpha, char *CPDbeta, char *CPDrho, int method)
/* 0=Sliding 3F vol; 1=Converging 3F vol; 2=Sliding Cvg Sbeta; 3=Cvg Cvg Sbeta;
   4=3F; 5=FwdSmile */
{
  int i, j;
  double power;
  long FxNotValDate[2];
  double FxNotValTime[2];
  double value;

  char *CPDvol[4];

  SMILE_PARAMETERS smile_params = NULL;

  Err err = NULL;

  for (i = 0; i < 2; i++) {
    if (otc + i < cpd->num_calls) {
      FxNotValDate[i] = add_unit(cpd->call[otc + i].ex_date, 2, SRT_BDAY,
                                 MODIFIED_SUCCEEDING);
      FxNotValTime[i] = (FxNotValDate[i] - und->today) * YEARS_IN_DAY;
    }
  }

  if (method == 5) {
    CPDvol[0] = CPDsigma;
    CPDvol[1] = CPDalpha;
    CPDvol[2] = CPDbeta;
    CPDvol[3] = CPDrho;

    for (i = 0; i < noptions; i++) {
      for (j = 0; j < 4; j++) {
        err = swp_f_vol(CPDvol[j], cpd->call[otc].ex_date,
                        (i == (noptions - 1)) ? cpd->pd_leg->not_ref.fx_fix_date
                                              : cpd->pd_leg->cpn[i].fx_fix_date,
                        (i == (noptions - 1)) ? cpd->pd_leg->not_ref.fx_val_date
                                              : cpd->pd_leg->cpn[i].fx_val_date,
                        &value, &power);

        if (err)
          goto FREE_RETURN;

        if (j == 0)
          optionsSigma[i] = value;
        else if (j == 1)
          optionsAlpha[i] = value;
        else if (j == 2)
          optionsBeta[i] = value;
        else if (j == 3)
          optionsRho[i] = value;
      }

      err = srt_f_optsarbvol(
          forward[i], forward[i],
          (i == (noptions - 1)) ? cpd->pd_leg->not_ref.fx_fix_time
                                : cpd->pd_leg->cpn[i].fx_fix_time,
          optionsSigma[i], optionsAlpha[i], optionsBeta[i], optionsRho[i],
          SRT_LOGNORMAL, SRT_BETAVOL, &(optionsSigma[i]));

      if (err)
        goto FREE_RETURN;
    }

    for (i = 0; i < 2; i++) {
      if (otc + i < cpd->num_calls) {
        for (j = 0; j < 4; j++) {
          err = swp_f_vol(CPDvol[j], cpd->call[otc + i].ex_date,
                          cpd->call[otc + i].ex_date, FxNotValDate[i], &value,
                          &power);

          if (err)
            goto FREE_RETURN;

          if (j == 0)
            optionsSigma[noptions + i] = value;
          else if (j == 1)
            optionsAlpha[noptions + i] = value;
          else if (j == 2)
            optionsBeta[noptions + i] = value;
          else if (j == 3)
            optionsRho[noptions + i] = value;
        }

        err = srt_f_optsarbvol(
            forward[noptions + i], forward[noptions + i],
            cpd->call[otc + i].ex_time, optionsSigma[noptions + i],
            optionsAlpha[noptions + i], optionsBeta[noptions + i],
            optionsRho[noptions + i], SRT_LOGNORMAL, SRT_BETAVOL,
            &(optionsSigma[noptions + i]));

        if (err)
          goto FREE_RETURN;
      }
    }
  } else {
    smile_params = calloc(1, sizeof(smile_parameters));

    for (i = 0; i < noptions; i++) {
      err = OTCfwdSABR(forward ? forward[i] : 0.0, cpd->call[otc].ex_time,
                       (i == (noptions - 1)) ? cpd->pd_leg->not_ref.fx_fix_time
                                             : cpd->pd_leg->cpn[i].fx_fix_time,
                       (i == (noptions - 1)) ? cpd->pd_leg->not_ref.fx_val_time
                                             : cpd->pd_leg->cpn[i].fx_val_time,
                       smile_mkt, smile_params, mergeTimes, mergeNtimes, sigDom,
                       ldaDom, sigFor, ldaFor, sigFx, corDF, corDFx, corFFx,
                       method);

      if (err)
        goto FREE_RETURN;

      optionsSigma[i] = smile_params->sigma;
      optionsAlpha[i] = smile_params->alpha;
      optionsBeta[i] = smile_params->beta;
      optionsRho[i] = smile_params->rho;

      if (smile_params->smile_spec_type == 2 ||
          smile_params->smile_spec_type == 3)
        optionsPi[i] = smile_params->pi;
    }

    for (i = 0; i < 2; i++) {
      if (otc + i < cpd->num_calls) {
        err = OTCfwdSABR(forward ? forward[noptions + i] : 0.0,
                         cpd->call[otc].ex_time, cpd->call[otc + i].ex_time,
                         FxNotValTime[i], smile_mkt, smile_params, mergeTimes,
                         mergeNtimes, sigDom, ldaDom, sigFor, ldaFor, sigFx,
                         corDF, corDFx, corFFx, method);

        if (err)
          goto FREE_RETURN;

        optionsSigma[noptions + i] = smile_params->sigma;
        optionsAlpha[noptions + i] = smile_params->alpha;
        optionsBeta[noptions + i] = smile_params->beta;
        optionsRho[noptions + i] = smile_params->rho;

        if (smile_params->smile_spec_type == 2 ||
            smile_params->smile_spec_type == 3)
          optionsPi[noptions + i] = smile_params->pi;
      }
    }
  }

FREE_RETURN:

  if (smile_params)
    free(smile_params);

  return err;
}

Err OTCfwdSABR(double forward, double fwdStartTime, double optionsTimeFix,
               double optionsTimeSettle, SMILE_VOL_MARKET smile_mkt,
               SMILE_PARAMETERS smile_params, double *mergeTimes,
               int mergeNtimes, double *sigDom, double ldaDom, double *sigFor,
               double ldaFor, double *sigFx, double *corDF, double *corDFx,
               double *corFFx, int method) {
  double timeToMaturity, logVol;
  Err err = NULL;

  if (optionsTimeFix > fwdStartTime) {
    if (method == 0 || method == 2) // Sliding Alpha Beta and Rho
      timeToMaturity = optionsTimeFix - fwdStartTime;
    else // Converging Alpha Beta and Rho
      timeToMaturity = optionsTimeFix;

    if (method == 4) {
      smile_params->alpha = 0.0;
      smile_params->beta = 1.0;
      smile_params->rho = 0.0;
      smile_params->pi = 0.5;
    } else {
      err = cpd_vol_get_smile_params(1, 0, timeToMaturity, smile_mkt,
                                     smile_params);
      if (err)
        goto FREE_RETURN;
    }

    if (method < 2) {
      err = Fx3DtsImpliedVol_corr(optionsTimeSettle, fwdStartTime,
                                  optionsTimeFix, mergeTimes, (long)mergeNtimes,
                                  sigDom, ldaDom, sigFor, ldaFor, mergeTimes,
                                  sigFx, (long)mergeNtimes, mergeTimes, corDF,
                                  corDFx, corFFx, (long)mergeNtimes, &logVol);

      if (err)
        goto FREE_RETURN;

      smile_params->sigma = logVol;

      err = cpd_vol_get_vol(forward, optionsTimeFix, forward, smile_params,
                            SRT_BETAVOL, &(smile_params->sigmabeta));
      if (err)
        goto FREE_RETURN;

      smile_params->sigma = smile_params->sigmabeta;
    } else if (method == 2 || method == 3) {
      err = Fx3DtsImpliedVol_corr(
          optionsTimeSettle, 0.0, optionsTimeFix, mergeTimes, (long)mergeNtimes,
          sigDom, ldaDom, sigFor, ldaFor, mergeTimes, sigFx, (long)mergeNtimes,
          mergeTimes, corDF, corDFx, corFFx, (long)mergeNtimes, &logVol);

      if (err)
        goto FREE_RETURN;

      smile_params->sigma = logVol;

      err = cpd_vol_get_vol(forward, optionsTimeFix, forward, smile_params,
                            SRT_BETAVOL, &(smile_params->sigmabeta));
      if (err)
        goto FREE_RETURN;

      smile_params->sigma = smile_params->sigmabeta;
    } else if (method == 1234 || method == 4) {
      err = Fx3DtsImpliedVol_corr(optionsTimeSettle, fwdStartTime,
                                  optionsTimeFix, mergeTimes, (long)mergeNtimes,
                                  sigDom, ldaDom, sigFor, ldaFor, mergeTimes,
                                  sigFx, (long)mergeNtimes, mergeTimes, corDF,
                                  corDFx, corFFx, (long)mergeNtimes, &logVol);

      if (err)
        goto FREE_RETURN;

      smile_params->sigma = logVol;
    }
  } else {
    smile_params->sigma = 0.0;
    smile_params->alpha = 0.0;
    smile_params->beta = 1.0;
    smile_params->rho = 0.0;
  }

FREE_RETURN:
  return err;
}

Err OTCgetSMILEparams(CPD_STR cpd, CPD_UND und, double Tfix, double Tval,
                      long DateVal, SMILE_VOL_MARKET smile_mkt,
                      SMILE_PARAMETERS smile_params, int get_full_smile) {
  Err err = NULL;

  /* ==================================
  SABR parameters @ Tval
  ================================== */
  err = Fx3DtsImpliedVol_corr(
      Tval, 0.0, Tfix, und->sigma_time_rates, und->sigma_n_rates,
      und->sigma_dom, und->lda_dom, und->sigma_for, und->lda_for,
      und->sigma_time_fx, und->sigma_fx, und->sigma_n_fx, und->corr_times,
      und->correl_dom_for, und->correl_dom_fx, und->correl_for_fx,
      und->corr_n_times, &(smile_params->sigma));

  if (err)
    goto FREE_RETURN;

  if (get_full_smile) {
    err = cpd_vol_get_smile_params(0, DateVal, 0.0, smile_mkt, smile_params);
    if (err)
      goto FREE_RETURN;
  } else {
    smile_params->alpha = 0.0;
    smile_params->beta = 1.0;
    smile_params->rho = 0.0;
    smile_params->pi = 0.5;
  }

FREE_RETURN:
  return err;
}

Err OTCgetSABRparams(CPD_STR cpd, CPD_UND und, double Tfix, double Tval,
                     long DateVal, long *fx_mkt_vol_date,
                     double *fx_mkt_smile_alpha, double *fx_mkt_smile_beta,
                     double *fx_mkt_smile_rho, int num_fx_mkt_vol,
                     double *smileVol, double *smileAlpha, double *smileBeta,
                     double *smileRho, int smileFee) {
  double interpCoef;
  int j;
  Err err = NULL;

  /* ==================================
  SABR parameters @ Tval
  ================================== */
  err = Fx3DtsImpliedVol_corr(Tval, 0.0, Tfix, und->sigma_time_rates,
                              und->sigma_n_rates, und->sigma_dom, und->lda_dom,
                              und->sigma_for, und->lda_for, und->sigma_time_fx,
                              und->sigma_fx, und->sigma_n_fx, und->corr_times,
                              und->correl_dom_for, und->correl_dom_fx,
                              und->correl_for_fx, und->corr_n_times, smileVol);

  if (err)
    goto FREE_RETURN;

  if (smileFee) {
    if (DateVal <= fx_mkt_vol_date[0]) {
      *smileAlpha = fx_mkt_smile_alpha[0];
      *smileBeta = fx_mkt_smile_beta[0];
      *smileRho = fx_mkt_smile_rho[0];
    } else if (DateVal >= fx_mkt_vol_date[num_fx_mkt_vol - 1]) {
      *smileAlpha = fx_mkt_smile_alpha[num_fx_mkt_vol - 1];
      *smileBeta = fx_mkt_smile_beta[num_fx_mkt_vol - 1];
      *smileRho = fx_mkt_smile_rho[num_fx_mkt_vol - 1];
    } else {
      for (j = 0; DateVal >= fx_mkt_vol_date[j + 1]; j++)
        ;
      interpCoef = (double)(DateVal - fx_mkt_vol_date[j]) /
                   (double)(fx_mkt_vol_date[j + 1] - fx_mkt_vol_date[j]);

      *smileAlpha = interpCoef * fx_mkt_smile_alpha[j + 1] +
                    (1.0 - interpCoef) * fx_mkt_smile_alpha[j];
      *smileBeta = interpCoef * fx_mkt_smile_beta[j + 1] +
                   (1.0 - interpCoef) * fx_mkt_smile_beta[j];
      *smileRho = interpCoef * fx_mkt_smile_rho[j + 1] +
                  (1.0 - interpCoef) * fx_mkt_smile_rho[j];
    }
  } else {
    *smileAlpha = 0.0;
    *smileBeta = 1.0;
    *smileRho = 0.0;
  }

FREE_RETURN:
  return err;
}

/*
Calibration of a smile model
*/
Err OTCcalibSmileModel(double Tex,     /* Maturity in years */
                       double forward, /* Forward */
                       otcpd_params *OTC_params, SMILE_PARAMETERS smile_params,
                       double *error) {
  double StrikeUp, StrikeDown, VolUp, VolDown, alpha, rho, CalibError,
      SigmaBeta;
  Err err = NULL;

  if (OTC_params->smileModel == 0) // supposes that the market entered is SABR
  {
    if (smile_params->smile_spec_type > 1) {
      err = "SSL not supported with BMM MKT parameters";
      goto FREE_RETURN;
    }
    err = CalibSSLtoSABR(Tex,                /* Maturity in years */
                         forward,            /* Forward */
                         OTC_params->SSLstd, /* # std to calibrate skew */
                         /* SABR Parameters */
                         smile_params->sigma,
                         OTC_params->OTCsmile ? smile_params->alpha : 0.0,
                         OTC_params->OTCsmile ? smile_params->beta : 1.0,
                         OTC_params->OTCsmile ? smile_params->rho : 0.0,
                         /* Results */
                         &(OTC_params->SSLvolUp), &(OTC_params->SSLvolDown),
                         &(OTC_params->SSLshift), error,
                         /* Parameters */
                         OTC_params->SSLnIter);

    if (err)
      goto FREE_RETURN;

    if (fabs(*error) > 1.0) {
      err = "otc_caller: SSL calib failed";
      goto FREE_RETURN;
    }

  } else if (OTC_params->smileModel == 1 || OTC_params->smileModel == 2) {
    if (smile_params->smile_spec_type > 1) {
      err = srt_f_optbmmvol(forward, forward, Tex, smile_params->sigma,
                            smile_params->alpha, smile_params->beta,
                            smile_params->rho, smile_params->pi, SRT_LOGNORMAL,
                            SRT_BETAVOL, &SigmaBeta);

      if (err)
        goto FREE_RETURN;

      err = BMMGetStates(
          forward, Tex, SigmaBeta, smile_params->alpha, smile_params->beta,
          smile_params->rho, smile_params->pi, &(OTC_params->BMsig1),
          &(OTC_params->BMsig2), &(OTC_params->BMfwd1), &(OTC_params->BMfwd2));

      if (err)
        goto FREE_RETURN;

      OTC_params->BMbeta = smile_params->beta;
      OTC_params->BMpi = smile_params->pi;
    } else {
      if (OTC_params->OTCsmile && OTC_params->smileModel == 1) {
        alpha = smile_params->alpha;
        rho = smile_params->rho;

        StrikeDown = forward *
                     exp(-OTC_params->SSLstd * smile_params->sigma * sqrt(Tex));
        StrikeUp = forward *
                   exp(+OTC_params->SSLstd * smile_params->sigma * sqrt(Tex));

        err = srt_f_optsarbvol(forward, StrikeDown, Tex, smile_params->sigma,
                               smile_params->alpha, smile_params->beta,
                               smile_params->rho, SRT_LOGNORMAL, SRT_LOGNORMAL,
                               &VolDown);

        if (err)
          goto FREE_RETURN;

        err = srt_f_optsarbvol(forward, StrikeUp, Tex, smile_params->sigma,
                               smile_params->alpha, smile_params->beta,
                               smile_params->rho, SRT_LOGNORMAL, SRT_LOGNORMAL,
                               &VolUp);

        if (err)
          goto FREE_RETURN;

        err = calib_sabr_rr_bt_given_beta(
            forward, Tex, smile_params->sigma, StrikeDown, VolDown, StrikeUp,
            VolUp, &SigmaBeta, &alpha, 1.0, &rho, SRT_LOGNORMAL, 1, 0.0001, 20,
            &CalibError);

        if (!err) {
          smile_params->alpha = alpha;
          smile_params->rho = rho;
          smile_params->beta = 1.0;
        }
      }

      err = BMMCalibOnSabrStates(
          forward, Tex, smile_params->sigma,
          OTC_params->OTCsmile ? smile_params->beta : 1.0,
          OTC_params->OTCsmile ? smile_params->alpha : 0.0,
          OTC_params->OTCsmile ? smile_params->rho : 0.0, OTC_params->BMpi,
          OTC_params->SSLstd, &(OTC_params->BMfwd1), &(OTC_params->BMfwd2),
          &(OTC_params->BMsig1), &(OTC_params->BMsig2), &(OTC_params->BMpi),
          SRT_LOGNORMAL, error);

      if (OTC_params->OTCsmile)
        OTC_params->BMbeta = smile_params->beta;
      else
        OTC_params->BMbeta = 1.0;

      if (*error > 1.0e-3) {
        err = "OTCcalibSmileModel: BiBetaCalibOnSabr calibration failed";
      }

      if (err)
        goto FREE_RETURN;
    }
  } else {
  }

FREE_RETURN:
  return err;
}

Err OTCgetCumulative(double forward, double Tex, otcpd_params *OTC_params) {
  double *dX = NULL, *dCumulative = NULL;
  long l, lNbPoints;
  Err err = NULL;

  if (OTC_params->smileModel == 0) {
    err = SSLCumulative(forward, Tex, OTC_params->SSLvolUp,
                        OTC_params->SSLvolDown, OTC_params->SSLshift,
                        OTC_params->CUMULnPoints, OTC_params->CUMULnStd,
                        OTC_params->CUMUL);

    if (err)
      goto FREE_RETURN;
  } else if (OTC_params->smileModel == 1 || OTC_params->smileModel == 2) {
    lNbPoints = OTC_params->CUMULnPoints;
    dX = calloc(lNbPoints, sizeof(double));
    dCumulative = calloc(lNbPoints, sizeof(double));

    err = copula_gaussian_get_BMM_linterp_cumulative(
        Tex, forward, OTC_params->BMfwd1, OTC_params->BMsig1,
        OTC_params->BMfwd2, OTC_params->BMsig2, OTC_params->BMbeta,
        OTC_params->BMpi, OTC_params->CUMULPrecision, OTC_params->CUMULnStd,
        OTC_params->CUMULnStd, &(lNbPoints), &dX, &dCumulative);

    if (err)
      goto FREE_RETURN;

    err = OTCmatrixRealloc(&(OTC_params->CUMUL), 0, 2 - 1, 0,
                           OTC_params->CUMULnPoints - 1, 0, 2 - 1, 0,
                           lNbPoints - 1);

    OTC_params->CUMULnPoints = lNbPoints;

    for (l = 0; l < lNbPoints; l++) {
      OTC_params->CUMUL[0][l] = dX[l];
      OTC_params->CUMUL[1][l] = dCumulative[l];
    }
  } else {
  }

FREE_RETURN:
  if (dX)
    free(dX);
  if (dCumulative)
    free(dCumulative);
  return err;
}

Err OTCmatrixRealloc(double ***mat, long nrl, long nrh, long ncl, long nch,
                     long new_nrl, long new_nrh, long new_ncl, long new_nch) {
  Err err = NULL;

  if (*mat) {
    free_dmatrix(*mat, nrl, nrh, ncl, nch);
    *mat = dmatrix(new_nrl, new_nrh, new_ncl, new_nch);
  } else {
    err = "matrix not allocated";
    goto FREE_RETURN;
  }

FREE_RETURN:
  return err;
}

Err OTCgetChangedPayoff(CPD_STR cpd, CPD_UND und, double *NewCoupons) {
  int i, k, pdIndex, startIndex, endIndex;

  Err err = NULL;

  if (cpd->fund_leg->dom_for == 0) {
    for (k = 0; k < cpd->num_calls; k++) {
      pdIndex = cpd->call[k].pd_idx;
      startIndex = cpd->call[k].fund_idx;

      if (k < cpd->num_calls - 1)
        endIndex = cpd->call[k + 1].fund_idx;
      else
        endIndex = cpd->fund_leg->num_cpn;

      for (i = startIndex; i < endIndex; i++) {
        NewCoupons[pdIndex] +=
            cpd->fund_leg->cpn[i].cpn *
            swp_f_df(und->today, cpd->fund_leg->cpn[i].pay_date, und->dom_yc);
      }

      NewCoupons[pdIndex] /=
          swp_f_df(und->today, cpd->pd_leg->cpn[pdIndex].pay_date, und->dom_yc);
    }
  } else {
    for (k = 0; k < cpd->num_calls; k++) {
      pdIndex = cpd->call[k].pd_idx;
      startIndex = cpd->call[k].fund_idx;

      if (k < cpd->num_calls - 1)
        endIndex = cpd->call[k + 1].fund_idx;
      else
        endIndex = cpd->fund_leg->num_cpn;

      for (i = startIndex; i < endIndex; i++) {
        NewCoupons[pdIndex] +=
            cpd->fund_leg->cpn[i].cpn *
            swp_f_df(und->today, cpd->fund_leg->cpn[i].pay_date, und->for_yc);
      }

      NewCoupons[pdIndex] /= swp_f_df(
          und->today, cpd->pd_leg->cpn[pdIndex].fx_val_date, und->for_yc);
    }
  }
  return err;
}

Err OTCgetChangedPayoff2(CPD_STR cpd, CPD_UND und, double *NewCoupons) {
  int i, j, k, couponIndex;
  long midDate;

  Err err = NULL;

  for (i = 0; i < cpd->fund_leg->num_cpn; i++) {
    k = 0;
    while (i >= cpd->call[k].fund_idx) // check if >=
      k++;

    j = 0;
    while (cpd->fund_leg->cpn[i].pay_date > cpd->pd_leg->cpn[j].pay_date)
      j++;

    if (cpd->fund_leg->cpn[i].pay_date == cpd->pd_leg->cpn[j].pay_date) {
      NewCoupons[j] += cpd->fund_leg->cpn[i].cpn;
    } else {
      if (!k) {
        if (!j) {
          if (cpd->pd_leg->cpn[j].pay_date > cpd->call[k].ex_date) {
            // put the funding to today
          } else
            couponIndex = j; // put the funding on the first pwd coupon : j
        } else {
          if (cpd->pd_leg->cpn[j].pay_date > cpd->call[k].ex_date)
            couponIndex = j - 1; // put the funding of the j-1
          else {
            // put the funding on the nearest pwd coupon
            midDate = (long)((cpd->pd_leg->cpn[j - 1].pay_date +
                              cpd->pd_leg->cpn[j].pay_date) /
                             2.0);
            if (cpd->fund_leg->cpn[i].pay_date > midDate)
              couponIndex = j; // put the funding on the j
            else
              couponIndex = j - 1; // put the funding on the j - 1
          }
        }
      } else {
        if (cpd->pd_leg->cpn[j].start_date > cpd->call[k].ex_date)
          couponIndex = j; // put all on the j-1 (test if j-1_date > k-1_date)
        else {
          if (cpd->pd_leg->cpn[j - 1].pay_date > cpd->call[k - 1].ex_date) {
            // put the funding on the nearest pwd coupon
            midDate = (long)((cpd->pd_leg->cpn[j - 1].pay_date +
                              cpd->pd_leg->cpn[j].pay_date) /
                             2.0);
            if (cpd->fund_leg->cpn[i].pay_date > midDate)
              couponIndex = j; // put the funding on the j
            else
              couponIndex = j - 1; // put the funding on the j - 1
          } else
            couponIndex = j; // put the funding on the j
        }
      }

      NewCoupons[couponIndex] +=
          cpd->fund_leg->cpn[i].cpn *
          swp_f_df(und->today, cpd->fund_leg->cpn[i].pay_date, und->dom_yc) /
          swp_f_df(und->today, cpd->pd_leg->cpn[couponIndex].pay_date,
                   und->dom_yc);
    }
  }

  return err;
}

/* =================================================================================================
   =================================================================================================

                PAYOFF FUNCTIONS

=================================================================================================
=================================================================================================
*/
Err OTCfwdSABRCalib(CPD_STR cpd, CPD_UND und, double *pd_fwdsmile_vol,
                    double *pd_fwdsmile_alpha, double *pd_fwdsmile_beta,
                    double *pd_fwdsmile_rho, int otc, int nSimul,
                    double **matrix, double tstarTime, int tstarDate,
                    double fwdSmileVisuNstd, double *FxAdj, double **Results) {

  int i, j, realisation;

  /* For the payoff */
  double *strikeUp = NULL, *strikeATM = NULL, *strikeDown = NULL;
  double *priceUp = NULL, *priceATM = NULL, *priceDown = NULL;
  double *forward, vol, localForward;
  int type = 3; /* CALL */
  double smile_std;
  double smile_half_std;

  PD_EXO_CPN cpn;

  double *PdDfDom = NULL, *FwdDfDom = NULL, *FwdDfFor = NULL;
  double df, dfDom, dfFor;

  /* For the call */
  int pd_otc_start; /*	index of the first coupon to be called */

  /* Positions */
  int sUp = 0, sATM = 2, sDown = 4, vUp = 1, vATM = 3, vDown = 5;

  /* For the Df */
  double phiDom, phiFor;

  Err err = NULL;

  /* ==================================
  memory allocation
  ================================== */
  PdDfDom = calloc(cpd->pd_leg->num_cpn, sizeof(double));
  FwdDfDom = calloc(cpd->pd_leg->num_cpn, sizeof(double));
  FwdDfFor = calloc(cpd->pd_leg->num_cpn, sizeof(double));
  forward = calloc(cpd->pd_leg->num_cpn, sizeof(double));

  if (!PdDfDom || !FwdDfDom || !FwdDfFor || !forward) {
    err = "OTCcaller: memory allocation error";
    goto FREE_RETURN;
  }

  /* ==================================
  Precalculations
  ================================== */
  pd_otc_start = cpd->call[otc].pd_idx;

  dfDom = swp_f_df(und->today, cpd->call[otc].ex_date, und->dom_yc);
  dfFor = swp_f_df(und->today, cpd->call[otc].ex_date, und->for_yc);

  for (i = pd_otc_start; i < cpd->pd_leg->num_cpn; i++) {
    forward[i] =
        FxAdj[i] * und->spot_fx *
        swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->for_yc) /
        swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->dom_yc);

    err = Fx3DtsImpliedVol_corr(
        cpd->pd_leg->cpn[i].fx_val_time, 0.0, cpd->pd_leg->cpn[i].fx_val_time,
        und->sigma_time_rates, und->sigma_n_rates, und->sigma_dom, und->lda_dom,
        und->sigma_for, und->lda_for, und->sigma_time_fx, und->sigma_fx,
        und->sigma_n_fx, und->corr_times, und->correl_dom_for,
        und->correl_dom_fx, und->correl_for_fx, und->corr_n_times, &vol);

    if (err)
      goto FREE_RETURN;

    Results[i][sUp] = forward[i] * exp(fwdSmileVisuNstd * vol *
                                       sqrt(cpd->pd_leg->cpn[i].fx_val_time));
    Results[i][sATM] = forward[i];
    Results[i][sDown] = forward[i] * exp(-fwdSmileVisuNstd * vol *
                                         sqrt(cpd->pd_leg->cpn[i].fx_val_time));

    PdDfDom[i] =
        swp_f_df(und->today, cpd->pd_leg->cpn[i].pay_date, und->dom_yc) / dfDom;
    FwdDfDom[i] =
        swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->dom_yc) /
        dfDom;
    FwdDfFor[i] =
        swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->for_yc) /
        dfFor;
  }

  /* ==================================
  Value Fx Options for the 3 strikes
  ================================== */
  err = OTCdfPhi(und->lda_dom, cpd->call[otc].ex_time, und->sigma_time_rates,
                 und->sigma_dom, und->sigma_n_rates, &phiDom);

  if (err)
    goto FREE_RETURN;

  err = OTCdfPhi(und->lda_for, cpd->call[otc].ex_time, und->sigma_time_rates,
                 und->sigma_for, und->sigma_n_rates, &phiFor);

  if (err)
    goto FREE_RETURN;

  for (realisation = 0; realisation < nSimul; realisation++) {
    for (i = pd_otc_start; i < cpd->pd_leg->num_cpn; i++) {
      cpn = cpd->pd_leg->cpn + i;

      /* ==================================
      Forward @ Tval
      ================================== */
      err = OTCdf(matrix[realisation][0], und->lda_dom, cpd->call[otc].ex_time,
                  cpn->fx_val_time, tstarTime, phiDom, &df);

      if (err)
        goto FREE_RETURN;

      dfDom = FwdDfDom[i] * df;

      err = OTCdf(matrix[realisation][1], und->lda_for, cpd->call[otc].ex_time,
                  cpn->fx_val_time, tstarTime, phiFor, &df);

      if (err)
        goto FREE_RETURN;

      dfFor = FwdDfFor[i] * df;

      localForward = FxAdj[i] * matrix[realisation][2] * dfFor / dfDom;

      /* ==================================
      Domestic Df @ Tpay
      ================================== */
      err = OTCdf(matrix[realisation][0], und->lda_dom, cpd->call[otc].ex_time,
                  cpn->pay_time, tstarTime, phiDom, &df);

      if (err)
        goto FREE_RETURN;

      dfDom = PdDfDom[i] * df;

      /* ==================================
      UP
      ================================== */
      err = srt_f_optsarbvol(localForward, Results[i][sUp],
                             cpn->fx_fix_time - cpd->call[otc].ex_time,
                             pd_fwdsmile_vol[i], pd_fwdsmile_alpha[i],
                             pd_fwdsmile_beta[i], pd_fwdsmile_rho[i],
                             SRT_BETAVOL, SRT_LOGNORMAL, &smile_std);

      if (err)
        goto FREE_RETURN;

      smile_std *= sqrt(cpn->fx_fix_time - cpd->call[otc].ex_time);
      smile_half_std = 0.5 * smile_std;

      if (localForward < 0.0 && type >= 3) {
        Results[i][vUp] +=
            dfDom * OPT_VAL_MACRO(type - 2, localForward, Results[i][sUp],
                                  smile_std, smile_half_std);
      } else {
        Results[i][vUp] +=
            dfDom * OPT_VAL_MACRO(type, localForward, Results[i][sUp],
                                  smile_std, smile_half_std);
      }

      /* ==================================
      ATM
      ================================== */
      err = srt_f_optsarbvol(localForward, Results[i][sATM],
                             cpn->fx_fix_time - cpd->call[otc].ex_time,
                             pd_fwdsmile_vol[i], pd_fwdsmile_alpha[i],
                             pd_fwdsmile_beta[i], pd_fwdsmile_rho[i],
                             SRT_BETAVOL, SRT_LOGNORMAL, &smile_std);

      if (err)
        goto FREE_RETURN;

      smile_std *= sqrt(cpn->fx_fix_time - cpd->call[otc].ex_time);
      smile_half_std = 0.5 * smile_std;

      if (localForward < 0.0 && type >= 3) {
        Results[i][vATM] +=
            dfDom * OPT_VAL_MACRO(type - 2, localForward, Results[i][sATM],
                                  smile_std, smile_half_std);
      } else {
        Results[i][vATM] +=
            dfDom * OPT_VAL_MACRO(type, localForward, Results[i][sATM],
                                  smile_std, smile_half_std);
      }

      /* ==================================
      Down
      ================================== */
      err = srt_f_optsarbvol(localForward, Results[i][sDown],
                             cpn->fx_fix_time - cpd->call[otc].ex_time,
                             pd_fwdsmile_vol[i], pd_fwdsmile_alpha[i],
                             pd_fwdsmile_beta[i], pd_fwdsmile_rho[i],
                             SRT_BETAVOL, SRT_LOGNORMAL, &smile_std);

      if (err)
        goto FREE_RETURN;

      smile_std *= sqrt(cpn->fx_fix_time - cpd->call[otc].ex_time);
      smile_half_std = 0.5 * smile_std;

      if (localForward < 0.0 && type >= 3) {
        Results[i][vDown] +=
            dfDom * OPT_VAL_MACRO(type - 2, localForward, Results[i][sDown],
                                  smile_std, smile_half_std);
      } else {
        Results[i][vDown] +=
            dfDom * OPT_VAL_MACRO(type, localForward, Results[i][sDown],
                                  smile_std, smile_half_std);
      }
    }
  }

  df = swp_f_df(und->today, cpd->call[otc].ex_date, und->dom_yc);

  for (i = pd_otc_start; i < cpd->pd_leg->num_cpn; i++) {
    dfDom = swp_f_df(und->today, cpd->pd_leg->cpn[i].pay_date, und->dom_yc);

    for (j = 0; j < 3; j++) {
      Results[i][2 * j + 1] *= df / nSimul;

      err = srt_f_optimpvol(Results[i][2 * j + 1], forward[i],
                            Results[i][2 * j], cpd->pd_leg->cpn[i].fx_fix_time,
                            dfDom, SRT_CALL, SRT_LOGNORMAL, &vol);

      if (err)
        Results[i][2 * j + 1] = 0.0;
      else
        Results[i][2 * j + 1] = vol;
    }
  }

FREE_RETURN:

  if (PdDfDom)
    free(PdDfDom);
  if (FwdDfDom)
    free(FwdDfDom);
  if (FwdDfFor)
    free(FwdDfFor);
  if (forward)
    free(forward);

  return err;
}

Err OTCPDpayoffFeeCombinedSmileAndFlat(
    CPD_STR cpd, CPD_UND und, otcpd_params *OTC_params,
    otcpd_precalc *OTC_precalc, SMILE_VOL_MARKET smile_mkt,
    SMILE_PARAMETERS smile_params, double pd_not, double **matrix,
    double tstarTime, int tstarDate, double **Results) {

  int i, realisation;

  /* For the payoff */
  double opt_string;
  double opt_string_flat;

  int str_idx, type;
  double *strikeCap = NULL, *strikeFloor = NULL;
  double forward, forward_flat;
  double cap, cap_flat;
  double floor, floor_flat;
  PD_EXO_CPN cpn;
  double *fundDf = NULL, *PdDfDom = NULL, *FwdDfDom = NULL, *FwdDfFor = NULL;
  double discount, df, dfDom, dfFor, dfFund, dfFirstFunding, dfNotExchgePd,
      dfNotExchgeShortPd, dfNotExchgeFund, dfNotExchgeShortFund;

  /* For the call */
  double *fund_leg = NULL, *pd_leg = NULL;
  double *fund_legShort = NULL, *pd_legShort = NULL;

  double *fund_leg_flat = NULL, *pd_leg_flat = NULL;
  double *fund_legShort_flat = NULL, *pd_legShort_flat = NULL;

  int pd_otc_start,
      pd_otc_short_end; /*	index of the first coupon to be called */
  int fund_otc_start, fund_otc_short_end;

  /* For the fee adjustment */
  double callRealisation;
  double fund_fee, pd_fee;
  double fund_feeShort, pd_feeShort;

  double fund_fee_flat, pd_fee_flat;
  double fund_feeShort_flat, pd_feeShort_flat;

  double coupon;

  /* For the Df */
  double phiDom, phiFor;

  /* For foreign funding */
  char *fund_yc;
  double fund_mult;

  /* For the dollar notional exchange */
  double dfNotEx[4];

  /* Realised discount factors*/
  double *PDdfReal = NULL, *PDFxDomdfReal = NULL, *PDFxFordfReal = NULL;
  double *NewCoupons = NULL;

  /* results */
  double FundingValue = 0.0;
  double PdValue = 0.0;
  double FundingShortValue = 0.0;
  double PdShortValue = 0.0;

  double FundingValue_flat = 0.0;
  double PdValue_flat = 0.0;
  double FundingShortValue_flat = 0.0;
  double PdShortValue_flat = 0.0;

  int Pay_Rec_Mult = 0;

  int Calc_flat = OTC_params->OTCsmileAndFlat;

  SMILE_PARAMETERS smile_params_flat = NULL;

  double OptionValue, OptionValue_flat;

  Err err = NULL;
  /* =============================
  memory allocation
  ============================= */
  fund_leg = calloc(OTC_params->COPULAnSimul, sizeof(double));
  pd_leg = calloc(OTC_params->COPULAnSimul, sizeof(double));
  fund_legShort = calloc(OTC_params->COPULAnSimul, sizeof(double));
  pd_legShort = calloc(OTC_params->COPULAnSimul, sizeof(double));

  fund_leg_flat = calloc(OTC_params->COPULAnSimul, sizeof(double));
  pd_leg_flat = calloc(OTC_params->COPULAnSimul, sizeof(double));
  fund_legShort_flat = calloc(OTC_params->COPULAnSimul, sizeof(double));
  pd_legShort_flat = calloc(OTC_params->COPULAnSimul, sizeof(double));

  fundDf = calloc(cpd->fund_leg->num_cpn, sizeof(double));
  PdDfDom = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  FwdDfDom = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  FwdDfFor = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  strikeCap = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  strikeFloor = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));

  PDdfReal = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  PDFxDomdfReal = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  PDFxFordfReal = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));

  smile_params_flat = calloc(1, sizeof(smile_parameters));

  if (!fund_leg || !pd_leg || !fund_legShort || !pd_legShort || !fundDf ||
      !PdDfDom || !FwdDfDom || !FwdDfFor || !strikeCap || !strikeFloor ||
      !PDdfReal || !PDFxDomdfReal || !PDFxFordfReal) {
    err = "OTCcaller: memory allocation error";
    goto FREE_RETURN;
  }

  /* =============================
  Domestic/foreign
  ============================= */
  if (cpd->fund_leg->dom_for == 0) {
    fund_yc = (char *)(und->dom_yc);
    fund_mult = 1.0;
  } else {
    fund_yc = (char *)(und->for_yc);
    fund_mult = und->spot_fx;
  }

  /* =============================
  Precalculations
  ============================= */

  pd_otc_start = cpd->call[OTC_params->OTC].pd_idx;
  fund_otc_start = cpd->call[OTC_params->OTC].fund_idx;

  if (OTC_params->OTC < cpd->num_calls - 1) {
    pd_otc_short_end = cpd->call[OTC_params->OTC + 1].pd_idx;
    fund_otc_short_end = cpd->call[OTC_params->OTC + 1].fund_idx;
  }

  dfDom = swp_f_df(und->today, cpd->call[OTC_params->OTC].ex_date, und->dom_yc);
  dfFor = swp_f_df(und->today, cpd->call[OTC_params->OTC].ex_date, und->for_yc);
  dfFund = swp_f_df(und->today, cpd->call[OTC_params->OTC].ex_date, fund_yc);

  dfFirstFunding =
      swp_f_df(und->today, cpd->fund_leg->cpn[fund_otc_start].start_date,
               fund_yc) /
      dfFund;
  dfNotExchgePd =
      swp_f_df(und->today, cpd->call[OTC_params->OTC].set_date, und->dom_yc) /
      dfDom;
  dfNotExchgeFund =
      swp_f_df(und->today, cpd->call[OTC_params->OTC].set_date, fund_yc) /
      dfFund;

  if (OTC_params->OTC < cpd->num_calls - 1) {
    dfNotExchgeShortPd =
        swp_f_df(und->today, cpd->call[OTC_params->OTC + 1].set_date,
                 und->dom_yc) /
        dfDom;
    dfNotExchgeShortFund =
        swp_f_df(und->today, cpd->call[OTC_params->OTC + 1].set_date, fund_yc) /
        dfFund;
  }

  for (i = fund_otc_start; i < cpd->fund_leg->num_cpn; i++)
    fundDf[i] =
        swp_f_df(und->today, cpd->fund_leg->cpn[i].pay_date, fund_yc) / dfFund;

  for (i = pd_otc_start; i < cpd->pd_leg->num_cpn; i++) {
    PdDfDom[i] =
        swp_f_df(und->today, cpd->pd_leg->cpn[i].pay_date, und->dom_yc) / dfDom;
    FwdDfDom[i] =
        swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->dom_yc) /
        dfDom;
    FwdDfFor[i] =
        swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->for_yc) /
        dfFor;
    if (fabs(cpd->pd_leg->cpn[i].beta) > 1.0e-16) {
      strikeCap[i] = (cpd->pd_leg->cpn[i].cap - cpd->pd_leg->cpn[i].alpha) /
                     cpd->pd_leg->cpn[i].beta;
      strikeFloor[i] = (cpd->pd_leg->cpn[i].floor - cpd->pd_leg->cpn[i].alpha) /
                       cpd->pd_leg->cpn[i].beta;
    }
  }

  PdDfDom[cpd->pd_leg->num_cpn] =
      swp_f_df(und->today, cpd->pd_leg->not_ref.pay_date, und->dom_yc) / dfDom;
  FwdDfDom[cpd->pd_leg->num_cpn] =
      swp_f_df(und->today, cpd->pd_leg->not_ref.fx_val_date, und->dom_yc) /
      dfDom;
  FwdDfFor[cpd->pd_leg->num_cpn] =
      swp_f_df(und->today, cpd->pd_leg->not_ref.fx_val_date, und->for_yc) /
      dfFor;
  if (fabs(cpd->pd_leg->not_ref.beta) > 1.0e-16) {
    strikeCap[cpd->pd_leg->num_cpn] =
        (cpd->pd_leg->not_ref.cap - cpd->pd_leg->not_ref.alpha) /
        cpd->pd_leg->not_ref.beta;
    strikeFloor[cpd->pd_leg->num_cpn] =
        (cpd->pd_leg->not_ref.floor - cpd->pd_leg->not_ref.alpha) /
        cpd->pd_leg->not_ref.beta;
  }

  /* ==============
  Strike and DF for the Not Exch
  =============== */
  for (i = 0; i < 2; i++) {
    if (OTC_params->OTC + i < cpd->num_calls) {
      if (cpd->call[OTC_params->OTC + i].use_opt_str) {
        dfNotEx[i] =
            (swp_f_df(und->today, cpd->call[OTC_params->OTC + i].fx_val_date,
                      und->for_yc) /
             dfFor) /
            (swp_f_df(und->today, cpd->call[OTC_params->OTC + i].fx_val_date,
                      und->dom_yc) /
             dfDom);
      }
    }
  }

  if (OTC_params->SpeedFunding) {
    NewCoupons = calloc(cpd->pd_leg->num_cpn, sizeof(double));

    err = OTCgetChangedPayoff(cpd, und, NewCoupons);

    if (err)
      goto FREE_RETURN;
  }

  if (!cpd->call[OTC_params->OTC].pay_rec)
    Pay_Rec_Mult = 1;
  else
    Pay_Rec_Mult = -1;

  /* ==========================================================
  =============================================================
          FEE
  =============================================================
  ========================================================== */
  fund_fee = 0.0;
  fund_feeShort = 0.0;
  pd_fee = 0.0;
  pd_feeShort = 0.0;

  fund_fee_flat = 0.0;
  fund_feeShort_flat = 0.0;
  pd_fee_flat = 0.0;
  pd_feeShort_flat = 0.0;

  /* ==================================
  0)	Notional exchange at Tset and last coupon date
  ================================== */
  /* Initial Notional for funding */
  fund_fee -=
      swp_f_df(und->today, cpd->call[OTC_params->OTC].set_date, fund_yc) *
      cpd->fund_leg->notional;
  if (OTC_params->OTC < cpd->num_calls - 1) {
    fund_feeShort = fund_fee;

    /*	Final Notional For Short option*/
    fund_feeShort +=
        swp_f_df(und->today, cpd->call[OTC_params->OTC + 1].set_date, fund_yc) *
        cpd->fund_leg->notional;
  }

  /* Initial Notional for Exo leg */
  for (i = 0; i < 2; i++) {
    if (OTC_params->OTC + i < cpd->num_calls) {
      dfDom = swp_f_df(und->today, cpd->call[OTC_params->OTC + i].set_date,
                       und->dom_yc);

      if (cpd->call[OTC_params->OTC + i].use_opt_str) {

        /* ==================================
        SABR parameters @ Tval
        ================================== */
        err = OTCgetSMILEparams(
            cpd, und, cpd->call[OTC_params->OTC + i].fx_fix_time,
            cpd->call[OTC_params->OTC + i].fx_val_time,
            cpd->call[OTC_params->OTC + i].fx_val_date, smile_mkt, smile_params,
            Calc_flat ? 1 : OTC_params->OTCsmileFee);

        if (err)
          goto FREE_RETURN;

        if (Calc_flat) // As sigma is log  , only the ABR params have to be set
                       // to the log case.
        {
          err = OTCgetSMILEparams(
              cpd, und, cpd->call[OTC_params->OTC + i].fx_fix_time,
              cpd->call[OTC_params->OTC + i].fx_val_time,
              cpd->call[OTC_params->OTC + i].fx_val_date, smile_mkt,
              smile_params_flat, Calc_flat ? 1 : OTC_params->OTCsmileFee);

          if (err)
            goto FREE_RETURN;

          smile_params_flat->alpha = 0.0;
          smile_params_flat->beta = 1.0;
          smile_params_flat->rho = 0.0;
          smile_params_flat->smile_spec_type = 0;
          smile_params_flat->sigma = smile_params->sigma;
        }

        /* ==================================
        Forward @ Tval
        ================================== */
        forward =
            OTC_precalc->NotFeeFxAdj[OTC_params->OTC + i] * und->spot_fx *
            swp_f_df(und->today, cpd->call[OTC_params->OTC + i].fx_val_date,
                     und->for_yc) /
            swp_f_df(und->today, cpd->call[OTC_params->OTC + i].fx_val_date,
                     und->dom_yc);

        /* ==================================
        SABR parameters @ Tval
        ================================== */
        opt_string = 0.0;
        opt_string_flat = 0.0;

        for (str_idx = 0; str_idx < cpd->call[OTC_params->OTC + i].nstrikes;
             str_idx++) {
          if (fabs(cpd->call[OTC_params->OTC + i].weights[str_idx]) >
              1e-16 * pd_not) {
            if (smile_params->sigma *
                        sqrt(cpd->call[OTC_params->OTC + i].fx_fix_time) >
                    1e-16 &&
                cpd->call[OTC_params->OTC + i].strikes[str_idx] > 1e-16)
              type = 3;
            else
              type = 1;

            err = cpd_vol_get_price(
                type, forward, cpd->call[OTC_params->OTC + i].fx_fix_time,
                cpd->call[OTC_params->OTC + i].strikes[str_idx], smile_params,
                &OptionValue);
            if (err)
              goto FREE_RETURN;

            opt_string +=
                cpd->call[OTC_params->OTC + i].weights[str_idx] * OptionValue;

            if (Calc_flat) {
              err = cpd_vol_get_price(
                  type, forward, cpd->call[OTC_params->OTC + i].fx_fix_time,
                  cpd->call[OTC_params->OTC + i].strikes[str_idx],
                  smile_params_flat, &OptionValue_flat);
              if (err)
                goto FREE_RETURN;

              opt_string_flat +=
                  cpd->call[OTC_params->OTC + i].weights[str_idx] *
                  OptionValue_flat;
            }
          }
        }

        /*	Coupon pv */
        if (!i) {
          pd_fee -=
              dfDom * (cpd->call[OTC_params->OTC + i].wcst +
                       cpd->call[OTC_params->OTC + i].wspot * forward +
                       opt_string + cpd->call[OTC_params->OTC + i].orig_fee);
          pd_feeShort = pd_fee;

          if (Calc_flat) {
            pd_fee_flat -=
                dfDom *
                (cpd->call[OTC_params->OTC + i].wcst +
                 cpd->call[OTC_params->OTC + i].wspot * forward +
                 opt_string_flat + cpd->call[OTC_params->OTC + i].orig_fee);
            pd_feeShort_flat = pd_fee_flat;
          }
        } else {
          pd_feeShort +=
              dfDom * (cpd->call[OTC_params->OTC + i].wcst +
                       cpd->call[OTC_params->OTC + i].wspot * forward +
                       opt_string + cpd->call[OTC_params->OTC + i].orig_fee);

          if (Calc_flat) {
            pd_feeShort_flat +=
                dfDom *
                (cpd->call[OTC_params->OTC + i].wcst +
                 cpd->call[OTC_params->OTC + i].wspot * forward +
                 opt_string_flat + cpd->call[OTC_params->OTC + i].orig_fee);
          }
        }
      } else {
        if (!i) {
          pd_fee -= dfDom * (pd_not + cpd->call[OTC_params->OTC + i].orig_fee);
          pd_feeShort = pd_fee;
          if (Calc_flat)
            pd_fee_flat = pd_fee;
          if (Calc_flat)
            pd_feeShort_flat = pd_fee;
        } else {
          pd_feeShort +=
              dfDom * (pd_not + cpd->call[OTC_params->OTC + i].orig_fee);
          if (Calc_flat)
            pd_feeShort_flat +=
                dfDom * (pd_not + cpd->call[OTC_params->OTC + i].orig_fee);
        }
      }
    }
  }

  /*	Final Notional For the Long Option*/
  cpn = &(cpd->pd_leg->not_ref);
  /*	Discount */
  dfDom = swp_f_df(und->today, cpn->pay_date, und->dom_yc);

  /* ==================================
  Forward @ Tval
  ================================== */
  forward = OTC_precalc->NotFeeFxAdj[cpd->num_calls] * und->spot_fx *
            swp_f_df(und->today, cpn->fx_val_date, und->for_yc) /
            swp_f_df(und->today, cpn->fx_val_date, und->dom_yc);

  /* ==================================
  SABR parameters @ Tval
  ================================== */
  err = OTCgetSMILEparams(cpd, und, cpn->fx_fix_time, cpn->fx_val_time,
                          cpn->fx_val_date, smile_mkt, smile_params,
                          Calc_flat ? 1 : OTC_params->OTCsmileFee);

  if (err)
    goto FREE_RETURN;

  if (Calc_flat) {
    err = OTCgetSMILEparams(cpd, und, cpn->fx_fix_time, cpn->fx_val_time,
                            cpn->fx_val_date, smile_mkt, smile_params_flat,
                            Calc_flat ? 1 : OTC_params->OTCsmileFee);

    if (err)
      goto FREE_RETURN;

    smile_params_flat->alpha = 0.0;
    smile_params_flat->beta = 1.0;
    smile_params_flat->rho = 0.0;
    smile_params_flat->smile_spec_type = 0;
  }

  if (!cpn->use_opt_str) {
    if (fabs(cpn->beta) > 1.0e-16) {
      /* ==================================
      Floor
      ================================== */
      err = OTCgetOptValue(cpn->floored, 0, cpn->floor,
                           strikeFloor[cpd->pd_leg->num_cpn], cpn->beta,
                           forward, cpn->fx_fix_time, smile_params, 0, &floor);

      if (err)
        goto FREE_RETURN;

      if (Calc_flat) {

        err = OTCgetOptValue(cpn->floored, 0, cpn->floor,
                             strikeFloor[cpd->pd_leg->num_cpn], cpn->beta,
                             forward, cpn->fx_fix_time, smile_params_flat, 0,
                             &floor_flat);

        if (err)
          goto FREE_RETURN;
      }
      /* ==================================
      Cap
      ================================== */
      err = OTCgetOptValue(cpn->capped, 1, cpn->cap,
                           strikeCap[cpd->pd_leg->num_cpn], cpn->beta, forward,
                           cpn->fx_fix_time, smile_params, 0, &cap);

      if (err)
        goto FREE_RETURN;

      if (Calc_flat) {
        err =
            OTCgetOptValue(cpn->capped, 1, cpn->cap,
                           strikeCap[cpd->pd_leg->num_cpn], cpn->beta, forward,
                           cpn->fx_fix_time, smile_params_flat, 0, &cap_flat);

        if (err)
          goto FREE_RETURN;
      }
    } else {
      forward = floor = cap = floor_flat = cap_flat = 0.0;
    }

    /*	Coupon pv */
    pd_fee += dfDom * (cpn->alpha + cpn->beta * forward +
                       fabs(cpn->beta) * (floor - cap));
    if (Calc_flat)
      pd_fee_flat += dfDom * (cpn->alpha + cpn->beta * forward +
                              fabs(cpn->beta) * (floor_flat - cap_flat));
  } else {
    opt_string = 0.0;
    opt_string_flat = 0.0;

    for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
      if (fabs(cpn->weights[str_idx]) > 1e-16 * pd_not) {
        if (smile_params->sigma * sqrt(cpn->fx_fix_time) > 1e-16 &&
            cpn->strikes[str_idx] > 1e-16)
          type = 3;
        else
          type = 1;

        err = cpd_vol_get_price(type, forward, cpn->fx_fix_time,
                                cpn->strikes[str_idx], smile_params,
                                &OptionValue);
        if (err)
          goto FREE_RETURN;

        opt_string += cpn->weights[str_idx] * OptionValue;

        if (Calc_flat) {
          err = cpd_vol_get_price(type, forward, cpn->fx_fix_time,
                                  cpn->strikes[str_idx], smile_params_flat,
                                  &OptionValue_flat);
          if (err)
            goto FREE_RETURN;

          opt_string_flat += cpn->weights[str_idx] * OptionValue_flat;
        }
      }
    }

    /*	Coupon pv */
    pd_fee += dfDom * (cpn->wcst + cpn->wspot * forward + opt_string);
    if (Calc_flat)
      pd_fee_flat +=
          dfDom * (cpn->wcst + cpn->wspot * forward + opt_string_flat);
  }

  /* ==================================
  1)	Funding Fee
  ================================== */
  /* B(t  ,Ts) */

  fund_fee += swp_f_df(und->today,
                       cpd->fund_leg->cpn[fund_otc_start].start_date, fund_yc) *
              cpd->fund_leg->notional;

  if (OTC_params->OTC < cpd->num_calls - 1)
    fund_feeShort +=
        swp_f_df(und->today, cpd->fund_leg->cpn[fund_otc_start].start_date,
                 fund_yc) *
        cpd->fund_leg->notional;

  for (i = fund_otc_start; i < cpd->fund_leg->num_cpn; i++) {
    /* ==================================
    coupon: spread + margin
    ================================== */
    coupon = swp_f_df(und->today, cpd->fund_leg->cpn[i].pay_date, fund_yc);

    fund_fee += coupon * cpd->fund_leg->cpn[i].cpn;

    if (OTC_params->OTC < cpd->num_calls - 1) {
      if (i >= fund_otc_start && i < fund_otc_short_end) {
        fund_feeShort += coupon * cpd->fund_leg->cpn[i].cpn;
      }

      if (i == fund_otc_short_end - 1)
        fund_feeShort -= coupon * cpd->fund_leg->notional;
    }
  }

  /*	In order to get the funding PV in domestic units */
  fund_fee *= fund_mult;
  fund_feeShort *= fund_mult;

  /* ==================================
  2)	Power dual Fee
  ================================== */
  for (i = pd_otc_start; i < cpd->pd_leg->num_cpn; i++) {
    cpn = cpd->pd_leg->cpn + i;

    /* ==================================
    Forward @ Tval with Convexity adjustment in the 3F
    ================================== */
    forward = OTC_precalc->FeeFxAdj[i] * und->spot_fx *
              swp_f_df(und->today, cpn->fx_val_date, und->for_yc) /
              swp_f_df(und->today, cpn->fx_val_date, und->dom_yc);

    /* ==================================
    SABR parameters @ Tval
    ================================== */
    err = OTCgetSMILEparams(cpd, und, cpn->fx_fix_time, cpn->fx_val_time,
                            cpn->fx_val_date, smile_mkt, smile_params,
                            Calc_flat ? 1 : OTC_params->OTCsmileFee);

    if (err)
      goto FREE_RETURN;

    smile_params_flat->sigma = smile_params->sigma;

    /* ==================================
    domestic DF @ Tpay
    ================================== */
    dfDom = swp_f_df(und->today, cpn->pay_date, und->dom_yc);

    if (!cpn->use_opt_str) {
      if (fabs(cpn->beta) > 1.0e-16) {
        /* ==================================
        Floor
        ================================== */
        err = OTCgetOptValue(cpn->floored, 0, cpn->floor, strikeFloor[i],
                             cpn->beta, forward, cpn->fx_fix_time, smile_params,
                             0, &floor);

        if (err)
          goto FREE_RETURN;

        if (Calc_flat) {
          err = OTCgetOptValue(cpn->floored, 0, cpn->floor, strikeFloor[i],
                               cpn->beta, forward, cpn->fx_fix_time,
                               smile_params_flat, 0, &floor_flat);

          if (err)
            goto FREE_RETURN;
        }
        /* ==================================
        Cap
        ================================== */
        err = OTCgetOptValue(cpn->capped, 1, cpn->cap, strikeCap[i], cpn->beta,
                             forward, cpn->fx_fix_time, smile_params, 0, &cap);

        if (err)
          goto FREE_RETURN;

        if (Calc_flat) {
          smile_params_flat->sigma = smile_params->sigma;

          err = OTCgetOptValue(cpn->capped, 1, cpn->cap, strikeCap[i],
                               cpn->beta, forward, cpn->fx_fix_time,
                               smile_params_flat, 0, &cap_flat);

          if (err)
            goto FREE_RETURN;
        }
      } else {
        forward = floor = floor_flat = cap = cap_flat = 0.0;
      }

      pd_fee += dfDom * (cpn->alpha + cpn->beta * forward +
                         fabs(cpn->beta) * (floor - cap));
      if (Calc_flat)
        pd_fee_flat += dfDom * (cpn->alpha + cpn->beta * forward +
                                fabs(cpn->beta) * (floor_flat - cap_flat));

      if (OTC_params->OTC < cpd->num_calls - 1) {
        if (i >= pd_otc_start && i < pd_otc_short_end) {
          pd_feeShort += dfDom * (cpn->alpha + cpn->beta * forward +
                                  fabs(cpn->beta) * (floor - cap));
          if (Calc_flat)
            pd_feeShort_flat +=
                dfDom * (cpn->alpha + cpn->beta * forward +
                         fabs(cpn->beta) * (floor_flat - cap_flat));
        }
      }
    } else {
      opt_string = 0.0;
      opt_string_flat = 0.0;

      for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
        if (fabs(cpn->weights[str_idx]) > 1e-16 * pd_not) {
          if (smile_params->sigma * sqrt(cpn->fx_fix_time) > 1e-16 &&
              cpn->strikes[str_idx] > 1e-16)
            type = 3;
          else
            type = 1;

          err = cpd_vol_get_price(type, forward, cpn->fx_fix_time,
                                  cpn->strikes[str_idx], smile_params,
                                  &OptionValue);
          if (err)
            goto FREE_RETURN;

          opt_string += cpn->weights[str_idx] * OptionValue;

          if (Calc_flat) {
            err = cpd_vol_get_price(type, forward, cpn->fx_fix_time,
                                    cpn->strikes[str_idx], smile_params_flat,
                                    &OptionValue_flat);
            if (err)
              goto FREE_RETURN;

            opt_string_flat += cpn->weights[str_idx] * OptionValue_flat;
          }
        }
      }

      /*	Coupon pv */
      pd_fee += dfDom * (cpn->wcst + cpn->wspot * forward + opt_string);
      if (Calc_flat)
        pd_fee_flat +=
            dfDom * (cpn->wcst + cpn->wspot * forward + opt_string_flat);

      if (OTC_params->OTC < cpd->num_calls - 1) {
        if (i >= pd_otc_start && i < pd_otc_short_end) {
          pd_feeShort +=
              dfDom * (cpn->wcst + cpn->wspot * forward + opt_string);
          if (Calc_flat)
            pd_feeShort_flat +=
                dfDom * (cpn->wcst + cpn->wspot * forward + opt_string_flat);
        }
      }
    }
  }

  if (OTC_params->OTCcalculate[OTC_params->OTC]) {
    smile_params->Use_BetaQuick_in_MC = 1;
    smile_params_flat->Use_BetaQuick_in_MC = 1;
    /* ==========================================================
    =============================================================
            PAYOFF
    =============================================================
    ========================================================== */
    err = OTCdfPhi(und->lda_dom, cpd->call[OTC_params->OTC].ex_time,
                   und->sigma_time_rates, und->sigma_dom, und->sigma_n_rates,
                   &phiDom);

    if (err)
      goto FREE_RETURN;

    err = OTCdfPhi(und->lda_for, cpd->call[OTC_params->OTC].ex_time,
                   und->sigma_time_rates, und->sigma_for, und->sigma_n_rates,
                   &phiFor);

    if (err)
      goto FREE_RETURN;

    for (realisation = 0; realisation < OTC_params->COPULAnSimul;
         realisation++) {
      /* ==================================
      A) Calculates all the Dom DF at Pwd dates
      ================================== */
      for (i = pd_otc_start; i < cpd->pd_leg->num_cpn; i++) {
        cpn = cpd->pd_leg->cpn + i;

        df = exp(OTC_precalc->df_const_pd_pay[OTC_params->OTC][i] +
                 OTC_precalc->df_lin_pd_pay[OTC_params->OTC][i] *
                     matrix[realisation][0]);

        PDdfReal[i] = PdDfDom[i] * df;

        if (cpn->pay_date != cpn->fx_val_date) {
          df = exp(OTC_precalc->df_const_dom_fx_val[OTC_params->OTC][i] +
                   OTC_precalc->df_lin_dom_fx_val[OTC_params->OTC][i] *
                       matrix[realisation][0]);
        }

        PDFxDomdfReal[i] = FwdDfDom[i] * df;

        df = exp(OTC_precalc->df_const_for_fx_val[OTC_params->OTC][i] +
                 OTC_precalc->df_lin_for_fx_val[OTC_params->OTC][i] *
                     matrix[realisation][1]);
        PDFxFordfReal[i] = FwdDfFor[i] * df;
      }

      /* ==================================
      0)	Notional exchange at Tset and last coupon date
      ================================== */
      for (i = 0; i < 2; i++) {
        if (OTC_params->OTC + i < cpd->num_calls) {
          df = exp(OTC_precalc->df_const_pd_pay[OTC_params->OTC]
                                               [cpd->pd_leg->num_cpn + i + 1] +
                   OTC_precalc->df_lin_pd_pay[OTC_params->OTC]
                                             [cpd->pd_leg->num_cpn + i + 1] *
                       matrix[realisation][0]);

          if (!i)
            dfDom = df * dfNotExchgePd;
          else
            dfDom = df * dfNotExchgeShortPd;

          err = cpd_vol_get_smile_params_from_otc(
              smile_mkt->smile_spec_type + 1, OTC_params->OTC,
              cpd->pd_leg->num_cpn + i + 1, OTC_precalc, smile_params);
          if (err)
            goto FREE_RETURN;

          smile_params_flat->sigma =
              OTC_precalc
                  ->pd_fwd3F_vol[OTC_params->OTC][cpd->pd_leg->num_cpn + i + 1];

          if (cpd->call[OTC_params->OTC + i].use_opt_str) {
            /* ==================================
            Forward @ Tval
            ================================== */
            forward = OTC_precalc->NotFxAdj[OTC_params->OTC][i] *
                      matrix[realisation][2] * dfNotEx[i];
            if (Calc_flat)
              forward_flat = OTC_precalc->NotFxAdj[OTC_params->OTC][i] *
                             matrix[realisation][3] * dfNotEx[i];

            df = exp(
                OTC_precalc->df_const_dom_fx_val[OTC_params->OTC]
                                                [cpd->pd_leg->num_cpn + i + 1] +
                OTC_precalc->df_lin_dom_fx_val[OTC_params->OTC]
                                              [cpd->pd_leg->num_cpn + i + 1] *
                    matrix[realisation][0]);

            forward /= df;
            if (Calc_flat)
              forward_flat /= df;

            df = exp(
                OTC_precalc->df_const_for_fx_val[OTC_params->OTC]
                                                [cpd->pd_leg->num_cpn + i + 1] +
                OTC_precalc->df_lin_for_fx_val[OTC_params->OTC]
                                              [cpd->pd_leg->num_cpn + i + 1] *
                    matrix[realisation][1]);

            forward *= df;
            if (Calc_flat)
              forward_flat *= df;

            opt_string = 0.0;
            opt_string_flat = 0.0;
            for (str_idx = 0; str_idx < cpd->call[OTC_params->OTC + i].nstrikes;
                 str_idx++) {
              if (fabs(cpd->call[OTC_params->OTC + i].weights[str_idx]) >
                  1e-16 * pd_not) {
                if (OTC_precalc->pd_fwdsmile_vol[OTC_params->OTC]
                                                [cpd->pd_leg->num_cpn + i + 1] *
                            sqrt(cpd->call[OTC_params->OTC + i].fx_fix_time -
                                 cpd->call[OTC_params->OTC].ex_time) >
                        1e-16 &&
                    cpd->call[OTC_params->OTC + i].strikes[str_idx] > 1e-16)
                  type = 3;
                else
                  type = 1;

                err = cpd_vol_get_price(
                    type, forward,
                    cpd->call[OTC_params->OTC + i].fx_fix_time -
                        cpd->call[OTC_params->OTC].ex_time,
                    cpd->call[OTC_params->OTC + i].strikes[str_idx],
                    smile_params, &OptionValue);
                if (err)
                  goto FREE_RETURN;

                opt_string += cpd->call[OTC_params->OTC + i].weights[str_idx] *
                              OptionValue;

                if (Calc_flat) {
                  err = cpd_vol_get_price(
                      type, forward_flat,
                      cpd->call[OTC_params->OTC + i].fx_fix_time -
                          cpd->call[OTC_params->OTC].ex_time,
                      cpd->call[OTC_params->OTC + i].strikes[str_idx],
                      smile_params_flat, &OptionValue_flat);
                  if (err)
                    goto FREE_RETURN;

                  opt_string_flat +=
                      cpd->call[OTC_params->OTC + i].weights[str_idx] *
                      OptionValue_flat;
                }
              }
            }

            if (!i) {
              pd_leg[realisation] -=
                  dfDom *
                  (cpd->call[OTC_params->OTC + i].wcst +
                   cpd->call[OTC_params->OTC + i].wspot * forward + opt_string +
                   cpd->call[OTC_params->OTC + i].orig_fee);
              pd_legShort[realisation] -=
                  dfDom *
                  (cpd->call[OTC_params->OTC + i].wcst +
                   cpd->call[OTC_params->OTC + i].wspot * forward + opt_string +
                   cpd->call[OTC_params->OTC + i].orig_fee);

              if (Calc_flat) {
                pd_leg_flat[realisation] -=
                    dfDom *
                    (cpd->call[OTC_params->OTC + i].wcst +
                     cpd->call[OTC_params->OTC + i].wspot * forward_flat +
                     opt_string_flat + cpd->call[OTC_params->OTC + i].orig_fee);
                pd_legShort_flat[realisation] -=
                    dfDom *
                    (cpd->call[OTC_params->OTC + i].wcst +
                     cpd->call[OTC_params->OTC + i].wspot * forward_flat +
                     opt_string_flat + cpd->call[OTC_params->OTC + i].orig_fee);
              }
            } else {
              pd_legShort[realisation] +=
                  dfDom *
                  (cpd->call[OTC_params->OTC + i].wcst +
                   cpd->call[OTC_params->OTC + i].wspot * forward + opt_string +
                   cpd->call[OTC_params->OTC + i].orig_fee);

              if (Calc_flat) {
                pd_legShort_flat[realisation] +=
                    dfDom *
                    (cpd->call[OTC_params->OTC + i].wcst +
                     cpd->call[OTC_params->OTC + i].wspot * forward_flat +
                     opt_string_flat + cpd->call[OTC_params->OTC + i].orig_fee);
              }
            }

          } else {
            if (!i) {
              pd_leg[realisation] -=
                  dfDom * (pd_not + cpd->call[OTC_params->OTC + i].orig_fee);
              pd_legShort[realisation] = pd_leg[realisation];
              if (Calc_flat) {
                pd_leg_flat[realisation] -=
                    dfDom * (pd_not + cpd->call[OTC_params->OTC + i].orig_fee);
                pd_legShort_flat[realisation] = pd_leg_flat[realisation];
              }
            } else {
              pd_legShort[realisation] +=
                  dfDom * (pd_not + cpd->call[OTC_params->OTC + i].orig_fee);
              if (Calc_flat)
                pd_legShort_flat[realisation] +=
                    dfDom * (pd_not + cpd->call[OTC_params->OTC + i].orig_fee);
            }
          }
        }
      }

      /*	Final Notional For Long option*/
      cpn = &(cpd->pd_leg->not_ref);

      /*	Discount */
      df = exp(
          OTC_precalc->df_const_pd_pay[OTC_params->OTC][cpd->pd_leg->num_cpn] +
          OTC_precalc->df_lin_pd_pay[OTC_params->OTC][cpd->pd_leg->num_cpn] *
              matrix[realisation][0]);

      discount = PdDfDom[cpd->pd_leg->num_cpn] * df;

      /*	Fwd fx */

      if ((fabs(cpn->beta) > 1.0e-16 && !cpn->use_opt_str) ||
          cpn->use_opt_str) {
        /* ==================================
        Forward @ Tval
        ================================== */
        df = exp(
            OTC_precalc
                ->df_const_dom_fx_val[OTC_params->OTC][cpd->pd_leg->num_cpn] +
            OTC_precalc
                    ->df_lin_dom_fx_val[OTC_params->OTC][cpd->pd_leg->num_cpn] *
                matrix[realisation][0]);

        dfDom = FwdDfDom[cpd->pd_leg->num_cpn] * df;

        df = exp(
            OTC_precalc
                ->df_const_for_fx_val[OTC_params->OTC][cpd->pd_leg->num_cpn] +
            OTC_precalc
                    ->df_lin_for_fx_val[OTC_params->OTC][cpd->pd_leg->num_cpn] *
                matrix[realisation][1]);

        dfFor = FwdDfFor[cpd->pd_leg->num_cpn] * df;

        forward = OTC_precalc->NotFxAdj[OTC_params->OTC][2] *
                  matrix[realisation][2] * dfFor / dfDom;
        if (Calc_flat)
          forward_flat = OTC_precalc->NotFxAdj[OTC_params->OTC][2] *
                         matrix[realisation][3] * dfFor / dfDom;
      }

      err = cpd_vol_get_smile_params_from_otc(
          smile_mkt->smile_spec_type + 1, OTC_params->OTC, cpd->pd_leg->num_cpn,
          OTC_precalc, smile_params);
      if (err)
        goto FREE_RETURN;

      smile_params_flat->sigma =
          OTC_precalc->pd_fwd3F_vol[OTC_params->OTC][cpd->pd_leg->num_cpn];

      if (!cpn->use_opt_str) {
        if (fabs(cpn->beta) > 1.0e-16) {
          /* ==================================
          Floor
          ================================== */
          err = OTCgetOptValue(
              cpn->floored, 0, cpn->floor, strikeFloor[cpd->pd_leg->num_cpn],
              cpn->beta, forward,
              cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
              smile_params, 1, &floor);

          if (err)
            goto FREE_RETURN;

          if (Calc_flat) {
            err = OTCgetOptValue(
                cpn->floored, 0, cpn->floor, strikeFloor[cpd->pd_leg->num_cpn],
                cpn->beta, forward_flat,
                cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
                smile_params_flat, 1, &floor_flat);

            if (err)
              goto FREE_RETURN;
          }
          /* ==================================
          Cap
          ================================== */
          err = OTCgetOptValue(
              cpn->capped, 1, cpn->cap, strikeCap[cpd->pd_leg->num_cpn],
              cpn->beta, forward,
              cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
              smile_params, 1, &cap);

          if (err)
            goto FREE_RETURN;

          if (Calc_flat) {
            err = OTCgetOptValue(
                cpn->capped, 1, cpn->cap, strikeCap[cpd->pd_leg->num_cpn],
                cpn->beta, forward_flat,
                cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
                smile_params_flat, 1, &cap_flat);

            if (err)
              goto FREE_RETURN;
          }
        } else {
          forward = forward_flat = floor = cap = floor_flat = cap_flat = 0.0;
        }

        /*	Coupon pv */
        pd_leg[realisation] += discount * (cpn->alpha + cpn->beta * forward +
                                           fabs(cpn->beta) * (floor - cap));
        if (Calc_flat)
          pd_leg_flat[realisation] +=
              discount * (cpn->alpha + cpn->beta * forward_flat +
                          fabs(cpn->beta) * (floor_flat - cap_flat));
      } else {
        opt_string = 0.0;
        opt_string_flat = 0.0;

        for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
          if (fabs(cpn->weights[str_idx]) > 1e-16 * pd_not) {
            if (OTC_precalc->pd_fwdsmile_vol[OTC_params->OTC]
                                            [cpd->pd_leg->num_cpn] *
                        sqrt(cpn->fx_fix_time -
                             cpd->call[OTC_params->OTC].ex_time) >
                    1e-16 &&
                cpn->strikes[str_idx] > 1e-16)
              type = 3;
            else
              type = 1;

            err = cpd_vol_get_price(
                type, forward,
                cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
                cpn->strikes[str_idx], smile_params, &OptionValue);
            if (err)
              goto FREE_RETURN;

            opt_string += cpn->weights[str_idx] * OptionValue;

            if (Calc_flat) {
              err = cpd_vol_get_price(
                  type, forward_flat,
                  cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
                  cpn->strikes[str_idx], smile_params_flat, &OptionValue_flat);
              if (err)
                goto FREE_RETURN;

              opt_string += cpn->weights[str_idx] * OptionValue_flat;
            }
          }
        }

        /*	Coupon pv */
        pd_leg[realisation] +=
            discount * (cpn->wcst + cpn->wspot * forward + opt_string);
        if (Calc_flat)
          pd_leg_flat[realisation] +=
              discount *
              (cpn->wcst + cpn->wspot * forward_flat + opt_string_flat);
      }

      /* ==================================
      1)	Value funding leg starting at call date
      ================================== */
      if (!OTC_params->SpeedFunding) {
        err = OTCfundingValuation(
            matrix[realisation][0], matrix[realisation][1],
            matrix[realisation][2], cpd, und, OTC_precalc, tstarTime, phiDom,
            phiFor, OTC_params->OTC, dfNotExchgeFund, dfNotExchgeShortFund,
            dfFirstFunding, fundDf, &(fund_leg[realisation]),
            &(fund_legShort[realisation]));
      } else {
        err = OTCNewfundingValuation(
            matrix[realisation][0], matrix[realisation][1],
            matrix[realisation][2], cpd, und, OTC_precalc, tstarTime, phiDom,
            phiFor, OTC_params->OTC, dfNotExchgeFund, dfNotExchgeShortFund,
            dfFirstFunding, NewCoupons, PDdfReal, PDFxFordfReal,
            &(fund_leg[realisation]), &(fund_legShort[realisation]));
      }
      if (err)
        goto FREE_RETURN;

      if (Calc_flat) {
        if (cpd->fund_leg->dom_for == 0) {
          fund_leg_flat[realisation] = fund_leg[realisation];
          fund_legShort_flat[realisation] = fund_legShort[realisation];
        } else {
          if (fabs(matrix[realisation][2]) > 1.0e-10) {
            fund_leg_flat[realisation] = fund_leg[realisation] *
                                         matrix[realisation][3] /
                                         matrix[realisation][2];
            fund_legShort_flat[realisation] = fund_legShort[realisation] *
                                              matrix[realisation][3] /
                                              matrix[realisation][2];
          } else {
            err = "Funding Conversion impossible (realisation == 0.0)";
            goto FREE_RETURN;
          }
        }
      }
      /* ==================================
      2)	Value power dual leg starting at call date
      ================================== */
      for (i = pd_otc_start; i < cpd->pd_leg->num_cpn; i++) {
        cpn = cpd->pd_leg->cpn + i;

        /* ==================================
        Forward @ Tval with Convexity adjustment
        ================================== */
        forward = OTC_precalc->FxAdj[OTC_params->OTC][i] *
                  matrix[realisation][2] * PDFxFordfReal[i] / PDFxDomdfReal[i];
        if (Calc_flat)
          forward_flat = OTC_precalc->FxAdj[OTC_params->OTC][i] *
                         matrix[realisation][3] * PDFxFordfReal[i] /
                         PDFxDomdfReal[i];

        /* ==================================
        domestic DF @ Tpay
        ================================== */
        // the forward vols are of beta type and vol mkt of log type
        err = cpd_vol_get_smile_params_from_otc(smile_mkt->smile_spec_type + 1,
                                                OTC_params->OTC, i, OTC_precalc,
                                                smile_params);
        if (err)
          goto FREE_RETURN;

        smile_params_flat->sigma =
            OTC_precalc->pd_fwd3F_vol[OTC_params->OTC][i];

        if (!cpn->use_opt_str) {
          if (fabs(cpn->beta) > 1.0e-16) {
            /* ==================================
            Floor
            ================================== */
            err = OTCgetOptValue(
                cpn->floored, 0, cpn->floor, strikeFloor[i], cpn->beta, forward,
                cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
                smile_params, 1, &floor);

            if (err)
              goto FREE_RETURN;

            if (Calc_flat) {
              err = OTCgetOptValue(cpn->floored, 0, cpn->floor, strikeFloor[i],
                                   cpn->beta, forward_flat,
                                   cpn->fx_fix_time -
                                       cpd->call[OTC_params->OTC].ex_time,
                                   smile_params_flat, 1, &floor_flat);

              if (err)
                goto FREE_RETURN;
            }
            /* ==================================
            Cap
            ================================== */
            err = OTCgetOptValue(
                cpn->capped, 1, cpn->cap, strikeCap[i], cpn->beta, forward,
                cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
                smile_params, 1, &cap);

            if (err)
              goto FREE_RETURN;

            if (Calc_flat) {
              err = OTCgetOptValue(cpn->capped, 1, cpn->cap, strikeCap[i],
                                   cpn->beta, forward_flat,
                                   cpn->fx_fix_time -
                                       cpd->call[OTC_params->OTC].ex_time,
                                   smile_params_flat, 1, &cap_flat);

              if (err)
                goto FREE_RETURN;
            }
          } else {
            forward = forward_flat = floor = floor_flat = cap = cap_flat = 0.0;
          }

          pd_leg[realisation] +=
              PDdfReal[i] * (cpn->alpha + cpn->beta * forward +
                             fabs(cpn->beta) * (floor - cap));
          if (Calc_flat)
            pd_leg_flat[realisation] +=
                PDdfReal[i] * (cpn->alpha + cpn->beta * forward_flat +
                               fabs(cpn->beta) * (floor_flat - cap_flat));
          if (OTC_params->OTC < cpd->num_calls - 1) {
            if (i >= pd_otc_start && i < pd_otc_short_end) {
              pd_legShort[realisation] +=
                  PDdfReal[i] * (cpn->alpha + cpn->beta * forward +
                                 fabs(cpn->beta) * (floor - cap));
              if (Calc_flat)
                pd_legShort_flat[realisation] +=
                    PDdfReal[i] * (cpn->alpha + cpn->beta * forward_flat +
                                   fabs(cpn->beta) * (floor_flat - cap_flat));
            }
          }
        } else {
          opt_string = 0.0;
          opt_string_flat = 0.0;

          for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
            if (fabs(cpn->weights[str_idx]) > 1e-16 * pd_not) {

              if (OTC_precalc->pd_fwdsmile_vol[OTC_params->OTC][i] *
                          sqrt(cpn->fx_fix_time -
                               cpd->call[OTC_params->OTC].ex_time) >
                      1e-16 &&
                  cpn->strikes[str_idx] > 1e-16)
                type = 3;
              else
                type = 1;

              err = cpd_vol_get_price(
                  type, forward,
                  cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
                  cpn->strikes[str_idx], smile_params, &OptionValue);
              if (err)
                goto FREE_RETURN;

              opt_string += cpn->weights[str_idx] * OptionValue;

              if (Calc_flat) {
                err = cpd_vol_get_price(type, forward_flat,
                                        cpn->fx_fix_time -
                                            cpd->call[OTC_params->OTC].ex_time,
                                        cpn->strikes[str_idx],
                                        smile_params_flat, &OptionValue_flat);
                if (err)
                  goto FREE_RETURN;

                opt_string_flat += cpn->weights[str_idx] * OptionValue_flat;
              }
            }
          }

          /*	Coupon pv */
          pd_leg[realisation] +=
              PDdfReal[i] * (cpn->wcst + cpn->wspot * forward + opt_string);
          if (Calc_flat)
            pd_leg_flat[realisation] +=
                PDdfReal[i] *
                (cpn->wcst + cpn->wspot * forward_flat + opt_string_flat);
          if (OTC_params->OTC < cpd->num_calls - 1) {
            if (i >= pd_otc_start && i < pd_otc_short_end) {
              pd_legShort[realisation] +=
                  PDdfReal[i] * (cpn->wcst + cpn->wspot * forward + opt_string);
              if (Calc_flat)
                pd_legShort_flat[realisation] +=
                    PDdfReal[i] *
                    (cpn->wcst + cpn->wspot * forward_flat + opt_string_flat);
            }
          }
        }
      }

      if (OTC_params->OTC == cpd->num_calls - 1) {
        fund_legShort[realisation] = fund_leg[realisation];
        pd_legShort[realisation] = pd_leg[realisation];

        if (Calc_flat) {
          fund_legShort_flat[realisation] = fund_leg_flat[realisation];
          pd_legShort_flat[realisation] = pd_leg_flat[realisation];
        }
      }
    }

    if (OTC_params->OTC == cpd->num_calls - 1) {
      fund_feeShort = fund_fee;
      pd_feeShort = pd_fee;
      if (Calc_flat)
        pd_feeShort_flat = pd_fee_flat;
    }

    df = swp_f_df(und->today, tstarDate, und->dom_yc);

    if (Calc_flat)
      Results[OTC_params->OTC + 1][0] = tstarTime;

    for (realisation = 0; realisation < OTC_params->COPULAnSimul;
         realisation++) {
      FundingValue += fund_leg[realisation];
      PdValue += pd_leg[realisation];
      FundingShortValue += fund_legShort[realisation];
      PdShortValue += pd_legShort[realisation];

      if (Calc_flat) {
        FundingValue_flat += fund_leg_flat[realisation];
        PdValue_flat += pd_leg_flat[realisation];
        FundingShortValue_flat += fund_legShort_flat[realisation];
        PdShortValue_flat += pd_legShort_flat[realisation];
      }
    }

    FundingValue /= OTC_params->COPULAnSimul;
    PdValue /= OTC_params->COPULAnSimul;
    FundingShortValue /= OTC_params->COPULAnSimul;
    PdShortValue /= OTC_params->COPULAnSimul;

    if (Calc_flat) {
      FundingValue_flat /= OTC_params->COPULAnSimul;
      PdValue_flat /= OTC_params->COPULAnSimul;
      FundingShortValue_flat /= OTC_params->COPULAnSimul;
      PdShortValue_flat /= OTC_params->COPULAnSimul;
    }
    /* ============================
    Get the adjustment at Tex
    ============================ */
    fund_fee /= df;
    pd_fee /= df;
    fund_feeShort /= df;
    pd_feeShort /= df;
    if (Calc_flat) {
      fund_fee_flat = fund_fee;
      pd_fee_flat /= df;
      fund_feeShort_flat = fund_feeShort;
      pd_feeShort_flat /= df;
    }

    if (Calc_flat) {
      Results[OTC_params->OTC + 1][0] = tstarTime;
      Results[OTC_params->OTC + 1][1] = Pay_Rec_Mult * (pd_fee - fund_fee);
      Results[OTC_params->OTC + 1][5] =
          Pay_Rec_Mult * (pd_fee_flat - fund_fee_flat);
    } else {
      Results[0][0] = fund_fee;
      Results[1][0] = pd_fee;
      Results[0][1] = fund_feeShort;
      Results[1][1] = pd_feeShort;

      Results[0][3] = fund_fee - FundingValue;
      Results[1][3] = pd_fee - PdValue;
      Results[0][4] = fund_feeShort - FundingShortValue;
      Results[1][4] = pd_feeShort - PdShortValue;
    }

    fund_fee -= FundingValue;
    fund_feeShort -= FundingShortValue;
    pd_fee -= PdValue;
    pd_feeShort -= PdShortValue;

    if (Calc_flat) {
      fund_fee_flat -= FundingValue_flat;
      fund_feeShort_flat -= FundingShortValue_flat;
      pd_fee_flat -= PdValue_flat;
      pd_feeShort_flat -= PdShortValue_flat;
    }

    if (OTC_params->OTCsmileFee == -1) {
      fund_fee = 0.0;
      fund_feeShort = 0.0;
      pd_fee = 0.0;
      pd_feeShort = 0.0;
      fund_fee_flat = 0.0;
      fund_feeShort_flat = 0.0;
      pd_fee_flat = 0.0;
      pd_feeShort_flat = 0.0;
    }

    if (Calc_flat) {
      Results[OTC_params->OTC + 1][9] = fund_fee;
      Results[OTC_params->OTC + 1][10] = pd_fee;
      Results[OTC_params->OTC + 1][11] = pd_fee_flat;
    } else {
    }

    for (realisation = 0; realisation < OTC_params->COPULAnSimul;
         realisation++) {
      /* Smile */
      callRealisation = Pay_Rec_Mult * ((pd_leg[realisation] + pd_fee) -
                                        (fund_leg[realisation] + fund_fee));
      if (callRealisation > 0.0) {
        if (Calc_flat)
          Results[OTC_params->OTC + 1][2] += callRealisation;
        else {
          Results[2][0] += callRealisation;
          Results[3][0] += callRealisation * callRealisation;
        }
      }

      callRealisation =
          callRealisation -
          Pay_Rec_Mult * (((pd_legShort[realisation] + pd_feeShort) -
                           (fund_legShort[realisation] + fund_feeShort)));
      if (callRealisation > 0.0) {
        if (Calc_flat)
          Results[OTC_params->OTC + 1][4] += callRealisation;
        else {
          Results[2][2] += callRealisation;
          Results[3][2] += callRealisation * callRealisation;
        }
      }

      callRealisation =
          Pay_Rec_Mult * ((pd_legShort[realisation] + pd_feeShort) -
                          (fund_legShort[realisation] + fund_feeShort));
      if (callRealisation > 0.0) {
        if (Calc_flat)
          Results[OTC_params->OTC + 1][3] += callRealisation;
        else {
          Results[2][1] += callRealisation;
          Results[3][1] += callRealisation * callRealisation;
        }
      }

      if (Calc_flat) {
        /* 3F Smile */
        callRealisation =
            Pay_Rec_Mult * ((pd_leg_flat[realisation] + pd_fee_flat) -
                            (fund_leg_flat[realisation] + fund_fee_flat));
        if (callRealisation > 0.0)
          Results[OTC_params->OTC + 1][6] += callRealisation;

        callRealisation =
            callRealisation -
            Pay_Rec_Mult *
                (((pd_legShort_flat[realisation] + pd_feeShort_flat) -
                  (fund_legShort_flat[realisation] + fund_feeShort_flat)));
        if (callRealisation > 0.0)
          Results[OTC_params->OTC + 1][8] += callRealisation;

        callRealisation =
            Pay_Rec_Mult *
            ((pd_legShort_flat[realisation] + pd_feeShort_flat) -
             (fund_legShort_flat[realisation] + fund_feeShort_flat));
        if (callRealisation > 0.0)
          Results[OTC_params->OTC + 1][7] += callRealisation;
      }
    }

    if (Calc_flat) {
      Results[OTC_params->OTC + 1][1] *= df;
      Results[OTC_params->OTC + 1][2] *= df / OTC_params->COPULAnSimul;
      Results[OTC_params->OTC + 1][3] *= df / OTC_params->COPULAnSimul;
      Results[OTC_params->OTC + 1][4] *= df / OTC_params->COPULAnSimul;
      Results[OTC_params->OTC + 1][5] *= df;
      Results[OTC_params->OTC + 1][6] *= df / OTC_params->COPULAnSimul;
      Results[OTC_params->OTC + 1][7] *= df / OTC_params->COPULAnSimul;
      Results[OTC_params->OTC + 1][8] *= df / OTC_params->COPULAnSimul;
      Results[OTC_params->OTC + 1][9] *= df;
      Results[OTC_params->OTC + 1][10] *= df;
      Results[OTC_params->OTC + 1][11] *= df;
    } else {
      Results[0][0] *= df;
      Results[1][0] *= df;
      Results[2][0] *= df / OTC_params->COPULAnSimul;
      Results[3][0] *=
          df / OTC_params->COPULAnSimul / sqrt(OTC_params->COPULAnSimul);

      Results[0][1] *= df;
      Results[1][1] *= df;
      Results[2][1] *= df / OTC_params->COPULAnSimul;
      Results[3][1] *=
          df / OTC_params->COPULAnSimul / sqrt(OTC_params->COPULAnSimul);

      Results[0][2] = Results[0][0] - Results[0][1];
      Results[1][2] = Results[1][0] - Results[1][1];
      Results[2][2] *= df / OTC_params->COPULAnSimul;
      Results[3][2] *=
          df / OTC_params->COPULAnSimul / sqrt(OTC_params->COPULAnSimul);
    }
  } else {
    if (Calc_flat) {
      Results[OTC_params->OTC + 1][0] = tstarTime;
      Results[OTC_params->OTC + 1][1] = Pay_Rec_Mult * (pd_fee - fund_fee);
      Results[OTC_params->OTC + 1][2] = 0.0;
      Results[OTC_params->OTC + 1][3] = 0.0;
      Results[OTC_params->OTC + 1][4] = 0.0;
      Results[OTC_params->OTC + 1][5] = Pay_Rec_Mult * (pd_fee_flat - fund_fee);
      Results[OTC_params->OTC + 1][6] = 0.0;
      Results[OTC_params->OTC + 1][7] = 0.0;
      Results[OTC_params->OTC + 1][8] = 0.0;
      Results[OTC_params->OTC + 1][9] = 0.0;
      Results[OTC_params->OTC + 1][10] = 0.0;
      Results[OTC_params->OTC + 1][11] = 0.0;
    } else {
      Results[0][0] += fund_fee;
      Results[1][0] += pd_fee;
      Results[0][1] += fund_feeShort;
      Results[1][1] += pd_feeShort;

      Results[2][0] = 0.0;
      Results[3][0] = 0.0;

      Results[2][1] = 0.0;
      Results[3][1] = 0.0;

      Results[0][2] = 0.0;
      Results[1][2] = 0.0;
      Results[2][2] = 0.0;
      Results[3][2] = 0.0;
    }
  }

  /* To avoid numerical imprecision */
  if (!Calc_flat) {
    if (Results[2][0] < Pay_Rec_Mult * (Results[1][0] - Results[0][0]))
      Results[2][0] =
          Pay_Rec_Mult * (Results[1][0] - Results[0][0]) * (1.0 + 1.0e-10);

    if (Results[2][1] < Pay_Rec_Mult * (Results[1][1] - Results[0][1]))
      Results[2][1] =
          Pay_Rec_Mult * (Results[1][1] - Results[0][1]) * (1.0 + 1.0e-10);

    if (Results[2][2] < Pay_Rec_Mult * (Results[1][2] - Results[0][2]))
      Results[2][2] =
          Pay_Rec_Mult * (Results[1][2] - Results[0][2]) * (1.0 + 1.0e-10);
  }

FREE_RETURN:
  if (smile_params_flat)
    free(smile_params_flat);

  if (fund_leg)
    free(fund_leg);
  if (pd_leg)
    free(pd_leg);
  if (fund_legShort)
    free(fund_legShort);
  if (pd_legShort)
    free(pd_legShort);

  if (fund_leg_flat)
    free(fund_leg_flat);
  if (pd_leg_flat)
    free(pd_leg_flat);
  if (fund_legShort_flat)
    free(fund_legShort_flat);
  if (pd_legShort_flat)
    free(pd_legShort_flat);

  if (fundDf)
    free(fundDf);
  if (PdDfDom)
    free(PdDfDom);
  if (FwdDfDom)
    free(FwdDfDom);
  if (FwdDfFor)
    free(FwdDfFor);

  if (strikeCap)
    free(strikeCap);
  if (strikeFloor)
    free(strikeFloor);

  if (PDdfReal)
    free(PDdfReal);
  if (PDFxDomdfReal)
    free(PDFxDomdfReal);
  if (PDFxFordfReal)
    free(PDFxFordfReal);

  if (NewCoupons)
    free(NewCoupons);

  return err;
}

Err OTCfundingValuation(double domReal, double forReal, double fxReal,
                        CPD_STR cpd, CPD_UND und, otcpd_precalc *OTC_precalc,
                        double tstarTime, double phiDom, double phiFor, int otc,
                        double dfNotExchge, double dfNotExchgeShort,
                        double dfFirstFunding, double *fundDf, double *fund_leg,
                        double *fund_legShort) {
  double df;
  int fund_otc_start, fund_otc_short_end;
  int i;

  double real, fund_mult, lda, phi;

  Err err = NULL;

  if (cpd->fund_leg->dom_for == 0) {
    real = domReal;
    fund_mult = 1.0;
    lda = und->lda_dom;
    phi = phiDom;
  } else {
    real = forReal;
    fund_mult = fxReal;
    lda = und->lda_for;
    phi = phiFor;
  }

  fund_otc_start = cpd->call[otc].fund_idx;
  if (otc < cpd->num_calls - 1) {
    fund_otc_short_end = cpd->call[otc + 1].fund_idx;
  }

  /* ==================================
  0)	Notional exchange at Tset and last coupon date
  ================================== */
  /* Initial Notional */
  df =
      exp(OTC_precalc->df_const_fund_pay[otc][cpd->fund_leg->num_cpn + 1] +
          OTC_precalc->df_lin_fund_pay[otc][cpd->fund_leg->num_cpn + 1] * real);

  *fund_leg = -dfNotExchge * df * cpd->fund_leg->notional;

  if (otc < cpd->num_calls - 1) {
    *fund_legShort = *fund_leg;

    /*	Final Notional For Short option*/
    df = exp(OTC_precalc->df_const_fund_pay[otc][cpd->fund_leg->num_cpn + 2] +
             OTC_precalc->df_lin_fund_pay[otc][cpd->fund_leg->num_cpn + 2] *
                 real);

    *fund_legShort += dfNotExchgeShort * df * cpd->fund_leg->notional;
  }

  /* ==================================
  1)	Value funding leg starting at call date
  ================================== */
  /* B(t  ,Ts) */
  df = exp(OTC_precalc->df_const_fund_start[otc] +
           OTC_precalc->df_lin_fund_start[otc] * real);

  *fund_leg += dfFirstFunding * df * cpd->fund_leg->notional;

  if (otc < cpd->num_calls - 1)
    *fund_legShort += dfFirstFunding * df * cpd->fund_leg->notional;

  for (i = fund_otc_start; i < cpd->fund_leg->num_cpn; i++) {
    /* ==================================
    coupon: spread + margin
    ================================== */
    df = exp(OTC_precalc->df_const_fund_pay[otc][i] +
             OTC_precalc->df_lin_fund_pay[otc][i] * real);

    /*  B(t  ,Te) */
    *fund_leg += fundDf[i] * df * cpd->fund_leg->cpn[i].cpn;

    if (otc < cpd->num_calls - 1) {
      if (i >= fund_otc_start && i < fund_otc_short_end) {
        *fund_legShort += fundDf[i] * df * cpd->fund_leg->cpn[i].cpn;
      }
      if (i == fund_otc_short_end - 1)
        *fund_legShort -= fundDf[i] * df * cpd->fund_leg->notional;
    }
  }

  if (otc == cpd->num_calls - 1)
    *fund_legShort = *fund_leg;

  *fund_leg *= fund_mult;
  *fund_legShort *= fund_mult;

  return err;
}

Err OTCNewfundingValuation(double domReal, double forReal, double fxReal,
                           CPD_STR cpd, CPD_UND und, otcpd_precalc *OTC_precalc,
                           double tstarTime, double phiDom, double phiFor,
                           int otc, double dfNotExchge, double dfNotExchgeShort,
                           double dfFirstFunding, double *NewCoupons,
                           double *fundDfDom, double *fundDfFor,
                           double *fund_leg, double *fund_legShort) {
  double df;
  int fund_otc_start;
  int pd_otc_start, pd_otc_short_end;
  int i;

  double real, fund_mult, lda, phi;

  Err err = NULL;

  if (cpd->fund_leg->dom_for == 0) {
    real = domReal;
    fund_mult = 1.0;
    lda = und->lda_dom;
    phi = phiDom;
  } else {
    real = forReal;
    fund_mult = fxReal;
    lda = und->lda_for;
    phi = phiFor;
  }

  pd_otc_start = cpd->call[otc].pd_idx;
  if (otc < cpd->num_calls - 1) {
    pd_otc_short_end = cpd->call[otc + 1].pd_idx;
  }

  fund_otc_start = cpd->call[otc].fund_idx;

  /* ==================================
  0)	Notional exchange at Tset and last coupon date
  ================================== */
  /* Initial Notional */
  df =
      exp(OTC_precalc->df_const_fund_pay[otc][cpd->fund_leg->num_cpn + 1] +
          OTC_precalc->df_lin_fund_pay[otc][cpd->fund_leg->num_cpn + 1] * real);

  *fund_leg = -dfNotExchge * df * cpd->fund_leg->notional;

  if (otc < cpd->num_calls - 1) {
    *fund_legShort = *fund_leg;

    /*	Final Notional For Short option*/
    df = exp(OTC_precalc->df_const_fund_pay[otc][cpd->fund_leg->num_cpn + 2] +
             OTC_precalc->df_lin_fund_pay[otc][cpd->fund_leg->num_cpn + 2] *
                 real);

    *fund_legShort += dfNotExchgeShort * df * cpd->fund_leg->notional;
  }

  /* ==================================
  1)	Value funding leg starting at call date
  ================================== */
  /* B(t  ,Ts) */
  df = exp(OTC_precalc->df_const_fund_start[otc] +
           OTC_precalc->df_lin_fund_start[otc] * real);

  *fund_leg += dfFirstFunding * df * cpd->fund_leg->notional;

  if (otc < cpd->num_calls - 1)
    *fund_legShort += dfFirstFunding * df * cpd->fund_leg->notional;

  if (cpd->fund_leg->dom_for == 0) {
    for (i = pd_otc_start; i < cpd->pd_leg->num_cpn; i++) {
      /* ==================================
      coupon: spread + margin
      ================================== */
      *fund_leg += fundDfDom[i] * NewCoupons[i];

      if (otc < cpd->num_calls - 1) {
        if (i >= pd_otc_start && i < pd_otc_short_end) {
          *fund_legShort += fundDfDom[i] * NewCoupons[i];
        }
        if (i == pd_otc_short_end - 1)
          *fund_legShort -= fundDfDom[i] * cpd->fund_leg->notional;
      }
    }
  } else {
    for (i = pd_otc_start; i < cpd->pd_leg->num_cpn; i++) {
      /* ==================================
      coupon: spread + margin
      ================================== */
      *fund_leg += fundDfFor[i] * NewCoupons[i];

      if (otc < cpd->num_calls - 1) {
        if (i >= pd_otc_start && i < pd_otc_short_end) {
          *fund_legShort += fundDfFor[i] * NewCoupons[i];
        }
        if (i == pd_otc_short_end - 1)
          *fund_legShort -= fundDfFor[i] * cpd->fund_leg->notional;
      }
    }
  }

  if (otc == cpd->num_calls - 1)
    *fund_legShort = *fund_leg;

  *fund_leg *= fund_mult;
  *fund_legShort *= fund_mult;

  return err;
}

/* ========================================================
Calculates the cumulative correlation between two short/long IV
========================================================== */
Err OTCPDcorrel(CPD_STR cpd, CPD_UND und, double pd_not,
                double *pd_fwdsmile_vol, double *pd_fwdsmile_alpha,
                double *pd_fwdsmile_beta, double *pd_fwdsmile_rho, int smileFee,
                int otc, int nSimul, double **matrix, double tstarTime,
                int tstarDate, double *FxAdj, double *FeeFxAdj, int firstIndex,
                int secondIndex, int firstLong, int secondLong,
                long correlTstar, double **Results) {
  int i, ii, realisation;

  /* For the payoff */
  double *strikeCap = NULL, *strikeFloor = NULL;
  double forward;
  int type;
  double smile_std;
  double smile_half_std;
  double cap;
  double floor;
  PD_EXO_CPN cpn;
  SrtCallPutType InstCallType;
  double *fundDf = NULL, *PdDfDom = NULL, *FwdDfDom = NULL, *FwdDfFor = NULL;
  double discount, df, dfDom, dfFor, dfFirstFunding, dfNotExchge,
      dfFinalNotExchge, df0T, df0Tstar;

  /* For the correl */
  double *fund_leg = NULL, *pd_leg = NULL;
  double *fund_leg1 = NULL, *pd_leg1 = NULL;
  int pd_otc_start,
      pd_otc_end; /*	index of the first coupon to be called */
  int fund_otc_start, fund_otc_end;
  int Index, Long;
  double correlTstarTime, dfTstar, *dfTstarVec = NULL;
  double IV1, IV2;

  /* For the Df */
  double phiDom, phiFor;

  /* For The output */
  double Vol1, Vol2, Cov12;

  Err err = NULL;

  /* =============================
  memory allocation
  ============================= */
  fund_leg = calloc(nSimul, sizeof(double));
  pd_leg = calloc(nSimul, sizeof(double));
  fund_leg1 = calloc(nSimul, sizeof(double));
  pd_leg1 = calloc(nSimul, sizeof(double));
  dfTstarVec = calloc(nSimul, sizeof(double));
  fundDf = calloc(cpd->fund_leg->num_cpn, sizeof(double));
  PdDfDom = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  FwdDfDom = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  FwdDfFor = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  strikeCap = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  strikeFloor = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));

  if (!fund_leg || !pd_leg || !fundDf || !PdDfDom || !FwdDfDom || !FwdDfFor ||
      !strikeCap || !strikeFloor) {
    err = "OTCcaller: memory allocation error";
    goto FREE_RETURN;
  }

  /* =============================
  Precalculations
  ============================= */
  /* Tests */
  if (firstIndex < otc || secondIndex < otc ||
      firstIndex > cpd->num_calls - 1 || secondIndex > cpd->num_calls - 1) {
    err = "Invalid first or second Index";
    goto FREE_RETURN;
  }

  correlTstarTime = (double)(correlTstar - und->today) / DAYS_IN_YEAR;

  dfTstar = swp_f_df(und->today, correlTstar, und->dom_yc) /
            swp_f_df(und->today, cpd->call[otc].ex_date, und->dom_yc);

  err = OTCdfPhi(und->lda_dom, cpd->call[otc].ex_time, und->sigma_time_rates,
                 und->sigma_dom, und->sigma_n_rates, &phiDom);

  if (err)
    goto FREE_RETURN;

  err = OTCdfPhi(und->lda_for, cpd->call[otc].ex_time, und->sigma_time_rates,
                 und->sigma_for, und->sigma_n_rates, &phiFor);

  if (err)
    goto FREE_RETURN;

  for (ii = 0; ii < 2; ii++) {
    if (ii == 0) {
      Index = firstIndex;
      Long = firstLong;
    } else {
      Index = secondIndex;
      Long = secondLong;
    }

    pd_otc_start = cpd->call[Index].pd_idx;
    fund_otc_start = cpd->call[Index].fund_idx;

    if (Long || Index == cpd->num_calls - 1) {
      pd_otc_end = cpd->pd_leg->num_cpn;
      fund_otc_end = cpd->fund_leg->num_cpn;
    } else {
      pd_otc_end = cpd->call[Index + 1].pd_idx;
      fund_otc_end = cpd->call[Index + 1].fund_idx;
    }

    dfDom = swp_f_df(und->today, cpd->call[otc].ex_date, und->dom_yc);
    dfFor = swp_f_df(und->today, cpd->call[otc].ex_date, und->for_yc);

    dfFirstFunding =
        swp_f_df(und->today, cpd->fund_leg->cpn[fund_otc_start].start_date,
                 und->dom_yc) /
        dfDom;
    dfNotExchge =
        swp_f_df(und->today, cpd->call[Index].set_date, und->dom_yc) / dfDom;

    if (Index < cpd->num_calls - 1)
      dfFinalNotExchge =
          swp_f_df(und->today, cpd->call[Index + 1].set_date, und->dom_yc) /
          dfDom;
    else
      dfFinalNotExchge =
          swp_f_df(und->today, cpd->pd_leg->cpn[pd_otc_end - 1].pay_date,
                   und->dom_yc) /
          dfDom;

    for (i = fund_otc_start; i < fund_otc_end; i++)
      fundDf[i] =
          swp_f_df(und->today, cpd->fund_leg->cpn[i].pay_date, und->dom_yc) /
          dfDom;

    for (i = pd_otc_start; i < pd_otc_end; i++) {
      PdDfDom[i] =
          swp_f_df(und->today, cpd->pd_leg->cpn[i].pay_date, und->dom_yc) /
          dfDom;
      FwdDfDom[i] =
          swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->dom_yc) /
          dfDom;
      FwdDfFor[i] =
          swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->for_yc) /
          dfFor;
      if (fabs(cpd->pd_leg->cpn[i].beta) > 1.0e-16) {
        strikeCap[i] = (cpd->pd_leg->cpn[i].cap - cpd->pd_leg->cpn[i].alpha) /
                       cpd->pd_leg->cpn[i].beta;
        strikeFloor[i] =
            (cpd->pd_leg->cpn[i].floor - cpd->pd_leg->cpn[i].alpha) /
            cpd->pd_leg->cpn[i].beta;
      }
    }

    PdDfDom[cpd->pd_leg->num_cpn] =
        swp_f_df(und->today, cpd->pd_leg->not_ref.pay_date, und->dom_yc) /
        dfDom;
    FwdDfDom[cpd->pd_leg->num_cpn] =
        swp_f_df(und->today, cpd->pd_leg->not_ref.fx_val_date, und->dom_yc) /
        dfDom;
    FwdDfFor[cpd->pd_leg->num_cpn] =
        swp_f_df(und->today, cpd->pd_leg->not_ref.fx_val_date, und->for_yc) /
        dfFor;
    if (fabs(cpd->pd_leg->not_ref.beta) > 1.0e-16) {
      strikeCap[cpd->pd_leg->num_cpn] =
          (cpd->pd_leg->not_ref.cap - cpd->pd_leg->not_ref.alpha) /
          cpd->pd_leg->not_ref.beta;
      strikeFloor[cpd->pd_leg->num_cpn] =
          (cpd->pd_leg->not_ref.floor - cpd->pd_leg->not_ref.alpha) /
          cpd->pd_leg->not_ref.beta;
    }

    /* ==========================================================
            PAYOFF
    ========================================================== */
    for (realisation = 0; realisation < nSimul; realisation++) {
      /* ==================================
      0)	Notional exchange at Tset and last coupon date
      ================================== */
      /* Initial Notional */
      err = OTCdf(matrix[realisation][0], und->lda_dom, cpd->call[otc].ex_time,
                  cpd->call[Index].set_time, tstarTime, phiDom, &df);

      if (err)
        goto FREE_RETURN;

      fund_leg[realisation] = -dfNotExchge * df * cpd->fund_leg->notional;
      pd_leg[realisation] = -dfNotExchge * df * pd_not;

      if (pd_otc_end == cpd->pd_leg->num_cpn) {
        /*	Final Notional For Long option*/
        cpn = &(cpd->pd_leg->not_ref);

        /*	Discount */
        err =
            OTCdf(matrix[realisation][0], und->lda_dom, cpd->call[otc].ex_time,
                  cpn->pay_time, tstarTime, phiDom, &df);

        if (err)
          goto FREE_RETURN;

        discount = PdDfDom[cpd->pd_leg->num_cpn] * df;

        /*	Fwd fx */

        if (fabs(cpn->beta) > 1.0e-16) {
          /* ==================================
          Forward @ Tval
          ================================== */
          err = OTCdf(matrix[realisation][0], und->lda_dom,
                      cpd->call[otc].ex_time, cpn->fx_val_time, tstarTime,
                      phiDom, &df);

          if (err)
            goto FREE_RETURN;

          dfDom = FwdDfDom[cpd->pd_leg->num_cpn] * df;

          err = OTCdf(matrix[realisation][1], und->lda_for,
                      cpd->call[otc].ex_time, cpn->fx_val_time, tstarTime,
                      phiFor, &df);

          if (err)
            goto FREE_RETURN;

          dfFor = FwdDfFor[cpd->pd_leg->num_cpn] * df;

          forward = FxAdj[cpd->pd_leg->num_cpn] * matrix[realisation][2] *
                    dfFor / dfDom;

          /* ==================================
          Floor
          ================================== */
          if (cpn->floored && fabs(cpn->beta) > 1.0e-16 &&
              cpn->floor > -1.0e-16) {
            if (cpn->beta > 0.0) {
              if (strikeFloor[cpd->pd_leg->num_cpn] > 1.0e-16)
                type = 4; /*	Put */
              else
                type = 2; /*	Put IV */

              InstCallType = SRT_PUT;
            } else {
              if (strikeFloor[cpd->pd_leg->num_cpn] > 1.0e-16)
                type = 3; /*	Call */
              else
                type = 1; /*	Call IV */

              InstCallType = SRT_CALL;
            }
          } else
            type = 0; /*	No floor */

          if (type >= 3) {
            err = srt_f_optsarbvol(forward, strikeFloor[cpd->pd_leg->num_cpn],
                                   cpn->fx_fix_time - cpd->call[otc].ex_time,
                                   pd_fwdsmile_vol[cpd->pd_leg->num_cpn],
                                   pd_fwdsmile_alpha[cpd->pd_leg->num_cpn],
                                   pd_fwdsmile_beta[cpd->pd_leg->num_cpn],
                                   pd_fwdsmile_rho[cpd->pd_leg->num_cpn],
                                   SRT_BETAVOL, SRT_LOGNORMAL, &smile_std);

            if (err)
              goto FREE_RETURN;

            smile_std *= sqrt(cpn->fx_fix_time - cpd->call[otc].ex_time);
            smile_half_std = 0.5 * smile_std;
          }

          if ((forward < 0.0 && type >= 3) || (smile_std < 1.0e-6 && type >= 3))
            type -= 2;

          floor =
              OPT_VAL_MACRO(type, forward, strikeFloor[cpd->pd_leg->num_cpn],
                            smile_std, smile_half_std);

          /* ==================================
          Cap
          ================================== */
          if (cpn->capped && fabs(cpn->beta) > 1.0e-16 && cpn->cap > -1.0e-16) {
            if (cpn->beta > 0.0) {
              if (strikeCap[cpd->pd_leg->num_cpn] > 1.0e-16)
                type = 3; /*	Call */
              else
                type = 1; /*	Call IV */

              InstCallType = SRT_CALL;
            } else {
              if (strikeCap[cpd->pd_leg->num_cpn] > 1.0e-16)
                type = 4; /*	Put */
              else
                type = 2; /*	Put IV */

              InstCallType = SRT_PUT;
            }
          } else
            type = 0; /*	No cap */

          if (type >= 3) {
            err = srt_f_optsarbvol(forward, strikeCap[cpd->pd_leg->num_cpn],
                                   cpn->fx_fix_time - cpd->call[otc].ex_time,
                                   pd_fwdsmile_vol[cpd->pd_leg->num_cpn],
                                   pd_fwdsmile_alpha[cpd->pd_leg->num_cpn],
                                   pd_fwdsmile_beta[cpd->pd_leg->num_cpn],
                                   pd_fwdsmile_rho[cpd->pd_leg->num_cpn],
                                   SRT_BETAVOL, SRT_LOGNORMAL, &smile_std);

            if (err)
              goto FREE_RETURN;

            smile_std *= sqrt(cpn->fx_fix_time - cpd->call[otc].ex_time);
            smile_half_std = 0.5 * smile_std;
          }

          if ((forward < 0.0 && type >= 3) || (smile_std < 1.0e-6 && type >= 3))
            type -= 2;

          cap = OPT_VAL_MACRO(type, forward, strikeCap[cpd->pd_leg->num_cpn],
                              smile_std, smile_half_std);

        } else {
          forward = floor = cap = 0.0;
        }

        /*	Coupon pv */
        pd_leg[realisation] += discount * (cpn->alpha + cpn->beta * forward +
                                           fabs(cpn->beta) * (floor - cap));
      } else {
        /*	Discount */
        err = OTCdf(
            matrix[realisation][0], und->lda_dom, cpd->call[otc].ex_time,
            cpd->pd_leg->cpn[pd_otc_end - 1].pay_time, tstarTime, phiDom, &df);

        if (err)
          goto FREE_RETURN;

        discount = PdDfDom[pd_otc_end - 1] * df;

        pd_leg[realisation] += discount * pd_not;
      }

      /* ==================================
      1)	Value funding leg starting at call date
      ================================== */
      /* B(t  ,Ts) */
      err = OTCdf(matrix[realisation][0], und->lda_dom, cpd->call[otc].ex_time,
                  cpd->fund_leg->cpn[fund_otc_start].start_time, tstarTime,
                  phiDom, &df);

      if (err)
        goto FREE_RETURN;

      fund_leg[realisation] += dfFirstFunding * df * cpd->fund_leg->notional;

      for (i = fund_otc_start; i < fund_otc_end; i++) {
        /* ==================================
        coupon: spread + margin
        ================================== */
        err =
            OTCdf(matrix[realisation][0], und->lda_dom, cpd->call[otc].ex_time,
                  cpd->fund_leg->cpn[i].pay_time, tstarTime, phiDom, &df);

        if (err)
          goto FREE_RETURN;

        /*  B(t  ,Te) */
        fund_leg[realisation] += fundDf[i] * df * cpd->fund_leg->cpn[i].cpn;
      }

      /* ==================================
      2)	Value power dual leg starting at call date
      ================================== */
      for (i = pd_otc_start; i < pd_otc_end; i++) {
        cpn = cpd->pd_leg->cpn + i;

        /* ==================================
        Forward @ Tval with Convexity adjustment
        ================================== */
        err =
            OTCdf(matrix[realisation][0], und->lda_dom, cpd->call[otc].ex_time,
                  cpn->fx_val_time, tstarTime, phiDom, &df);

        if (err)
          goto FREE_RETURN;

        dfDom = FwdDfDom[i] * df;

        err =
            OTCdf(matrix[realisation][1], und->lda_for, cpd->call[otc].ex_time,
                  cpn->fx_val_time, tstarTime, phiFor, &df);

        if (err)
          goto FREE_RETURN;

        dfFor = FwdDfFor[i] * df;

        forward = FxAdj[i] * matrix[realisation][2] * dfFor / dfDom;

        /* ==================================
        domestic DF @ Tpay
        ================================== */
        err =
            OTCdf(matrix[realisation][0], und->lda_dom, cpd->call[otc].ex_time,
                  cpn->pay_time, tstarTime, phiDom, &df);

        if (err)
          goto FREE_RETURN;

        dfDom = PdDfDom[i] * df;

        /* ==================================
        Floor
        ================================== */
        if (cpn->floored && fabs(cpn->beta) > 1.0e-16 &&
            cpn->floor > -1.0e-16) {
          if (cpn->beta > 0.0) {
            if (strikeFloor[i] > 1.0e-16)
              type = 4; /*	Put */
            else
              type = 2; /*	Put IV */

            InstCallType = SRT_PUT;
          } else {
            if (strikeFloor[i] > 1.0e-16)
              type = 3; /*	Call */
            else
              type = 1; /*	Call IV */

            InstCallType = SRT_CALL;
          }
        } else {
          if (cpn->floor <= -1.0e-16) {
            err = "Coupon floor must be positive";
            goto FREE_RETURN;
          } else
            type = 0; /*	No floor */
        }

        if (type >= 3) {
          err = srt_f_optsarbvol(forward, strikeFloor[i],
                                 cpn->fx_fix_time - cpd->call[otc].ex_time,
                                 pd_fwdsmile_vol[i], pd_fwdsmile_alpha[i],
                                 pd_fwdsmile_beta[i], pd_fwdsmile_rho[i],
                                 SRT_BETAVOL, SRT_LOGNORMAL, &smile_std);

          if (err)
            goto FREE_RETURN;

          smile_std *= sqrt(cpn->fx_fix_time - cpd->call[otc].ex_time);
          smile_half_std = 0.5 * smile_std;
        }

        if ((forward < 0.0 && type >= 3) || (smile_std < 1.0e-6 && type >= 3))
          type -= 2;

        floor = OPT_VAL_MACRO(type, forward, strikeFloor[i], smile_std,
                              smile_half_std);
        /* ==================================
        Cap
        ================================== */
        if (cpn->capped && fabs(cpn->beta) > 1.0e-16 && cpn->cap > -1.0e-16) {
          if (cpn->beta > 0.0) {
            if (strikeCap[i] > 1.0e-16)
              type = 3; /*	Call */
            else
              type = 1; /*	Call IV */

            InstCallType = SRT_CALL;
          } else {
            if (strikeCap[i] > 1.0e-16)
              type = 4; /*	Put */
            else
              type = 2; /*	Put IV */

            InstCallType = SRT_PUT;
          }
        } else {
          if (cpn->cap <= -1.0e-16) {
            err = "Coupon cap must be positive";
            goto FREE_RETURN;
          } else
            type = 0; /*	No cap */
        }

        if (type >= 3) {
          err = srt_f_optsarbvol(
              forward, strikeCap[i], cpn->fx_fix_time - cpd->call[otc].ex_time,
              pd_fwdsmile_vol[i], pd_fwdsmile_alpha[i], pd_fwdsmile_beta[i],
              pd_fwdsmile_rho[i], SRT_BETAVOL, SRT_LOGNORMAL, &smile_std);

          if (err)
            goto FREE_RETURN;

          smile_std *= sqrt(cpn->fx_fix_time - cpd->call[otc].ex_time);
          smile_half_std = 0.5 * smile_std;
        }

        if ((forward < 0.0 && type >= 3) || (smile_std < 1.0e-6 && type >= 3))
          type -= 2;

        cap = OPT_VAL_MACRO(type, forward, strikeCap[i], smile_std,
                            smile_half_std);

        pd_leg[realisation] += dfDom * (cpn->alpha + cpn->beta * forward +
                                        fabs(cpn->beta) * (floor - cap));
      }

      /* B(T  ,T*) */
      if (ii == 0) {
        err =
            OTCdf(matrix[realisation][0], und->lda_dom, cpd->call[otc].ex_time,
                  correlTstarTime, tstarTime, phiDom, &df);

        if (err)
          goto FREE_RETURN;

        dfTstarVec[realisation] = dfTstar * df;
      }
    }

    if (ii == 0) {
      for (realisation = 0; realisation < nSimul; realisation++) {
        fund_leg1[realisation] = fund_leg[realisation];
        pd_leg1[realisation] = pd_leg[realisation];
      }
    }
  }

  df0T = swp_f_df(und->today, tstarDate, und->dom_yc);
  df0Tstar = swp_f_df(und->today, correlTstar, und->dom_yc);
  df = df0T / df0Tstar;

  Vol1 = Vol2 = Cov12 = 0.0;

  for (realisation = 0; realisation < nSimul; realisation++) {
    Results[3][0] += fund_leg1[realisation];
    Results[4][0] += pd_leg1[realisation];
    Results[5][0] += fund_leg[realisation];
    Results[6][0] += pd_leg[realisation];

    Vol1 += (fund_leg1[realisation] - pd_leg1[realisation]) *
            (fund_leg1[realisation] - pd_leg1[realisation]) /
            dfTstarVec[realisation];
    Vol2 += (fund_leg[realisation] - pd_leg[realisation]) *
            (fund_leg[realisation] - pd_leg[realisation]) /
            dfTstarVec[realisation];
    Cov12 += (fund_leg1[realisation] - pd_leg1[realisation]) *
             (fund_leg[realisation] - pd_leg[realisation]) /
             dfTstarVec[realisation];
  }

  Vol1 /= nSimul;
  Vol2 /= nSimul;
  Cov12 /= nSimul;

  Results[3][0] /= nSimul;
  Results[4][0] /= nSimul;
  Results[5][0] /= nSimul;
  Results[6][0] /= nSimul;

  IV1 = df * (Results[3][0] - Results[4][0]);
  IV2 = df * (Results[5][0] - Results[6][0]);

  Vol1 *= df;
  Vol2 *= df;
  Cov12 *= df;

  Vol1 -= IV1 * IV1;
  Vol2 -= IV2 * IV2;
  Cov12 -= IV1 * IV2;

  Results[3][0] *= df0T;
  Results[4][0] *= df0T;
  Results[5][0] *= df0T;
  Results[6][0] *= df0T;

  Results[0][0] = Cov12 / sqrt(Vol1 * Vol2);
  Results[1][0] = sqrt(Vol1 / cpd->call[otc].ex_time);
  Results[2][0] = sqrt(Vol2 / cpd->call[otc].ex_time);

FREE_RETURN:
  if (fund_leg)
    free(fund_leg);
  if (pd_leg)
    free(pd_leg);

  if (fund_leg1)
    free(fund_leg1);
  if (pd_leg1)
    free(pd_leg1);

  if (dfTstarVec)
    free(dfTstarVec);

  if (fundDf)
    free(fundDf);
  if (PdDfDom)
    free(PdDfDom);
  if (FwdDfDom)
    free(FwdDfDom);
  if (FwdDfFor)
    free(FwdDfFor);

  if (strikeCap)
    free(strikeCap);
  if (strikeFloor)
    free(strikeFloor);

  return err;
}

Err OTCtestPayoff(CPD_UND und, double **matrix, int nSimul, double tstarTime,
                  int tstarDate, int nStrikes, double **Results) {

  int i, realisation;
  double eval;
  Err err = NULL;

  /* ==================================
  Value Fx Options for the 3 strikes
  ================================== */
  for (realisation = 0; realisation < nSimul; realisation++) {
    for (i = 0; i < nStrikes; i++) {
      eval = matrix[realisation][2] - Results[i][0];
      if (eval < 0.0)
        eval = 0.0;

      Results[i][1] += eval;
    }
  }

  for (i = 0; i < nStrikes; i++)
    Results[i][1] /= nSimul;

  return err;
}

Err fwdSABRpayoff(long today, long forwardFixDate, double forwardFixTime,
                  long valDate, double fixTime, double valTime, double forward,
                  double fwdsmile_vol, double fwdsmile_alpha,
                  double fwdsmile_beta, double fwdsmile_rho, double **matrix,
                  int nSimul, double tstarTime, int tstarDate, char *dom_yc,
                  char *for_yc, double *sigma_time_rates, int sigma_n_rates,
                  double *sigma_dom, double lda_dom, double *sigma_for,
                  double lda_for, double *sigma_time_fx, double *sigma_fx,
                  int sigma_n_fx, double *corr_times, double *correl_dom_for,
                  double *correl_dom_fx, double *correl_for_fx,
                  int corr_n_times, int nStrikes, double *Strikes, int isCall,
                  int OutputVol, double *Results) {

  int j, realisation;

  /* For the payoff */
  double fwdSmileVisuNstd = 0.5;
  double vol, localForward;
  int type; /* CALL */
  double smile_std;
  double smile_half_std;

  double PdDfDom, FwdDfDom, FwdDfFor;
  double df, dfDom, dfFor;

  /* For the Df */
  double phiDom, phiFor;

  Err err = NULL;

  if (isCall)
    type = 3;
  else
    type = 4;

  /* ==================================
  Precalculations
  ================================== */
  dfDom = swp_f_df(today, forwardFixDate, dom_yc);
  dfFor = swp_f_df(today, forwardFixDate, for_yc);

  PdDfDom = swp_f_df(today, valDate, dom_yc) / dfDom;
  FwdDfDom = swp_f_df(today, valDate, dom_yc) / dfDom;
  FwdDfFor = swp_f_df(today, valDate, for_yc) / dfFor;

  /* ==================================
  Value Fx Options for the 3 strikes
  ================================== */
  err = OTCdfPhi(lda_dom, forwardFixTime, sigma_time_rates, sigma_dom,
                 sigma_n_rates, &phiDom);

  if (err)
    goto FREE_RETURN;

  err = OTCdfPhi(lda_for, forwardFixTime, sigma_time_rates, sigma_for,
                 sigma_n_rates, &phiFor);

  if (err)
    goto FREE_RETURN;

  for (realisation = 0; realisation < nSimul; realisation++) {
    /* ==================================
    Forward @ Tval
    ================================== */
    err = OTCdf(matrix[realisation][0], lda_dom, forwardFixTime, valTime,
                tstarTime, phiDom, &df);

    if (err)
      goto FREE_RETURN;

    dfDom = FwdDfDom * df;

    err = OTCdf(matrix[realisation][1], lda_for, forwardFixTime, valTime,
                tstarTime, phiFor, &df);

    if (err)
      goto FREE_RETURN;

    dfFor = FwdDfFor * df;

    localForward = matrix[realisation][2] * dfFor / dfDom;

    /* ==================================
    Domestic Df @ Tpay
    ================================== */
    err = OTCdf(matrix[realisation][0], lda_dom, forwardFixTime, valTime,
                tstarTime, phiDom, &df);

    if (err)
      goto FREE_RETURN;

    dfDom = PdDfDom * df;

    /* ==================================
    LOOP ON THE STRIKES
    ================================== */
    for (j = 0; j < nStrikes; j++) {
      err = srt_f_optsarbvol(localForward, Strikes[j], fixTime - forwardFixTime,
                             fwdsmile_vol, fwdsmile_alpha, fwdsmile_beta,
                             fwdsmile_rho, SRT_BETAVOL, SRT_LOGNORMAL,
                             &smile_std);

      if (err)
        goto FREE_RETURN;

      smile_std *= sqrt(fixTime - forwardFixTime);
      smile_half_std = 0.5 * smile_std;

      if (localForward < 0.0 && type >= 3) {
        Results[j] += dfDom * OPT_VAL_MACRO(type - 2, localForward, Strikes[j],
                                            smile_std, smile_half_std);
      } else {
        Results[j] += dfDom * OPT_VAL_MACRO(type, localForward, Strikes[j],
                                            smile_std, smile_half_std);
      }
    }
  }

  df = swp_f_df(today, forwardFixDate, dom_yc);

  dfDom = swp_f_df(today, valDate, dom_yc);

  for (j = 0; j < nStrikes; j++)
    Results[j] *= df / nSimul;

  if (OutputVol) {
    for (j = 0; j < nStrikes; j++) {
      err = srt_f_optimpvol(Results[j], forward, Strikes[j], fixTime, dfDom,
                            isCall ? SRT_CALL : SRT_PUT, SRT_LOGNORMAL, &vol);

      if (err)
        Results[j] = 0.0;
      else
        Results[j] = vol;
    }
  }

FREE_RETURN:

  return err;
}

Err fwdSABRpayoffAlphaFudge(
    long today, long forwardFixDate, double forwardFixTime, long valDate,
    double fixTime, double valTime, double forward, double fwdsmile_vol,
    double fwdsmile_alpha, double fwdsmile_beta, double fwdsmile_rho,
    double **matrix, int nSimul, double tstarTime, int tstarDate, char *dom_yc,
    char *for_yc, double *sigma_time_rates, int sigma_n_rates,
    double *sigma_dom, double lda_dom, double *sigma_for, double lda_for,
    double *sigma_time_fx, double *sigma_fx, int sigma_n_fx, double *corr_times,
    double *correl_dom_for, double *correl_dom_fx, double *correl_for_fx,
    int corr_n_times, int nStrikes, double *Strikes, int isCall, int OutputVol,
    double SwitchLevel, double Proba, double AlphaFudge, double *Results) {

  int j, realisation;

  /* For the payoff */
  double fwdSmileVisuNstd = 0.5;
  double vol, localForward;
  int type; /* CALL */
  double smile_std;
  double smile_half_std;

  double PdDfDom, FwdDfDom, FwdDfFor;
  double df, dfDom, dfFor;

  /* For the Df */
  double phiDom, phiFor;
  double betaVolUp, betaVolDown, p;

  Err err = NULL;

  if (isCall)
    type = 3;
  else
    type = 4;

  /* ==================================
  Precalculations
  ================================== */
  dfDom = swp_f_df(today, forwardFixDate, dom_yc);
  dfFor = swp_f_df(today, forwardFixDate, for_yc);

  PdDfDom = swp_f_df(today, valDate, dom_yc) / dfDom;
  FwdDfDom = swp_f_df(today, valDate, dom_yc) / dfDom;
  FwdDfFor = swp_f_df(today, valDate, for_yc) / dfFor;

  /* ==================================
  Value Fx Options for the 3 strikes
  ================================== */
  err = OTCdfPhi(lda_dom, forwardFixTime, sigma_time_rates, sigma_dom,
                 sigma_n_rates, &phiDom);

  if (err)
    goto FREE_RETURN;

  err = OTCdfPhi(lda_for, forwardFixTime, sigma_time_rates, sigma_for,
                 sigma_n_rates, &phiFor);

  if (err)
    goto FREE_RETURN;

  for (realisation = 0; realisation < nSimul; realisation++) {
    /* ==================================
    Forward @ Tval
    ================================== */
    err = OTCdf(matrix[realisation][0], lda_dom, forwardFixTime, valTime,
                tstarTime, phiDom, &df);

    if (err)
      goto FREE_RETURN;

    dfDom = FwdDfDom * df;

    err = OTCdf(matrix[realisation][1], lda_for, forwardFixTime, valTime,
                tstarTime, phiFor, &df);

    if (err)
      goto FREE_RETURN;

    dfFor = FwdDfFor * df;

    localForward = matrix[realisation][2] * dfFor / dfDom;

    /* ==================================
    Domestic Df @ Tpay
    ================================== */
    err = OTCdf(matrix[realisation][0], lda_dom, forwardFixTime, valTime,
                tstarTime, phiDom, &df);

    if (err)
      goto FREE_RETURN;

    dfDom = PdDfDom * df;

    /* ==================================
    LOOP ON THE STRIKES
    ================================== */
    for (j = 0; j < nStrikes; j++) {
      if (matrix[realisation][2] > SwitchLevel)
        p = Proba;
      else
        p = 1.0 - Proba;

      betaVolUp = fwdsmile_vol * exp(AlphaFudge * sqrt(forwardFixTime));
      betaVolDown = fwdsmile_vol * exp(-AlphaFudge * sqrt(forwardFixTime));

      err = srt_f_optsarbvol(localForward, Strikes[j], fixTime - forwardFixTime,
                             betaVolUp, fwdsmile_alpha, fwdsmile_beta,
                             fwdsmile_rho, SRT_BETAVOL, SRT_LOGNORMAL,
                             &smile_std);

      if (err)
        goto FREE_RETURN;

      smile_std *= sqrt(fixTime - forwardFixTime);
      smile_half_std = 0.5 * smile_std;

      if (localForward < 0.0 && type >= 3) {
        Results[j] += p * dfDom *
                      OPT_VAL_MACRO(type - 2, localForward, Strikes[j],
                                    smile_std, smile_half_std);
      } else {
        Results[j] += p * dfDom *
                      OPT_VAL_MACRO(type, localForward, Strikes[j], smile_std,
                                    smile_half_std);
      }

      err = srt_f_optsarbvol(localForward, Strikes[j], fixTime - forwardFixTime,
                             betaVolDown, fwdsmile_alpha, fwdsmile_beta,
                             fwdsmile_rho, SRT_BETAVOL, SRT_LOGNORMAL,
                             &smile_std);

      if (err)
        goto FREE_RETURN;

      smile_std *= sqrt(fixTime - forwardFixTime);
      smile_half_std = 0.5 * smile_std;

      if (localForward < 0.0 && type >= 3) {
        Results[j] += (1.0 - p) * dfDom *
                      OPT_VAL_MACRO(type - 2, localForward, Strikes[j],
                                    smile_std, smile_half_std);
      } else {
        Results[j] += (1.0 - p) * dfDom *
                      OPT_VAL_MACRO(type, localForward, Strikes[j], smile_std,
                                    smile_half_std);
      }
    }
  }

  df = swp_f_df(today, forwardFixDate, dom_yc);

  dfDom = swp_f_df(today, valDate, dom_yc);

  for (j = 0; j < nStrikes; j++)
    Results[j] *= df / nSimul;

  if (OutputVol) {
    for (j = 0; j < nStrikes; j++) {
      err = srt_f_optimpvol(Results[j], forward, Strikes[j], fixTime, dfDom,
                            isCall ? SRT_CALL : SRT_PUT, SRT_LOGNORMAL, &vol);

      if (err)
        Results[j] = 0.0;
      else
        Results[j] = vol;
    }
  }

FREE_RETURN:

  return err;
}

/* ====================================================
====================================================
KO
====================================================
==================================================== */
Err KO_get_coupon_paid(int Coupon_idx, CPD_STR cpd, double fx,
                       double *Coupon_value) {
  int str_idx;
  double call, cap, floor, opt_string;
  PD_EXO_CPN cpn;

  Err err = NULL;

  cpn = cpd->pd_leg->cpn + Coupon_idx;

  // Sum of past PD Coupon
  if (!cpn->use_opt_str) {
    if (cpn->floored && fabs(cpn->beta) > 1.0e-10) {
      if (cpn->beta > 0.0)
        floor = (cpn->floor - cpn->alpha) / cpn->beta - fx; /*	Put IV */
      else
        floor = fx - (cpn->floor - cpn->alpha) / cpn->beta; /*	Call IV */

      if (floor < 0.0)
        floor = 0.0;
    } else
      floor = 0.0;

    if (cpn->capped && fabs(cpn->beta) > 1.0e-10) {
      if (cpn->beta < 0.0)
        cap = (cpn->cap - cpn->alpha) / cpn->beta - fx; /*	Put IV */
      else
        cap = fx - (cpn->cap - cpn->alpha) / cpn->beta; /*	Call IV */

      if (cap < 0.0)
        cap = 0.0;
    } else
      cap = 0.0;

    *Coupon_value =
        cpn->alpha + cpn->beta * fx + fabs(cpn->beta) * (floor - cap);
  } else {
    opt_string = 0.0;
    for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
      if (fabs(cpn->weights[str_idx]) >
          1e-16 * cpd->call[Coupon_idx].pd_not_amt) {
        call = fx - cpn->strikes[str_idx];
        if (call > 0.0)
          opt_string += cpn->weights[str_idx] * call;
      }
    }
    *Coupon_value = cpn->wcst + cpn->wspot * fx + opt_string;
  }

  return err;
}

Err KOPDpayoffFee(CPD_STR cpd, CPD_UND und, otcpd_params *OTC_params,
                  otcpd_precalc *OTC_precalc, SMILE_VOL_MARKET smile_mkt,
                  SMILE_PARAMETERS smile_params, long start_date, double pd_not,
                  double tstarTime, int tstarDate, double **Results,
                  double ***save_values) {

  int i, j, realisation;
  int *KnockedOut = NULL;
  double *Payoff_fund = NULL, *Payoff_pd = NULL;
  double *KO_fund = NULL, *KO_pd = NULL;
  double phiDom, phiFor;

  /* For the payoff */
  double *strikeCap = NULL, *strikeFloor = NULL;
  double forward;
  double cap;
  double floor;
  PD_EXO_CPN cpn;
  double *fundDf = NULL, *PdDfDom = NULL, *FwdDfDom = NULL, *FwdDfFor = NULL;
  double discount, df, dfDom, dfFor, dfFund, dfFirstFunding, dfNotExchgePd,
      dfNotExchgeShortPd, dfNotExchgeFund, dfNotExchgeShortFund,
      dfStartNotExchgePd, dfStartNotExchgeFund, dfNotExchgeEndPd,
      dfNotExchgeEndFund;
  double start_time;

  /* For the call */
  double *fund_leg = NULL, *pd_leg = NULL;
  double *fund_legShort = NULL, *pd_legShort = NULL;
  int pd_otc_start,
      pd_otc_short_end; /*	index of the first coupon to be called */
  int fund_otc_start, fund_otc_short_end;

  /* For the fee adjustment */
  double fund_fee, pd_fee, fund_fee_mc, pd_fee_mc;
  double ko_fund_fee, ko_pd_fee, ko_fund_fee_mc, ko_pd_fee_mc;
  double coupon, coupon_fund, coupon_pd;
  double ko_coupon_fund, ko_coupon_pd;
  double exchange_fund, exchange_pd;

  /* For the Std */
  double PriceStd;

  /* For foreign funding */
  char *fund_yc;
  double fund_mult;

  /* For the dollat notional exchange */
  double dfNotEx[4];

  /* Number of paths for which the coupons will be paid */
  long Ncoupons = OTC_params->COPULAnSimul;

  /* Barrier adjustment */
  double Barrier[256];

  /* For string of options */
  int str_idx, type;
  double opt_string;

  /* To calculate until the std is under a certain level */
  double fmc_current_std = 1.0;
  double fmc_current_mean = 0.0;
  double fmc_current_var = 0.0;
  int fmc_nb_realisations = 0;
  int fmc_nb_coupons = 0;

  /* For the TARN */
  double Sum_coupon_paid, fx;

  double OptionValue;

  Err err = NULL;

  /* =============================
  Populates a vector with 0-1 => 0 : Pay Coupons; 1 : Exchange Notionals; 2 :
  allready KO in the past
  ============================= */
  KnockedOut = calloc(OTC_params->COPULAnSimul, sizeof(int));
  Payoff_fund = calloc(OTC_params->COPULAnSimul, sizeof(double));
  Payoff_pd = calloc(OTC_params->COPULAnSimul, sizeof(double));

  if (!KnockedOut || !Payoff_fund || !Payoff_pd) {
    err = "KOPDpayoff: memory allocation failure";
    if (err)
      goto FREE_RETURN;
  }

  if (!cpd->call[0].TARN_Do) {
    /* Precalc Barrier */
    for (j = 0; j < OTC_params->OTC; j++)
      Barrier[j] = cpd->call[j].orig_barrier /
                   OTC_precalc->KO_barrier[OTC_params->OTC][j];

    Barrier[OTC_params->OTC] = cpd->call[OTC_params->OTC].orig_barrier;

    /* Precalc KO */
    for (realisation = 0; realisation < OTC_params->COPULAnSimul;
         realisation++) {
      for (j = 0; j < OTC_params->OTC && !KnockedOut[realisation]; j++) {
        if (OTC_params->KO_FxMatrix[realisation][j] >
            Barrier[j]) // change the barrier /* To do */
        {
          if (!OTC_params->KO_do_optim ||
              (OTC_params->KO_do_optim && cpd->call[j].call_type)) {
            KnockedOut[realisation] = 2;
            Ncoupons--;
          }
        }
      }

      if (!KnockedOut[realisation]) {
        if (OTC_params->KO_FxMatrix[realisation][OTC_params->OTC] >
            Barrier[OTC_params->OTC]) {
          if (!OTC_params->KO_do_optim ||
              (OTC_params->KO_do_optim &&
               cpd->call[OTC_params->OTC].call_type)) {
            KnockedOut[realisation] = 1;
            Ncoupons--; /* To calculate the proba to pay a coupon */
          }
        } else {
          // Ncoupons --;			/* To calculate the proba of an exchange of
          // notional */
        }
      }
    }
  } else {
    /* Precalc TARN */
    for (realisation = 0; realisation < OTC_params->COPULAnSimul;
         realisation++) {
      Sum_coupon_paid = 0.0;
      for (j = 0; j < OTC_params->OTC && !KnockedOut[realisation]; j++) {
        fx = OTC_params->KO_FxMatrix[realisation][j] *
             OTC_precalc->KO_barrier[OTC_params->OTC][j];

        err = KO_get_coupon_paid(j, cpd, fx, &coupon);
        if (err)
          goto FREE_RETURN;

        Sum_coupon_paid += coupon;

        if (Sum_coupon_paid > cpd->call[j].orig_barrier) {
          KnockedOut[realisation] = 2;
          Ncoupons--;
        }
      }

      if (!KnockedOut[realisation]) {
        fx = OTC_params->KO_FxMatrix[realisation][OTC_params->OTC];

        err = KO_get_coupon_paid(OTC_params->OTC, cpd, fx, &coupon);
        if (err)
          goto FREE_RETURN;

        Sum_coupon_paid += coupon;

        if (Sum_coupon_paid > cpd->call[OTC_params->OTC].orig_barrier) {
          KnockedOut[realisation] = 1;
          Ncoupons--; /* To calculate the proba not to have hitten the TARN
                         Level */
        } else {
        }
      }
    }
  }

  /* =============================
  memory allocation
  ============================= */
  fund_leg = calloc(OTC_params->COPULAnSimul, sizeof(double));
  pd_leg = calloc(OTC_params->COPULAnSimul, sizeof(double));
  fund_legShort = calloc(OTC_params->COPULAnSimul, sizeof(double));
  pd_legShort = calloc(OTC_params->COPULAnSimul, sizeof(double));
  fundDf = calloc(cpd->fund_leg->num_cpn, sizeof(double));
  PdDfDom = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  FwdDfDom = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  FwdDfFor = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  strikeCap = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));
  strikeFloor = calloc(cpd->pd_leg->num_cpn + 1, sizeof(double));

  if (!fund_leg || !pd_leg || !fund_legShort || !pd_legShort || !fundDf ||
      !PdDfDom || !FwdDfDom || !FwdDfFor || !strikeCap || !strikeFloor) {
    err = "OTC_params->OTCcaller: memory allocation error";
    goto FREE_RETURN;
  }

  /* =============================
  Domestic/foreign
  ============================= */
  if (cpd->fund_leg->dom_for == 0) {
    fund_yc = (char *)(und->dom_yc);
    fund_mult = 1.0;
  } else {
    fund_yc = (char *)(und->for_yc);
    fund_mult = und->spot_fx;
  }

  /* =============================
  Precalculations
  ============================= */

  pd_otc_start = cpd->call[OTC_params->OTC].pd_idx;
  fund_otc_start = cpd->call[OTC_params->OTC].fund_idx;

  if (OTC_params->OTC < cpd->num_calls - 1) {
    pd_otc_short_end = cpd->call[OTC_params->OTC + 1].pd_idx;
    fund_otc_short_end = cpd->call[OTC_params->OTC + 1].fund_idx;
  } else {
    pd_otc_short_end = cpd->pd_leg->num_cpn;
    fund_otc_short_end = cpd->fund_leg->num_cpn;
  }

  dfDom = swp_f_df(und->today, cpd->call[OTC_params->OTC].ex_date, und->dom_yc);
  dfFor = swp_f_df(und->today, cpd->call[OTC_params->OTC].ex_date, und->for_yc);
  dfFund = swp_f_df(und->today, cpd->call[OTC_params->OTC].ex_date, fund_yc);

  dfFirstFunding =
      swp_f_df(und->today, cpd->fund_leg->cpn[fund_otc_start].start_date,
               fund_yc) /
      dfFund;
  dfNotExchgePd =
      swp_f_df(und->today, cpd->call[OTC_params->OTC].set_date, und->dom_yc) /
      dfDom;
  dfNotExchgeFund =
      swp_f_df(und->today, cpd->call[OTC_params->OTC].set_date, fund_yc) /
      dfFund;

  if (OTC_params->OTC < cpd->num_calls - 1) {
    dfNotExchgeEndPd =
        swp_f_df(und->today, cpd->call[OTC_params->OTC + 1].set_date,
                 und->dom_yc) /
        dfDom;
    dfNotExchgeEndFund =
        swp_f_df(und->today, cpd->call[OTC_params->OTC + 1].set_date, fund_yc) /
        dfFund;
  }

  if (und->today < start_date) {
    dfStartNotExchgePd = swp_f_df(und->today, start_date, und->dom_yc) / dfDom;
    dfStartNotExchgeFund = swp_f_df(und->today, start_date, fund_yc) / dfFund;
  } else {
    dfStartNotExchgePd = 1.0;
    dfStartNotExchgeFund = 1.0;
  }

  if (OTC_params->OTC < cpd->num_calls - 1) {
    dfNotExchgeShortPd =
        swp_f_df(und->today, cpd->call[OTC_params->OTC + 1].set_date,
                 und->dom_yc) /
        dfDom;
    dfNotExchgeShortFund =
        swp_f_df(und->today, cpd->call[OTC_params->OTC + 1].set_date, fund_yc) /
        dfFund;
  }

  for (i = fund_otc_start; i < cpd->fund_leg->num_cpn; i++)
    fundDf[i] =
        swp_f_df(und->today, cpd->fund_leg->cpn[i].pay_date, fund_yc) / dfFund;

  for (i = pd_otc_start; i < cpd->pd_leg->num_cpn; i++) {
    PdDfDom[i] =
        swp_f_df(und->today, cpd->pd_leg->cpn[i].pay_date, und->dom_yc) / dfDom;
    FwdDfDom[i] =
        swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->dom_yc) /
        dfDom;
    FwdDfFor[i] =
        swp_f_df(und->today, cpd->pd_leg->cpn[i].fx_val_date, und->for_yc) /
        dfFor;
    if (fabs(cpd->pd_leg->cpn[i].beta) > 1.0e-16) {
      strikeCap[i] = (cpd->pd_leg->cpn[i].cap - cpd->pd_leg->cpn[i].alpha) /
                     cpd->pd_leg->cpn[i].beta;
      strikeFloor[i] = (cpd->pd_leg->cpn[i].floor - cpd->pd_leg->cpn[i].alpha) /
                       cpd->pd_leg->cpn[i].beta;
    }
  }

  PdDfDom[cpd->pd_leg->num_cpn] =
      swp_f_df(und->today, cpd->pd_leg->not_ref.pay_date, und->dom_yc) / dfDom;
  FwdDfDom[cpd->pd_leg->num_cpn] =
      swp_f_df(und->today, cpd->pd_leg->not_ref.fx_val_date, und->dom_yc) /
      dfDom;
  FwdDfFor[cpd->pd_leg->num_cpn] =
      swp_f_df(und->today, cpd->pd_leg->not_ref.fx_val_date, und->for_yc) /
      dfFor;
  if (fabs(cpd->pd_leg->not_ref.beta) > 1.0e-16) {
    strikeCap[cpd->pd_leg->num_cpn] =
        (cpd->pd_leg->not_ref.cap - cpd->pd_leg->not_ref.alpha) /
        cpd->pd_leg->not_ref.beta;
    strikeFloor[cpd->pd_leg->num_cpn] =
        (cpd->pd_leg->not_ref.floor - cpd->pd_leg->not_ref.alpha) /
        cpd->pd_leg->not_ref.beta;
  }
  /* ==============
  Strike and DF for the Not Exch
  =============== */
  for (i = 0; i < 2; i++) {
    if (OTC_params->OTC + i < cpd->num_calls) {
      if (cpd->call[OTC_params->OTC + i].use_opt_str) {
        dfNotEx[i] =
            (swp_f_df(und->today, cpd->call[OTC_params->OTC + i].fx_val_date,
                      und->for_yc) /
             dfFor) /
            (swp_f_df(und->today, cpd->call[OTC_params->OTC + i].fx_val_date,
                      und->dom_yc) /
             dfDom);
      }
    }
  }

  if (!OTC_params->OTC && OTC_params->KO_do_optim) {
    KO_fund = calloc(OTC_params->COPULAnSimul, sizeof(double));
    KO_pd = calloc(OTC_params->COPULAnSimul, sizeof(double));

    if (!KO_fund || !KO_pd) {
      err = "KOPDpayoff: memory allocation failure";
      if (err)
        goto FREE_RETURN;
    }
  }
  /* ==========================================================
  =============================================================
  /
  /		FEE
  /
  =============================================================
  ========================================================== */
  fund_fee = 0.0;
  pd_fee = 0.0;
  fund_fee_mc = 0.0;
  pd_fee_mc = 0.0;

  ko_fund_fee = 0.0;
  ko_pd_fee = 0.0;
  ko_fund_fee_mc = 0.0;
  ko_pd_fee_mc = 0.0;

  /*Coupon Fee */

  /* ==================================
  1)	Funding Fee
  ================================== */
  if (start_date > cpd->call[OTC_params->OTC].ex_date)
    fund_fee -=
        swp_f_df(und->today, start_date, fund_yc) * cpd->fund_leg->notional;

  /* B(t  ,Ts) */
  fund_fee += swp_f_df(und->today,
                       cpd->fund_leg->cpn[fund_otc_start].start_date, fund_yc) *
              cpd->fund_leg->notional;

  for (i = fund_otc_start; i < fund_otc_short_end; i++) {
    /* ==================================
    coupon: spread + margin
    ================================== */
    coupon = swp_f_df(und->today, cpd->fund_leg->cpn[i].pay_date, fund_yc);

    fund_fee += coupon * cpd->fund_leg->cpn[i].cpn;

    if (OTC_params->OTC < cpd->num_calls - 1 && i == fund_otc_short_end - 1)
      fund_fee -= coupon * cpd->fund_leg->notional;
  }

  /*	In order to get the funding PV in domestic units */
  fund_fee *= fund_mult;

  /* ==================================
  2)	Power dual Fee
  ================================== */
  if (start_date > cpd->call[OTC_params->OTC].ex_date)
    pd_fee -= swp_f_df(und->today, start_date, und->dom_yc) * pd_not;

  for (i = pd_otc_start; i < pd_otc_short_end; i++) {
    cpn = cpd->pd_leg->cpn + i;

    /* ==================================
    Forward @ Tval with Convexity adjustment in the 3F
    ================================== */
    forward = OTC_precalc->FeeFxAdj[i] * und->spot_fx *
              swp_f_df(und->today, cpn->fx_val_date, und->for_yc) /
              swp_f_df(und->today, cpn->fx_val_date, und->dom_yc);

    /* ==================================
    SABR parameters @ Tval
    ================================== */
    err = OTCgetSMILEparams(cpd, und, cpn->fx_fix_time, cpn->fx_val_time,
                            cpn->fx_val_date, smile_mkt, smile_params,
                            OTC_params->OTCsmileFee);

    if (err)
      goto FREE_RETURN;

    /* ==================================
    domestic DF @ Tpay
    ================================== */
    dfDom = swp_f_df(und->today, cpn->pay_date, und->dom_yc);

    if (!cpn->use_opt_str) {
      if (fabs(cpn->beta) > 1.0e-16) {
        /* ==================================
        Floor
        ================================== */
        err = OTCgetOptValue(cpn->floored, 0, cpn->floor, strikeFloor[i],
                             cpn->beta, forward, cpn->fx_fix_time, smile_params,
                             0, &floor);

        if (err)
          goto FREE_RETURN;

        /* ==================================
        Cap
        ================================== */
        err = OTCgetOptValue(cpn->capped, 1, cpn->cap, strikeCap[i], cpn->beta,
                             forward, cpn->fx_fix_time, smile_params, 0, &cap);

        if (err)
          goto FREE_RETURN;
      } else {
        forward = floor = cap = 0.0;
      }

      pd_fee += dfDom * (cpn->alpha + cpn->beta * forward +
                         fabs(cpn->beta) * (floor - cap));
    } else {
      opt_string = 0.0;

      for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
        if (fabs(cpn->weights[str_idx]) > 1e-16 * pd_not) {
          if (smile_params->sigma * sqrt(cpn->fx_fix_time) > 1e-16 &&
              cpn->strikes[str_idx] > 1e-16)
            type = 3;
          else
            type = 1;

          err = cpd_vol_get_price(type, forward, cpn->fx_fix_time,
                                  cpn->strikes[str_idx], smile_params,
                                  &OptionValue);
          if (err)
            goto FREE_RETURN;

          opt_string += cpn->weights[str_idx] * OptionValue;
        }
      }

      /*	Coupon pv */
      pd_fee += dfDom * (cpn->wcst + cpn->wspot * forward + opt_string);
    }
  }

  if (OTC_params->OTC == cpd->num_calls - 1) {
    /*	Final Notional if we value the last KO date*/
    cpn = &(cpd->pd_leg->not_ref);
    /*	Discount */
    dfDom = swp_f_df(und->today, cpn->pay_date, und->dom_yc);

    /* ==================================
    Forward @ Tval
    ================================== */
    forward = OTC_precalc->NotFeeFxAdj[cpd->num_calls] * und->spot_fx *
              swp_f_df(und->today, cpn->fx_val_date, und->for_yc) /
              swp_f_df(und->today, cpn->fx_val_date, und->dom_yc);

    /* ==================================
    SABR parameters @ Tval
    ================================== */
    err = OTCgetSMILEparams(cpd, und, cpn->fx_fix_time, cpn->fx_val_time,
                            cpn->fx_val_date, smile_mkt, smile_params,
                            OTC_params->OTCsmileFee);

    if (err)
      goto FREE_RETURN;

    if (!cpn->use_opt_str) {
      if (fabs(cpn->beta) > 1.0e-16) {
        /* ==================================
        Floor
        ================================== */
        err = OTCgetOptValue(
            cpn->floored, 0, cpn->floor, strikeFloor[cpd->pd_leg->num_cpn],
            cpn->beta, forward, cpn->fx_fix_time, smile_params, 0, &floor);

        if (err)
          goto FREE_RETURN;

        /* ==================================
        Cap
        ================================== */
        err = OTCgetOptValue(cpn->capped, 1, cpn->cap,
                             strikeCap[cpd->pd_leg->num_cpn], cpn->beta,
                             forward, cpn->fx_fix_time, smile_params, 0, &cap);

        if (err)
          goto FREE_RETURN;

      } else {
        forward = floor = cap = 0.0;
      }

      /*	Coupon pv */
      pd_fee += dfDom * (cpn->alpha + cpn->beta * forward +
                         fabs(cpn->beta) * (floor - cap));
    } else {
      opt_string = 0.0;

      for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
        if (fabs(cpn->weights[str_idx]) > 1e-16 * pd_not) {
          if (smile_params->sigma * sqrt(cpn->fx_fix_time) > 1e-16 &&
              cpn->strikes[str_idx] > 1e-16)
            type = 3;
          else
            type = 1;

          err = cpd_vol_get_price(type, forward, cpn->fx_fix_time,
                                  cpn->strikes[str_idx], smile_params,
                                  &OptionValue);
          if (err)
            goto FREE_RETURN;

          opt_string += cpn->weights[str_idx] * OptionValue;
        }
      }

      /*	Coupon pv */
      pd_fee += dfDom * (cpn->wcst + cpn->wspot * forward + opt_string);
    }
  }

  /*Ko Fee */

  // Exchange the notionals if the deal has started
  for (i = 0; i < 2; i++) {
    if (OTC_params->OTC < cpd->num_calls - 1 || !i) {
      if (cpd->call[OTC_params->OTC + i].ex_date >
          start_date) //&& OTC_params->OTC == 0)
      {
        exchange_fund =
            fund_mult *
            swp_f_df(und->today, cpd->call[OTC_params->OTC + i].set_date,
                     fund_yc) *
            cpd->fund_leg->notional;

        if (!OTC_params->OTC && !i)
          ko_fund_fee = exchange_fund;

        if (OTC_params->KO_do_optim ||
            (!OTC_params->KO_do_optim && !(!OTC_params->OTC && !i))) {
          if (!i)
            fund_fee -= exchange_fund;
          else
            fund_fee += exchange_fund;
        }
      } else if (cpd->call[OTC_params->OTC + i].ex_date <= start_date &&
                 OTC_params->OTC == 0)
        ko_fund_fee = 0.0;

      dfDom = swp_f_df(und->today, cpd->call[OTC_params->OTC + i].set_date,
                       und->dom_yc);

      if (cpd->call[OTC_params->OTC + i].use_opt_str) {
        forward =
            OTC_precalc->NotFeeFxAdj[OTC_params->OTC + i] * und->spot_fx *
            swp_f_df(und->today, cpd->call[OTC_params->OTC + i].fx_val_date,
                     und->for_yc) /
            swp_f_df(und->today, cpd->call[OTC_params->OTC + i].fx_val_date,
                     und->dom_yc);

        /* ==================================
        SABR parameters @ Tval
        ================================== */
        err =
            OTCgetSMILEparams(cpd, und, cpd->call[OTC_params->OTC + i].ex_time,
                              cpd->call[OTC_params->OTC + i].fx_val_time,
                              cpd->call[OTC_params->OTC + i].fx_val_date,
                              smile_mkt, smile_params, OTC_params->OTCsmileFee);

        if (err)
          goto FREE_RETURN;

        opt_string = 0.0;

        for (str_idx = 0; str_idx < cpd->call[OTC_params->OTC + i].nstrikes;
             str_idx++) {
          if (fabs(cpd->call[OTC_params->OTC + i].weights[str_idx]) >
              1e-16 * pd_not) {
            if (smile_params->sigma *
                        sqrt(cpd->call[OTC_params->OTC + i].fx_fix_time) >
                    1e-16 &&
                cpd->call[OTC_params->OTC + i].strikes[str_idx] > 1e-16)
              type = 3;
            else
              type = 1;

            err = cpd_vol_get_price(
                type, forward, cpd->call[OTC_params->OTC + i].fx_fix_time,
                cpd->call[OTC_params->OTC + i].strikes[str_idx], smile_params,
                &OptionValue);
            if (err)
              goto FREE_RETURN;

            opt_string +=
                cpd->call[OTC_params->OTC + i].weights[str_idx] * OptionValue;
          }
        }

        /*	Coupon pv */
        exchange_pd =
            dfDom * (cpd->call[OTC_params->OTC + i].wcst +
                     cpd->call[OTC_params->OTC + i].wspot * forward +
                     opt_string + cpd->call[OTC_params->OTC + i].orig_fee);

        if (cpd->call[OTC_params->OTC + i].ex_date >
            start_date) //&& OTC_params->OTC == 0)
        {
          if (!OTC_params->OTC && !i)
            ko_pd_fee += exchange_pd;

          if (OTC_params->KO_do_optim ||
              (!OTC_params->KO_do_optim && !(!OTC_params->OTC && !i))) {
            if (!i)
              pd_fee -= exchange_pd;
            else
              pd_fee += exchange_pd;
          }
        } else if (cpd->call[OTC_params->OTC + i].ex_date <= start_date &&
                   OTC_params->OTC == 0)
          ko_pd_fee =
              exchange_pd - dfDom * cpd->call[OTC_params->OTC + i].pd_not_amt;

      } else {
        if (cpd->call[OTC_params->OTC + i].ex_date >
            start_date) //&& OTC_params->OTC == 0)
        {
          exchange_pd =
              dfDom * (pd_not + cpd->call[OTC_params->OTC + i].orig_fee);

          if (!OTC_params->OTC && !i)
            ko_pd_fee += exchange_pd;

          if (OTC_params->KO_do_optim ||
              (!OTC_params->KO_do_optim && !(!OTC_params->OTC && !i))) {
            if (!i)
              pd_fee -= exchange_pd;
            else
              pd_fee += exchange_pd;
          }
        } else if (cpd->call[OTC_params->OTC + i].ex_date <= start_date &&
                   OTC_params->OTC == 0)
          ko_pd_fee += dfDom * cpd->call[OTC_params->OTC + i].orig_fee;
      }
    }
  }

  /* ==========================================================
  =============================================================
  /
  /		PAYOFF
  /
  =============================================================
  ========================================================== */
  err = OTCdfPhi(und->lda_dom, cpd->call[OTC_params->OTC].ex_time,
                 und->sigma_time_rates, und->sigma_dom, und->sigma_n_rates,
                 &phiDom);

  if (err)
    goto FREE_RETURN;

  err = OTCdfPhi(und->lda_for, cpd->call[OTC_params->OTC].ex_time,
                 und->sigma_time_rates, und->sigma_for, und->sigma_n_rates,
                 &phiFor);

  if (err)
    goto FREE_RETURN;

  // Coupon evaluation
  smile_params->Use_BetaQuick_in_MC = 1;
  for (realisation = 0; realisation < OTC_params->COPULAnSimul &&
                        (OTC_params->KO_do_optim ||
                         fmc_current_std > OTC_params->KO_FMC_precision);
       realisation++) {
    coupon_pd = 0.0;
    coupon_fund = 0.0;
    ko_coupon_pd = 0.0;
    ko_coupon_fund = 0.0;

    // Evaluate the coupons until next KO date

    /* ==================================
    0)	Notional Exchange if we have a KO date before the start of the deal
    ================================== */
    if (start_date > cpd->call[OTC_params->OTC].ex_date) {
      start_time = (start_date - und->today) * YEARS_IN_DAY;

      err = OTCdf(OTC_params->KO_RatesMatrix[realisation][0], und->lda_dom,
                  cpd->call[OTC_params->OTC].ex_time, start_time, tstarTime,
                  phiDom, &df);

      if (err)
        goto FREE_RETURN;

      coupon_pd -= dfStartNotExchgePd * df * pd_not;
    }

    for (i = pd_otc_start; i < pd_otc_short_end; i++) {
      cpn = cpd->pd_leg->cpn + i;

      /* ==================================
      Forward @ Tval with Convexity adjustment
      ================================== */
      df = exp(OTC_precalc->df_const_for_fx_val[OTC_params->OTC][i] +
               OTC_precalc->df_lin_for_fx_val[OTC_params->OTC][i] *
                   OTC_params->KO_RatesMatrix[realisation][1]);

      dfFor = FwdDfFor[i] * df;

      df = exp(OTC_precalc->df_const_dom_fx_val[OTC_params->OTC][i] +
               OTC_precalc->df_lin_dom_fx_val[OTC_params->OTC][i] *
                   OTC_params->KO_RatesMatrix[realisation][0]);

      dfDom = FwdDfDom[i] * df;

      forward = OTC_precalc->FxAdj[OTC_params->OTC][i] *
                OTC_params->KO_FxMatrix[realisation][OTC_params->OTC] * dfFor /
                dfDom;

      /* ==================================
      domestic DF @ Tpay
      ================================== */
      if (cpn->pay_date != cpn->fx_val_date) {
        df = exp(OTC_precalc->df_const_pd_pay[OTC_params->OTC][i] +
                 OTC_precalc->df_lin_pd_pay[OTC_params->OTC][i] *
                     OTC_params->KO_RatesMatrix[realisation][0]);
      }

      dfDom = PdDfDom[i] * df;

      err = cpd_vol_get_smile_params_from_otc(smile_mkt->smile_spec_type + 1,
                                              OTC_params->OTC, i, OTC_precalc,
                                              smile_params);
      if (err)
        goto FREE_RETURN;

      if (!cpn->use_opt_str) {
        if (fabs(cpn->beta) > 1.0e-16) {
          /* ==================================
          Floor
          ================================== */
          err = OTCgetOptValue(
              cpn->floored, 0, cpn->floor, strikeFloor[i], cpn->beta, forward,
              cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
              smile_params, 1, &floor);

          if (err)
            goto FREE_RETURN;

          /* ==================================
          Cap
          ================================== */
          err = OTCgetOptValue(
              cpn->capped, 1, cpn->cap, strikeCap[i], cpn->beta, forward,
              cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
              smile_params, 1, &cap);

          if (err)
            goto FREE_RETURN;
        } else {
          forward = floor = cap = 0.0;
        }

        coupon_pd += dfDom * (cpn->alpha + cpn->beta * forward +
                              fabs(cpn->beta) * (floor - cap));
      } else {
        opt_string = 0.0;

        for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
          if (fabs(cpn->weights[str_idx]) > 1e-16 * pd_not) {
            if (OTC_precalc->pd_fwdsmile_vol[OTC_params->OTC][i] *
                        sqrt(cpn->fx_fix_time) >
                    1e-16 &&
                cpn->strikes[str_idx] > 1e-16)
              type = 3;
            else
              type = 1;

            err = cpd_vol_get_price(
                type, forward,
                cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
                cpn->strikes[str_idx], smile_params, &OptionValue);
            if (err)
              goto FREE_RETURN;

            opt_string += cpn->weights[str_idx] * OptionValue;
          }
        }

        /*	Coupon pv */
        coupon_pd += dfDom * (cpn->wcst + cpn->wspot * forward + opt_string);
      }
    }

    // Last Notional if we value the last KO date

    if (OTC_params->OTC == cpd->num_calls - 1) {
      cpn = &(cpd->pd_leg->not_ref);

      /*	Discount */
      df = exp(
          OTC_precalc->df_const_pd_pay[OTC_params->OTC][cpd->pd_leg->num_cpn] +
          OTC_precalc->df_lin_pd_pay[OTC_params->OTC][cpd->pd_leg->num_cpn] *
              OTC_params->KO_RatesMatrix[realisation][0]);

      discount = PdDfDom[cpd->pd_leg->num_cpn] * df;

      /*	Fwd fx */

      /* ==================================
      Forward @ Tval
      ================================== */

      df = exp(
          OTC_precalc
              ->df_const_dom_fx_val[OTC_params->OTC][cpd->pd_leg->num_cpn] +
          OTC_precalc
                  ->df_lin_dom_fx_val[OTC_params->OTC][cpd->pd_leg->num_cpn] *
              OTC_params->KO_RatesMatrix[realisation][0]);

      dfDom = FwdDfDom[cpd->pd_leg->num_cpn] * df;

      df = exp(
          OTC_precalc
              ->df_const_for_fx_val[OTC_params->OTC][cpd->pd_leg->num_cpn] +
          OTC_precalc
                  ->df_lin_for_fx_val[OTC_params->OTC][cpd->pd_leg->num_cpn] *
              OTC_params->KO_RatesMatrix[realisation][1]);

      dfFor = FwdDfFor[cpd->pd_leg->num_cpn] * df;

      forward = OTC_precalc->NotFxAdj[OTC_params->OTC][2] *
                OTC_params->KO_FxMatrix[realisation][OTC_params->OTC] * dfFor /
                dfDom;

      err = cpd_vol_get_smile_params_from_otc(
          smile_mkt->smile_spec_type + 1, OTC_params->OTC, cpd->pd_leg->num_cpn,
          OTC_precalc, smile_params);
      if (err)
        goto FREE_RETURN;

      if (!cpn->use_opt_str) {
        if (fabs(cpn->beta) > 1.0e-16) {
          /* ==================================
          Floor
          ================================== */
          err = OTCgetOptValue(
              cpn->floored, 0, cpn->floor, strikeFloor[cpd->pd_leg->num_cpn],
              cpn->beta, forward,
              cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
              smile_params, 1, &floor);

          if (err)
            goto FREE_RETURN;

          /* ==================================
          Cap
          ================================== */
          err = OTCgetOptValue(
              cpn->capped, 1, cpn->cap, strikeCap[cpd->pd_leg->num_cpn],
              cpn->beta, forward,
              cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
              smile_params, 1, &cap);

          if (err)
            goto FREE_RETURN;

        } else {
          forward = floor = cap = 0.0;
        }

        /*	Coupon pv */
        coupon_pd += discount * (cpn->alpha + cpn->beta * forward +
                                 fabs(cpn->beta) * (floor - cap));

      } else {
        opt_string = 0.0;

        for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
          if (fabs(cpn->weights[str_idx]) > 1e-16 * pd_not) {
            if (OTC_precalc->pd_fwdsmile_vol[OTC_params->OTC][i] *
                        sqrt(cpn->fx_fix_time) >
                    1e-16 &&
                cpn->strikes[str_idx] > 1e-16)
              type = 3;
            else
              type = 1;

            err = cpd_vol_get_price(
                type, forward,
                cpn->fx_fix_time - cpd->call[OTC_params->OTC].ex_time,
                cpn->strikes[str_idx], smile_params, &OptionValue);
            if (err)
              goto FREE_RETURN;

            opt_string += cpn->weights[str_idx] * OptionValue;
          }
        }

        /*	Coupon pv */
        coupon_pd += discount * (cpn->wcst + cpn->wspot * forward + opt_string);
      }
    }

    // Notional exchange + basis spreads
    if (start_date > cpd->call[OTC_params->OTC].ex_date) {
      start_time = (start_date - und->today) * YEARS_IN_DAY;

      err = OTCdf(cpd->fund_leg->dom_for
                      ? OTC_params->KO_RatesMatrix[realisation][1]
                      : OTC_params->KO_RatesMatrix[realisation][0],
                  cpd->fund_leg->dom_for ? und->lda_for : und->lda_dom,
                  cpd->call[OTC_params->OTC].ex_time, start_time, tstarTime,
                  cpd->fund_leg->dom_for ? phiFor : phiDom, &df);

      if (err)
        goto FREE_RETURN;

      coupon_fund -= dfStartNotExchgeFund * df * cpd->fund_leg->notional;
    }

    if (cpd->fund_leg->dom_for)
      df = exp(OTC_precalc->df_const_fund_start[OTC_params->OTC] +
               OTC_precalc->df_lin_fund_start[OTC_params->OTC] *
                   OTC_params->KO_RatesMatrix[realisation][1]);
    else
      df = exp(OTC_precalc->df_const_fund_start[OTC_params->OTC] +
               OTC_precalc->df_lin_fund_start[OTC_params->OTC] *
                   OTC_params->KO_RatesMatrix[realisation][0]);

    coupon_fund += dfFirstFunding * df * cpd->fund_leg->notional;

    for (i = fund_otc_start; i < fund_otc_short_end; i++) {
      /* ==================================
      coupon: spread + margin
      ================================== */
      if (cpd->fund_leg->dom_for)
        df = exp(OTC_precalc->df_const_fund_pay[OTC_params->OTC][i] +
                 OTC_precalc->df_lin_fund_pay[OTC_params->OTC][i] *
                     OTC_params->KO_RatesMatrix[realisation][1]);
      else
        df = exp(OTC_precalc->df_const_fund_pay[OTC_params->OTC][i] +
                 OTC_precalc->df_lin_fund_pay[OTC_params->OTC][i] *
                     OTC_params->KO_RatesMatrix[realisation][0]);

      /*  B(t  ,Te) */
      coupon_fund += fundDf[i] * df * cpd->fund_leg->cpn[i].cpn;

      if (i == fund_otc_short_end - 1 && OTC_params->OTC < cpd->num_calls - 1)
        coupon_fund -= fundDf[i] * df * cpd->fund_leg->notional;
    }

    if (cpd->fund_leg->dom_for)
      coupon_fund *= OTC_params->KO_FxMatrix[realisation][OTC_params->OTC];

    // Exchange the notionals
    for (i = 0; i < 2; i++) {
      if (OTC_params->OTC < cpd->num_calls - 1 || !i) {
        // Funding
        if (cpd->call[OTC_params->OTC + i].ex_date > start_date) {
          /* ==================================
          0)	Notional exchange at Tset and last coupon date
          ================================== */
          if (cpd->fund_leg->dom_for)
            df = exp(
                OTC_precalc->df_const_fund_pay[OTC_params->OTC]
                                              [cpd->fund_leg->num_cpn + 1 + i] +
                OTC_precalc->df_lin_fund_pay[OTC_params->OTC]
                                            [cpd->fund_leg->num_cpn + 1 + i] *
                    OTC_params->KO_RatesMatrix[realisation][1]);
          else
            df = exp(
                OTC_precalc->df_const_fund_pay[OTC_params->OTC]
                                              [cpd->fund_leg->num_cpn + 1 + i] +
                OTC_precalc->df_lin_fund_pay[OTC_params->OTC]
                                            [cpd->fund_leg->num_cpn + 1 + i] *
                    OTC_params->KO_RatesMatrix[realisation][0]);

          if (!i)
            exchange_fund = dfNotExchgeFund * df * cpd->fund_leg->notional;
          else
            exchange_fund = dfNotExchgeEndFund * df * cpd->fund_leg->notional;

          if (cpd->fund_leg->dom_for)
            exchange_fund *=
                OTC_params->KO_FxMatrix[realisation][OTC_params->OTC];

          if (!OTC_params->OTC && !i)
            ko_coupon_fund = exchange_fund;

          if (OTC_params->KO_do_optim ||
              (!OTC_params->KO_do_optim && !(!OTC_params->OTC && !i))) {
            if (!i)
              coupon_fund -= exchange_fund;
            else
              coupon_fund += exchange_fund;
          }

        } else if (cpd->call[OTC_params->OTC].ex_date <= start_date &&
                   OTC_params->OTC == 0)
          ko_coupon_fund = 0.0;

        // Pd
        df = exp(OTC_precalc->df_const_pd_pay[OTC_params->OTC]
                                             [cpd->pd_leg->num_cpn + 1 + i] +
                 OTC_precalc->df_lin_pd_pay[OTC_params->OTC]
                                           [cpd->pd_leg->num_cpn + 1 + i] *
                     OTC_params->KO_RatesMatrix[realisation][0]);

        if (!i)
          dfDom = df * dfNotExchgePd;
        else
          dfDom = df * dfNotExchgeEndPd;

        err = cpd_vol_get_smile_params_from_otc(
            smile_mkt->smile_spec_type + 1, OTC_params->OTC,
            cpd->pd_leg->num_cpn + i + 1, OTC_precalc, smile_params);
        if (err)
          goto FREE_RETURN;

        if (cpd->call[OTC_params->OTC + i].use_opt_str) {
          /* ==================================
          Forward @ Tval
          ================================== */
          forward = OTC_precalc->NotFxAdj[OTC_params->OTC][i] *
                    OTC_params->KO_FxMatrix[realisation][OTC_params->OTC] *
                    dfNotEx[i];

          df = exp(
              OTC_precalc->df_const_dom_fx_val[OTC_params->OTC]
                                              [cpd->pd_leg->num_cpn + i + 1] +
              OTC_precalc->df_lin_dom_fx_val[OTC_params->OTC]
                                            [cpd->pd_leg->num_cpn + i + 1] *
                  OTC_params->KO_RatesMatrix[realisation][0]);

          forward /= df;

          df = exp(
              OTC_precalc->df_const_for_fx_val[OTC_params->OTC]
                                              [cpd->pd_leg->num_cpn + i + 1] +
              OTC_precalc->df_lin_for_fx_val[OTC_params->OTC]
                                            [cpd->pd_leg->num_cpn + i + 1] *
                  OTC_params->KO_RatesMatrix[realisation][1]);

          forward *= df;

          opt_string = 0.0;
          for (str_idx = 0; str_idx < cpd->call[OTC_params->OTC + i].nstrikes;
               str_idx++) {
            if (fabs(cpd->call[OTC_params->OTC + i].weights[str_idx]) >
                1e-16 * pd_not) {
              if (OTC_precalc->pd_fwdsmile_vol[OTC_params->OTC]
                                              [cpd->pd_leg->num_cpn + i + 1] *
                          sqrt(cpd->call[OTC_params->OTC + i].fx_fix_time -
                               cpd->call[OTC_params->OTC].ex_time) >
                      1e-16 &&
                  cpd->call[OTC_params->OTC + i].strikes[str_idx] > 1e-16)
                type = 3;
              else
                type = 1;

              err = cpd_vol_get_price(
                  type, forward,
                  cpd->call[OTC_params->OTC + i].fx_fix_time -
                      cpd->call[OTC_params->OTC].ex_time,
                  cpd->call[OTC_params->OTC + i].strikes[str_idx], smile_params,
                  &OptionValue);
              if (err)
                goto FREE_RETURN;

              opt_string +=
                  cpd->call[OTC_params->OTC + i].weights[str_idx] * OptionValue;
            }
          }

          exchange_pd =
              dfDom * (cpd->call[OTC_params->OTC + i].wcst +
                       cpd->call[OTC_params->OTC + i].wspot * forward +
                       opt_string + cpd->call[OTC_params->OTC + i].orig_fee);

          if (cpd->call[OTC_params->OTC + i].ex_date > start_date) {
            if (!OTC_params->OTC && !i)
              ko_coupon_pd = exchange_pd;

            if (OTC_params->KO_do_optim ||
                (!OTC_params->KO_do_optim && !(!OTC_params->OTC && !i))) {
              if (!i)
                coupon_pd -= exchange_pd;
              else
                coupon_pd += exchange_pd;
            }
          } else if (cpd->call[OTC_params->OTC].ex_date <= start_date &&
                     OTC_params->OTC == 0)
            ko_coupon_pd += exchange_pd;
        } else {
          if (cpd->call[OTC_params->OTC + i].ex_date > start_date) {
            /* ==================================
            0)	Notional exchange at Tset and last coupon date
            ================================== */
            exchange_pd =
                dfDom * (pd_not + cpd->call[OTC_params->OTC + i].orig_fee);

            if (!OTC_params->OTC && !i)
              ko_coupon_pd = exchange_pd;

            if (OTC_params->KO_do_optim ||
                (!OTC_params->KO_do_optim && !(!OTC_params->OTC && !i))) {
              if (!i)
                coupon_pd -= exchange_pd;
              else
                coupon_pd += exchange_pd;
            }

          } else if (cpd->call[OTC_params->OTC].ex_date <= start_date &&
                     OTC_params->OTC == 0)
            ko_coupon_pd += dfDom * cpd->call[OTC_params->OTC + i].orig_fee;
        }
      }
    }

    ko_fund_fee_mc += ko_coupon_fund;
    ko_pd_fee_mc += ko_coupon_pd;
    fund_fee_mc += coupon_fund;
    pd_fee_mc += coupon_pd;

    if (!OTC_params->KO_do_optim) {
      if (KnockedOut[realisation] == 2) {
        // The deal has been called in the past
        Payoff_fund[realisation] = 0.0;
        Payoff_pd[realisation] = 0.0;

        if (OTC_params->KO_FMC && (realisation + 1) % 500 == 0 &&
            realisation > OTC_params->KO_FMC_min_paths) {
          fmc_current_std = fmc_current_var - fmc_current_mean *
                                                  fmc_current_mean /
                                                  (realisation + 1);
          fmc_current_std /= realisation + 1;
          fmc_current_std = sqrt(fmc_current_std / (realisation + 1));
          fmc_current_std /= pd_not;
        }
      } else if (KnockedOut[realisation] == 0) {
        fmc_nb_coupons++;
        // We pay the coupon
        Payoff_fund[realisation] = coupon_fund;
        Payoff_pd[realisation] = coupon_pd;

        if (OTC_params->KO_FMC) {
          fmc_current_var +=
              (coupon_pd - coupon_fund) * (coupon_pd - coupon_fund);
          fmc_current_mean += coupon_pd - coupon_fund;

          if ((realisation + 1) % 500 == 0 &&
              realisation > OTC_params->KO_FMC_min_paths) {
            fmc_current_std = fmc_current_var - fmc_current_mean *
                                                    fmc_current_mean /
                                                    (realisation + 1);
            fmc_current_std /= (realisation + 1);
            fmc_current_std = sqrt(fmc_current_std / (realisation + 1));
            fmc_current_std /= pd_not;
          }
        }
      } else if (KnockedOut[realisation] == 1) {
        // We pay the notional exchange
        Payoff_fund[realisation] = ko_coupon_fund;
        Payoff_pd[realisation] = ko_coupon_pd;

        if (OTC_params->KO_FMC) {
          fmc_current_var +=
              (ko_coupon_pd - ko_coupon_fund) * (ko_coupon_pd - ko_coupon_fund);
          fmc_current_mean += ko_coupon_pd - ko_coupon_fund;

          if ((realisation + 1) % 500 == 0 &&
              realisation > OTC_params->KO_FMC_min_paths) {
            fmc_current_std = fmc_current_var - fmc_current_mean *
                                                    fmc_current_mean /
                                                    (realisation + 1);
            fmc_current_std /= (realisation + 1);
            fmc_current_std = sqrt(fmc_current_std / (realisation + 1));
            fmc_current_std /= pd_not;
          }
        }
      } else {
        err = "invalid KnockedOut parameters";
        if (err)
          goto FREE_RETURN;
      }
    } else {
      Payoff_fund[realisation] = coupon_fund;
      Payoff_pd[realisation] = coupon_pd;

      if (!OTC_params->OTC) {
        KO_fund[realisation] = ko_coupon_fund;
        KO_pd[realisation] = ko_coupon_pd;
      }
    }
  }

  if (!OTC_params->KO_do_optim && OTC_params->KO_FMC) {
    fmc_nb_realisations = realisation;
    smessage("Phase 4 -Payoff evaluation Finished with %.2d simulations",
             fmc_nb_realisations);
  } else {
    fmc_nb_realisations = OTC_params->COPULAnSimul;
  }

  df = swp_f_df(und->today, tstarDate, und->dom_yc);

  if (!OTC_params->KO_do_optim) {
    PriceStd =
        (Payoff_pd[0] - Payoff_fund[0]) * (Payoff_pd[0] - Payoff_fund[0]);

    for (realisation = 1; realisation < fmc_nb_realisations; realisation++) {
      Payoff_pd[0] += Payoff_pd[realisation];
      Payoff_fund[0] += Payoff_fund[realisation];
      PriceStd += (Payoff_pd[realisation] - Payoff_fund[realisation]) *
                  (Payoff_pd[realisation] - Payoff_fund[realisation]);
    }

    PriceStd -= (Payoff_fund[0] - Payoff_pd[0]) *
                (Payoff_fund[0] - Payoff_pd[0]) / fmc_nb_realisations;
    PriceStd /= fmc_nb_realisations;
    PriceStd = sqrt(PriceStd / fmc_nb_realisations);

    if (OTC_params->OTCsmileFee > -1) {
      Payoff_pd[0] -=
          (pd_fee_mc / fmc_nb_realisations - pd_fee / df) * fmc_nb_coupons;
      Payoff_fund[0] -=
          (fund_fee_mc / fmc_nb_realisations - fund_fee / df) * fmc_nb_coupons;
    }

    Payoff_pd[0] /= fmc_nb_realisations;
    Payoff_fund[0] /= fmc_nb_realisations;

    /* ============================
    Get the adjustment at Tex
    ============================ */
    Payoff_fund[0] *= df;
    Payoff_pd[0] *= df;

    if (!cpd->call[OTC_params->OTC].pay_rec) {
      Results[0][0] += Payoff_fund[0] - Payoff_pd[0];
      Results[1][0] += fund_fee - pd_fee;
    } else {
      Results[0][0] -= Payoff_fund[0] - Payoff_pd[0];
      Results[1][0] -= fund_fee - pd_fee;
    }

    Results[2 + OTC_params->OTC][0] = cpd->call[OTC_params->OTC].ex_date;
    Results[2 + OTC_params->OTC][1] =
        (double)Ncoupons / OTC_params->COPULAnSimul;

    if (!cpd->call[OTC_params->OTC].pay_rec)
      Results[2 + OTC_params->OTC][2] = Payoff_fund[0] - Payoff_pd[0];
    else
      Results[2 + OTC_params->OTC][2] = -Payoff_fund[0] + Payoff_pd[0];

    Results[2 + OTC_params->OTC][3] = PriceStd;
    if (OTC_params->OTCsmileFee > -1)
      Results[2 + OTC_params->OTC][4] =
          (fund_fee_mc / fmc_nb_realisations - fund_fee / df) -
          (pd_fee_mc / fmc_nb_realisations - pd_fee / df);

    Results[2 + OTC_params->OTC][5] = fund_fee;
    Results[2 + OTC_params->OTC][6] = pd_fee;
  } else {
    if (!OTC_params->OTC) {
      ko_pd_fee_mc = (ko_pd_fee_mc / OTC_params->COPULAnSimul - ko_pd_fee / df);
      ko_fund_fee_mc =
          (ko_fund_fee_mc / OTC_params->COPULAnSimul - ko_fund_fee / df);
      for (realisation = 0; realisation < OTC_params->COPULAnSimul;
           realisation++) {
        save_values[0][1][realisation] = KO_fund[realisation] - ko_fund_fee_mc;
        save_values[0][1][realisation] -= KO_pd[realisation] - ko_pd_fee_mc;
        save_values[0][1][realisation] *= df;
      }
    }

    pd_fee = (pd_fee_mc / OTC_params->COPULAnSimul - pd_fee / df);
    fund_fee = (fund_fee_mc / OTC_params->COPULAnSimul - fund_fee / df);

    for (realisation = 0; realisation < OTC_params->COPULAnSimul;
         realisation++) {
      save_values[OTC_params->OTC + 1][0][realisation] =
          OTC_params->KO_FxMatrix[realisation][OTC_params->OTC];

      if (!KnockedOut[realisation]) {
        save_values[OTC_params->OTC + 1][1][realisation] =
            Payoff_fund[realisation] - fund_fee;
        save_values[OTC_params->OTC + 1][1][realisation] -=
            Payoff_pd[realisation] - pd_fee;
        save_values[OTC_params->OTC + 1][1][realisation] *= df;
      }
    }
  }

FREE_RETURN:

  if (KnockedOut)
    free(KnockedOut);
  if (Payoff_fund)
    free(Payoff_fund);
  if (Payoff_pd)
    free(Payoff_pd);

  if (fund_leg)
    free(fund_leg);
  if (pd_leg)
    free(pd_leg);

  if (fund_legShort)
    free(fund_legShort);
  if (pd_legShort)
    free(pd_legShort);

  if (KO_fund)
    free(KO_fund);
  if (KO_pd)
    free(KO_pd);

  if (fundDf)
    free(fundDf);
  if (PdDfDom)
    free(PdDfDom);
  if (FwdDfDom)
    free(FwdDfDom);
  if (FwdDfFor)
    free(FwdDfFor);

  if (strikeCap)
    free(strikeCap);
  if (strikeFloor)
    free(strikeFloor);

  return err;
}

/* ====================================================
UTILS
==================================================== */
Err OTCgetOptValue(int IsOption, int IsCap, double cap, double Strike,
                   double Notional, double Forward, double Maturity,
                   SMILE_PARAMETERS smile_params, int IsBetaVol,
                   double *Value) {
  int type;
  // double	smile_std  , smile_half_std;

  Err err = NULL;

  if (IsOption) {
    if (fabs(Notional) > 1.0e-16 && cap > -1.0e-16) {
      if (Notional > 0.0) {
        if (Strike > 1.0e-16)
          type = 3; /*	Call */
        else
          type = 1; /*	Call IV */
      } else {
        if (Strike > 1.0e-16)
          type = 4; /*	Put */
        else
          type = 2; /*	Put IV */
      }
    } else {
      if (cap <= -1.0e-16) {
        err = "Coupon cap/floor must be positive";
        goto FREE_RETURN;
      } else
        type = 0; /*	No cap */
    }

    if (!IsCap) {
      if (type == 1 || type == 3)
        type++;
      else if (type == 2 || type == 4)
        type--;
    }

    /*
    if (type >= 3)
    {
            err = cpd_vol_get_vol(Forward  ,Maturity  ,Strike  , smile_params  ,
    SRT_LOGNORMAL  , &smile_std); if(err) goto FREE_RETURN;

            smile_std *= sqrt (Maturity);
            smile_half_std = 0.5 * smile_std;
    }

    if ((Forward < 0.0 && type >= 3) || (smile_std < 1.0e-6 && type >= 3))
            type -= 2;

    *Value = OPT_VAL_MACRO(	type  ,
                                                    Forward  ,
                                                    Strike  ,
                                                    smile_std  ,
                                                    smile_half_std);
                                                    */

    err =
        cpd_vol_get_price(type, Forward, Maturity, Strike, smile_params, Value);
    if (err)
      goto FREE_RETURN;

  } else {
    *Value = 0.0;
  }

FREE_RETURN:
  return err;
}

Err OTCgetCopula(double forward, double forwardFixTime, double tstarTime,
                 double *mergeTimes, int mergeNtimes, double *sigDom,
                 double lda_dom, double *sigFor, double lda_for, double *sigFx,
                 double *corDF, double *corDFx, double *corFFx, double smilevol,
                 double smilealpha, double smilebeta, double smilerho,
                 int nSimul, int nPoints, int nIter, int nStd, double std,
                 int smileModel, double **matrix) {
  double stdDom, cumulDomForCor, stdFor, cumulForFxCor, meanFor, cumulDomFxCor;

  double **gauss = NULL, **cumulative = NULL;

  /* For the SSL */
  double volUp, volDown, shift, error;

  /* For the BM */
  // double			BMsigma  , BMalpha  , BMbeta  , BMlambda;
  double BMpi, BMerror, BMfwd1, BMsig1, BMfwd2, BMsig2;

  double *dX = NULL, *dCumulative = NULL;
  int l, lNbPoints;

  Err err = NULL;

  /* ==================================
  gaussian simulation
  ================================== */
  gauss = dmatrix(0, nSimul - 1, 0, 3 - 1);
  cumulative = dmatrix(0, 2 - 1, 0, nPoints - 1);

  if (!gauss || !cumulative) {
    err = "otc_caller: memory allocation error";
    goto FREE_RETURN;
  }

  err = balsam_generation(nSimul, 3, gauss);

  if (err)
    goto FREE_RETURN;

  adjust_emp_covar(gauss, nSimul, 3);

  /* ==================================
  get the 3F Cumulative correlations at the call date
  ================================== */
  err = OTCgetMoments(forwardFixTime, lda_dom, lda_for, tstarTime, mergeTimes,
                      mergeNtimes, sigDom, sigFor, sigFx, corDF, corDFx, corFFx,
                      &stdDom, &meanFor, &stdFor, &cumulDomForCor,
                      &cumulDomFxCor, &cumulForFxCor);

  if (err)
    goto FREE_RETURN;

  err = OTCcorrelateVariables(stdDom, stdFor, cumulDomForCor, cumulDomFxCor,
                              cumulForFxCor, nSimul, gauss, matrix);

  if (err)
    goto FREE_RETURN;

  /* ==================================
  Calibrate the Fx smile
  ================================== */
  if (smileModel == 0) {
    err = CalibSSLtoSABR(forwardFixTime, /* Maturity in years */
                         forward,        /* Forward */
                         std,            /* # std to calibrate skew */
                         /* smile Parameters */
                         smilevol, smilealpha, smilebeta, smilerho,
                         /* Results */
                         &volUp, &volDown, &shift, &error,
                         /* Parameters */
                         nIter);

    if (err)
      goto FREE_RETURN;

    if (fabs(error) > 1.0) {
      err = "otc_caller: SSL calib failed";
      goto FREE_RETURN;
    }

    /* ==================================
    Get the Cumulative Fx distribution under Qtstar
    ==================================*/
    err = SSLCumulative(forward, forwardFixTime, volUp, volDown, shift, nPoints,
                        nStd, cumulative);

    if (err)
      goto FREE_RETURN;

  } else if (smileModel == 1) {
    BMpi = 0.4;

    err =
        BMMCalibOnSabrStates(forward, forwardFixTime, smilevol, smilebeta,
                             smilealpha, smilerho, BMpi, std, &BMfwd1, &BMfwd2,
                             &BMsig1, &BMsig2, &BMpi, SRT_LOGNORMAL, &BMerror);

    if (BMerror > 1e-3) {
      err = "OTCcalibSmileModel: BiBetaCalibOnSabr calibration failed";
    }

    if (err)
      goto FREE_RETURN;

    lNbPoints = (long)(nPoints);
    dX = calloc(lNbPoints, sizeof(double));
    dCumulative = calloc(lNbPoints, sizeof(double));

    err = copula_gaussian_get_BMM_linterp_cumulative(
        forwardFixTime, forward, BMfwd1, BMsig1, BMfwd2, BMsig2, smilebeta,
        BMpi, 0.0005, nStd, nStd, &lNbPoints, &dX, &dCumulative);

    /*
    err = BMCumulative(	forward  ,
                                            forwardFixTime  ,
                                            SABRbeta  ,
                                            BMfwd1  ,
                                            BMfwd2  ,
                                            BMsig1  ,
                                            BMsig2  ,
                                            BMpi  ,
                                            nPoints  ,
                                            nStd  ,
                                            cumulative);
    */
    nPoints = lNbPoints;

    if (err)
      goto FREE_RETURN;

    for (l = 0; l < lNbPoints; l++) {
      cumulative[0][l] = dX[l];
      cumulative[1][l] = dCumulative[l];
    }

  } else {
  }

  err = GenerateRealisations(0.0, meanFor, cumulative, nSimul, nPoints, 1, 2, 0,
                             matrix);

  if (err)
    goto FREE_RETURN;

FREE_RETURN:
  if (gauss)
    free_dmatrix(gauss, 0, nSimul - 1, 0, 3 - 1);
  if (cumulative)
    free_dmatrix(cumulative, 0, 2 - 1, 0, nPoints - 1);
  if (dX)
    free(dX);
  if (dCumulative)
    free(dCumulative);

  return err;
}

Err OTCutilsInterpolate(double X, double XX, double XXX, double F1, double F2,
                        double F3, double Fseek, int method, double *XX1,
                        double *XX2) {
  double a, b, b1, b2, c, delta;
  int Islin = 0;
  int WhichSeg;

  Err err = NULL;

  if (method == 2) {
    if (fabs(X - XX) < 1.0e-15 && fabs(XX - XXX) < 1.0e-15) {
      err = "interpolation impossible";
      goto FREE_RETURN;
    } else if (fabs(X - XXX) < 1.0e-15) // lin interp X XX
    {
      Islin = 1;
      WhichSeg = 1;
    } else if (fabs(X - XX) < 1.0e-15 ||
               fabs(XX - XXX) < 1.0e-15) // lin interp X XXX
    {
      Islin = 1;
      WhichSeg = 2;
    }

    if (!Islin) {
      b1 = (F2 - F3) * (X - XX) * (X + XX);
      b1 -= (F1 - F2) * (XX - XXX) * (XX + XXX);

      b2 = (XX - XXX) * (X - XX) * (X + XX);
      b2 -= (X - XX) * (XX - XXX) * (XX + XXX);

      b = b1 / b2;

      a = (F1 - F2) - b * (X - XX);
      a /= (X - XX) * (X + XX);

      c = F3 - a * XXX * XXX - b * XXX;

      delta = b * b - 4.0 * a * (c - Fseek);

      if (delta < 0.0) {
        if (X < XX) {
          Islin = 1;
          WhichSeg = 1;
        } else {
          Islin = 1;
          WhichSeg = 3;
        }
      } else {
        *XX1 = 0.5 * (-b + sqrt(delta)) / a;
        *XX2 = 0.5 * (-b - sqrt(delta)) / a;
      }
    }

    if (Islin) {
      if (WhichSeg == 1) {
        b = (F1 - F2) / (X - XX);
        c = F1 - X * (F1 - F2) / (X - XX);
      } else if (WhichSeg == 2) {
        b = (F1 - F3) / (X - XXX);
        c = F1 - X * (F1 - F3) / (X - XXX);
      } else if (WhichSeg == 3) {
        b = (F2 - F3) / (XX - XXX);
        c = F2 - XX * (F2 - F3) / (XX - XXX);
      }

      if (fabs(b) < 1.0e-16)
        *XX1 = F1;
      else
        *XX1 = (F1 - c) / b;

      *XX2 = *XX1;
    }
  } else if (method == 1) {
    b = (F1 - F2) / (X - XX);
    c = F1 - X * (F1 - F2) / (X - XX);

    if (fabs(b) < 1.0e-16)
      *XX1 = F1;
    else
      *XX1 = (F1 - c) / b;

    *XX2 = *XX1;
  } else {
    *XX1 = F1;
    *XX2 = *XX1;
  }

FREE_RETURN:
  return err;
}

Err OTCutilsGetIndex(double *probas, int Nprobas, double prob, int isLinear,
                     int *index1, int *index2, int *index3, int *method) {
  // Gets the indexes for increasing probabilities
  // double step;
  double test, leftDist, rightDist;
  int dir = 1, stop = 0;

  Err err = NULL;

  if (Nprobas < 2) {
    err = "dimension error in OTCutilsGetIndex";
    goto FREE_RETURN;
  }

  if (isLinear) {
    *index2 = 0;
    while (probas[*index2] < prob && *index2 < Nprobas)
      (*index2)++;

    if ((*index1) < 0 || (*index2) == Nprobas) {
      if ((*index1) < 0 || (*index2) == Nprobas) {
        *index1 = 0;
        *index2 = 0;
        *index3 = 0;
        *method = 0;
      } else {
        *index1 = Nprobas - 1;
        *index2 = Nprobas - 1;
        *index3 = Nprobas - 1;
        *method = 0;
      }
    } else {
      *index1 = *index2 - 1;
      *index3 = *index2;
      *method = 1;
    }
  } else {
    if (prob < probas[0]) {
      *method = 2;
      *index1 = 0;
    } else if (prob > probas[Nprobas - 1]) {
      *method = 2;
      *index1 = Nprobas - 3;
    } else {
      // step = 0.3 * probas[(int)(Nprobas / 2.0)] - probas[(int)(Nprobas / 2.0)
      // - 1];
      *index1 = (int)(prob * (double)(Nprobas - 1));

      if (*index1 >= Nprobas - 1)
        *index1 = (int)((double)Nprobas / 2.0);

      if (&index1 <= 0)
        *index1 = 1;

      test = probas[*index1];

      if (test > prob)
        dir = -1;

      do {
        *index1 += dir;
        test = probas[*index1];
        if ((dir == 1 && test > prob) || (dir == -1 && test < prob) ||
            *index1 <= 0 || *index1 >= Nprobas - 1)
          stop = 1;

      } while (stop == 0);

      if (*index1 <= 0 || *index1 >= Nprobas - 1) {
        *method = 2;
        if (*index1 >= Nprobas - 1)
          *index1 = Nprobas - 3;
        else
          *index1 = 0;
      } else {
        *method = 2;
        if (*index1 == 1)
          *index1 = 0;
        else if (*index1 == Nprobas - 2)
          *index1 = Nprobas - 3;
        else if (dir == 1) {
          leftDist = fabs(probas[*index1 - 2] - probas[*index1]);
          rightDist = fabs(probas[*index1] - probas[*index1 + 1]);
          if (leftDist > rightDist)
            *index1 = *index1 - 1;
          else
            *index1 = *index1 - 2;
        } else {
          leftDist = fabs(probas[*index1 - 1] - probas[*index1]);
          rightDist = fabs(probas[*index1] - probas[*index1 + 2]);
          if (leftDist < rightDist)
            *index1 = *index1 - 1;
          else
            *index1 = *index1;
        }
      }
    }

    *index2 = *index1 + 1;
    *index3 = *index1 + 2;
  }
FREE_RETURN:
  return err;
}

Err OTCutilsGetFxRealisations(otcpd_params *otc_params) {
  Err err = NULL;
  double test1, test2;
  int i, index1, index2, index3, method;

  for (i = 0; i < otc_params->COPULAnSimul; i++) {
    otc_params->KO_FxMatrix[i][otc_params->OTC] =
        norm(otc_params->KO_FxGauss[i][otc_params->OTC]);

    err = OTCutilsGetIndex(otc_params->CUMUL[1], otc_params->CUMULnPoints,
                           otc_params->KO_FxMatrix[i][otc_params->OTC],
                           otc_params->CUMULlinear, &index1, &index2, &index3,
                           &method);

    if (err)
      goto FREE_RETURN;

    err = OTCutilsInterpolate(
        otc_params->CUMUL[0][index1], otc_params->CUMUL[0][index2],
        otc_params->CUMUL[0][index3], otc_params->CUMUL[1][index1],
        otc_params->CUMUL[1][index2], otc_params->CUMUL[1][index3],
        otc_params->KO_FxMatrix[i][otc_params->OTC], method, &test1, &test2);

    if (err)
      goto FREE_RETURN;

    if (otc_params->CUMULlinear == 0) {
      if (test1 > otc_params->CUMUL[0][index1] &&
          test1 < otc_params->CUMUL[0][index3])
        otc_params->KO_FxMatrix[i][otc_params->OTC] = test1;
      else if (test2 > otc_params->CUMUL[0][index1] &&
               test2 < otc_params->CUMUL[0][index3])
        otc_params->KO_FxMatrix[i][otc_params->OTC] = test2;
      else if (index1 == 0 && fabs(otc_params->KO_FxMatrix[i][otc_params->OTC] -
                                   test1) < 1.0e-10)
        otc_params->KO_FxMatrix[i][otc_params->OTC] = test1;
      else if (otc_params->KO_FxMatrix[i][otc_params->OTC] >
               otc_params->CUMUL[1][otc_params->CUMULnPoints - 1])
        otc_params->KO_FxMatrix[i][otc_params->OTC] =
            otc_params->CUMUL[0][otc_params->CUMULnPoints - 1];
      else if (otc_params->KO_FxMatrix[i][otc_params->OTC] <
               otc_params->CUMUL[1][0])
        otc_params->KO_FxMatrix[i][otc_params->OTC] = otc_params->CUMUL[0][0];
      else {
        err = "no solution in interpolation";
        if (err)
          goto FREE_RETURN;
      }
    } else {
      otc_params->KO_FxMatrix[i][otc_params->OTC] = test1;
    }
  }

FREE_RETURN:
  return err;
}

Err OTCutilsCorrelateFxAndRates(otcpd_params *OTC_params, double meanDom,
                                double stdDom, double meanFor, double stdFor,
                                double corDF, double corDFx, double corFFx) {
  long i;
  Err err = NULL;

  double b, c, d, e, f, x, y, z;

  b = corDFx;
  c = 1 - b * b;
  c = sqrt(c);

  d = corFFx;

  if (c != 0) {
    e = (corDF - b * d) / c;
    f = 1 - d * d - e * e;
    if (f < 0) {
      return "Check your correlations";
    }
    f = sqrt(f);
  } else {
    if (fabs(corDFx - corFFx) < 1e-10) {
      e = 1 - d * d;
      f = 0;
      if (e < 0) {
        return "Check your correlations";
      }
      e = sqrt(e);
    } else {
      return "Check your correlations";
    }
  }

  for (i = 0; i < OTC_params->COPULAnSimul; i++) {
    x = OTC_params->KO_FxGauss[i][OTC_params->OTC];
    y = OTC_params->KO_RatesGauss[i][0];
    z = OTC_params->KO_RatesGauss[i][1];

    OTC_params->KO_RatesMatrix[i][0] = meanDom + stdDom * (b * x + c * y);
    OTC_params->KO_RatesMatrix[i][1] =
        meanFor + stdFor * (d * x + e * y + f * z);
  }

  return err;
}

/* ==================================
 Calibration using NAG
================================== */
Err alloc_calib_sl_params(int iNbStrikes, CALIBSL_PARAMS sParams) {
  Err err = NULL;

  sParams->iNbStrikes = iNbStrikes;

  sParams->dStrikes = calloc(sParams->iNbStrikes, sizeof(double));
  sParams->dTargetPrices = calloc(sParams->iNbStrikes, sizeof(double));
  sParams->dTargetVols = calloc(sParams->iNbStrikes, sizeof(double));

  if (!sParams->dStrikes || !sParams->dStrikes || !sParams->dStrikes) {
    err = "Memory allocation faillure in alloc_calib_sl_params";
    return err;
  }

  return err;
}

void free_calib_sl_params(CALIBSL_PARAMS sParams) {
  if (sParams) {
    if (sParams->dStrikes)
      free(sParams->dStrikes);
    if (sParams->dTargetPrices)
      free(sParams->dTargetPrices);
    if (sParams->dTargetVols)
      free(sParams->dTargetVols);
  }
}

static void NAG_CALL MixedSL_Calib_Fonctionnal(Integer n, double x[],
                                               double *objf, double g[],
                                               Nag_Comm *comm) {
  CALIBSL_PARAMS sParams = (CALIBSL_PARAMS)(comm->p);
  int i;
  double dPremium, dOption, dVol, dErrorTemp, dError, dRescale;
  double dConv;
  double stdev, d1, d2, dDelta1, dDelta2, dVega1, dVega2;
  Err err = NULL;

  dError = 0.0;

  if (sParams->dForward + x[0] > 0.0) {
    dConv = 1.0;
  } else {
    dConv = -1.0;
  }

  /* Initialisation */
  if (comm->flag == 2) {
    g[0] = 0.0;
    g[1] = 0.0;
    g[2] = 0.0;

    dRescale = 2.0 * 10000.0 * 10000.0;
  }

  for (i = 0; i < sParams->iNbStrikes; i++) {
    /* Pricing */
    dPremium =
        srt_f_optblksch(dConv * (sParams->dForward + x[0]),
                        dConv * (sParams->dStrikes[i] + x[0]), x[1] + x[2],
                        sParams->dMaturity, 1.0, SRT_CALL, SRT_PREMIUM);

    dOption = sParams->dProba1 * dPremium;

    dPremium =
        srt_f_optblksch(dConv * (sParams->dForward + x[0]),
                        dConv * (sParams->dStrikes[i] + x[0]), x[1] - x[2],
                        sParams->dMaturity, 1.0, SRT_CALL, SRT_PREMIUM);

    dOption += sParams->dProba2 * dPremium;

    if (dConv < 0.0) {
      dOption += sParams->dForward - sParams->dStrikes[i];
    }

    if (sParams->iOptimiseOnVol) {
      err = srt_f_optimpvol(dOption, sParams->dForward, sParams->dStrikes[i],
                            sParams->dMaturity, 1.0, SRT_CALL,
                            sParams->eVolType, &dVol);

      if (err) {
        dVol = 0.0;
        err = NULL;
      }

      dErrorTemp = dVol - sParams->dTargetVols[i];
    } else {
      dErrorTemp = dOption - sParams->dTargetPrices[i];
    }

    dError += pow(dErrorTemp * 10000, 2.0);

    if (comm->flag == 2) {
      /* Calculates the first derivatives */
      stdev = (x[1] + x[2]) * sParams->dSqMaturity;
      d1 = log((sParams->dForward + x[0]) / (sParams->dStrikes[i] + x[0])) /
               stdev +
           0.5 * stdev;
      d2 = d1 - stdev;

      dDelta1 = norm(d1) - norm(d2);
      dVega1 =
          dConv * (sParams->dForward + x[0]) * gauss(d1) * sParams->dSqMaturity;

      stdev = (x[1] - x[2]) * sParams->dSqMaturity;
      d1 = log((sParams->dForward + x[0]) / (sParams->dStrikes[i] + x[0])) /
               stdev +
           0.5 * stdev;
      d2 = d1 - stdev;

      dDelta2 = norm(d1) - norm(d2);
      dVega2 =
          dConv * (sParams->dForward + x[0]) * gauss(d1) * sParams->dSqMaturity;

      g[0] += dConv *
              (sParams->dProba1 * dDelta1 + sParams->dProba2 * dDelta2) *
              dErrorTemp;
      g[1] +=
          (sParams->dProba1 * dVega1 + sParams->dProba2 * dVega2) * dErrorTemp;
      g[2] +=
          (sParams->dProba1 * dVega1 - sParams->dProba2 * dVega2) * dErrorTemp;
    }
  }

  /* Last rescale */
  if (comm->flag == 2) {
    g[0] *= dRescale;
    g[1] *= dRescale;
    g[2] *= dRescale;
  }

  if (dError < sParams->dTolerance) {
    *objf = 0.0;
    g[0] = 0.0;
    g[1] = 0.0;
    g[2] = 0.0;
  } else {
    *objf = dError;
  }
}

static void NAG_CALL MixedSL_Calib_LeastSquare(Integer m, Integer n, double x[],
                                               double fvec[], double fjac[],
                                               Integer tdj, Nag_Comm *comm) {
  CALIBSL_PARAMS sParams = (CALIBSL_PARAMS)(comm->p);
  int i;
  double dPremium, dOption, dVol, dErrorTemp;
  double dConv;
  double stdev, d1, d2, dDelta1, dDelta2, dVega1, dVega2;
  Err err = NULL;

  if (sParams->dForward + x[0] > 0.0) {
    dConv = 1.0;
  } else {
    dConv = -1.0;
  }

  /* Initialisation */

  for (i = 0; i < sParams->iNbStrikes; i++) {
    /* Pricing */
    dPremium =
        srt_f_optblksch(dConv * (sParams->dForward + x[0]),
                        dConv * (sParams->dStrikes[i] + x[0]), x[1] + x[2],
                        sParams->dMaturity, 1.0, SRT_CALL, SRT_PREMIUM);

    dOption = sParams->dProba1 * dPremium;

    dPremium =
        srt_f_optblksch(dConv * (sParams->dForward + x[0]),
                        dConv * (sParams->dStrikes[i] + x[0]), x[1] - x[2],
                        sParams->dMaturity, 1.0, SRT_CALL, SRT_PREMIUM);

    dOption += sParams->dProba2 * dPremium;

    if (dConv < 0.0) {
      dOption += sParams->dForward - sParams->dStrikes[i];
    }

    if (sParams->iOptimiseOnVol) {
      err = srt_f_optimpvol(dOption, sParams->dForward, sParams->dStrikes[i],
                            sParams->dMaturity, 1.0, SRT_CALL,
                            sParams->eVolType, &dVol);

      if (err) {
        dVol = 0.0;
        err = NULL;
      }

      dErrorTemp = dVol - sParams->dTargetVols[i];
    } else {
      dErrorTemp = dOption - sParams->dTargetPrices[i];
    }

    if (!comm->flag != 1) {
      fvec[i] = dErrorTemp;
    }

    if (!comm->flag == 0) {
      /* Calculates the first derivatives */
      stdev = (x[1] + x[2]) * sParams->dSqMaturity;
      d1 = log((sParams->dForward + x[0]) / (sParams->dStrikes[i] + x[0])) /
               stdev +
           0.5 * stdev;
      d2 = d1 - stdev;

      dDelta1 = norm(d1) - norm(d2);
      dVega1 =
          dConv * (sParams->dForward + x[0]) * gauss(d1) * sParams->dSqMaturity;

      stdev = (x[1] - x[2]) * sParams->dSqMaturity;
      d1 = log((sParams->dForward + x[0]) / (sParams->dStrikes[i] + x[0])) /
               stdev +
           0.5 * stdev;
      d2 = d1 - stdev;

      dDelta2 = norm(d1) - norm(d2);
      dVega2 =
          dConv * (sParams->dForward + x[0]) * gauss(d1) * sParams->dSqMaturity;

      fjac[tdj * i] =
          dConv * (sParams->dProba1 * dDelta1 + sParams->dProba2 * dDelta2);
      fjac[tdj * i + 1] = sParams->dProba1 * dVega1 + sParams->dProba2 * dVega2;
      fjac[tdj * i + 2] = sParams->dProba1 * dVega1 - sParams->dProba2 * dVega2;
    }
  }
}

Err Calib_MixedSL_ToSmile(CALIBSL_PARAMS sParams, double *dShift,
                          double *dVolDown, double *dVolUp,
                          double *dFittingError) {
  Err err = NULL;

  double dInputs[3], dLBound[3], dHBound[3], dGradient[3], dFuncValue[3],
      dJacobi[9];

  static NagError fail;
  Nag_BoundType bound; // NAG: says there are constraints
  Nag_Comm nagcomm;
  Nag_E04_Opt options; // optional parameters of the NAG optimizer

  nag_opt_init(&options);

  bound = Nag_Bounds;
  fail.print = FALSE;
  options.print_level = Nag_NoPrint;
  options.list = FALSE;
  options.optim_tol = sParams->dForward / 1000.0;
  options.init_state = Nag_Init_None;
  /*
  strcpy(options.outfile  , "C:\\Temp\\OptimDLM.txt");
  */

  dInputs[0] = *dShift;
  dInputs[1] = 0.5 * (*dVolUp + *dVolDown);
  dInputs[2] = 0.5 * (*dVolUp - *dVolDown);

  dLBound[0] = -sParams->dForward * 100.0;
  dHBound[0] = sParams->dForward * 100.0;

  dLBound[1] = 0.00001;
  dHBound[1] = sParams->dForward * 10000.0;
  dLBound[2] = 0.00001;
  dHBound[2] = sParams->dForward * 10000.0;

  nagcomm.p = (void *)(sParams);

  if (0) {
    nag_opt_bounds_deriv(3, MixedSL_Calib_Fonctionnal, bound, dLBound, dHBound,
                         dInputs, dFittingError, dGradient, &options, &nagcomm,
                         &fail);
  } else {
    nag_opt_lsq_deriv(3, 3, MixedSL_Calib_LeastSquare, dInputs, dFittingError,
                      dFuncValue, dJacobi, 3, &options, &nagcomm, &fail);
  }

  *dFittingError = sqrt(*dFittingError / sParams->iNbStrikes);
  *dShift = dInputs[0];
  *dVolDown = dInputs[1] - dInputs[2];
  *dVolUp = dInputs[1] + dInputs[2];

  return err;
}

/* -----------------------------------------------------------------------------------------------

          Save matrices or doubles to analyze them

   ------------------------------------------------------------------------------------------------
 */

void OTCutilsSaveMatrix(int iNRows, int iNCols, double **Matrix) {
  int i, j;
  FILE *fid;

  fid = fopen("C:\\srt_matrix.txt", "a");

  for (i = 0; i < iNRows; i++) {
    for (j = 0; j < iNCols; j++) {
      fprintf(fid, "%f//", Matrix[i][j]);
    }
    fprintf(fid, "\n");
  }

  fclose(fid);
}

void OTCutilsSaveDouble(char *fileName, double toSave) {
  FILE *fid;
  char tempString[256];

  /*
  fid = fopen("C:\\SaveDouble.txt"  ,"ab");
  fwrite(&toSave  ,1  ,sizeof(double)  ,fid);
  fclose(fid);
  */

  sprintf(tempString, "C:\\%s.txt", fileName);

  fid = fopen(tempString, "a");
  fprintf(fid, "%f\n", toSave);
  fclose(fid);
}

Err Srt_smooth_matrix(int method, int num_adj, int NRow, int NCol,
                      double **in_matrix, double **out_matrix) {
  int i, j, ii, jj, posLeft, posRight, posAbove, posBelow, num_points,
      num_max_points;
  double coef;
  Err err = NULL;

  num_max_points = (1 + 2 * num_adj) * (1 + 2 * num_adj);

  for (i = 0; i < NRow; i++) {
    for (j = 0; j < NCol; j++) {
      posLeft = j - num_adj;
      posRight = j + num_adj;
      posAbove = i - num_adj;
      posBelow = i + num_adj;

      if (posLeft < 0)
        posLeft = 0;
      if (posRight > NCol - 1)
        posRight = NCol - 1;
      if (posAbove < 0)
        posAbove = 0;
      if (posBelow > NRow - 1)
        posBelow = NRow - 1;

      num_points = (1 + posRight - posLeft) * (1 + posBelow - posAbove);

      for (ii = posAbove; ii <= posBelow; ii++) {
        for (jj = posLeft; jj <= posRight; jj++) {
          if (method == 0) {
            coef = 1.0 / (double)num_points;
            out_matrix[ii][jj] += coef * in_matrix[i][j];
          } else if (method == 1) {
            if (ii == i && jj == j)
              coef = (double)(1 + num_max_points - num_points) /
                     (double)(num_max_points);
            else
              coef = 1.0 / (double)num_max_points;

            out_matrix[ii][jj] += coef * in_matrix[i][j];
          } else if (method == 2) {
            coef = 1.0 / (double)num_points;
            out_matrix[i][j] += coef * in_matrix[ii][jj];
          }
        }
      }
    }
  }

  return err;
}
