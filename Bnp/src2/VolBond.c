
#include "VolBond.h"
#include "VolBondCalib.h"
#include "VolBondProdStruct.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"

double CMSAdjustementWithDelay(double fact1, double fact2, long today,
                               long eventDate, swap_for_lgm2f compSwap,
                               double comp_sabrBetaVol, double comp_sabrBeta,
                               swap_for_lgm2f fwdSwap, double fwd_sabrBetaVol,
                               double fwd_sabrBeta, double swapCorrel,
                               TermStruct *l, char *yc) {
  double compSwapValue;
  double fwdSwapValue;
  double temp;

  compSwapValue =
      max(0.000001,
          SWAP_FOR_LGM2F_SwapRate_From_Factors_Value(fact1, fact2, &compSwap));

  fwdSwapValue =
      max(0.000001,
          SWAP_FOR_LGM2F_SwapRate_From_Factors_Value(fact1, fact2, &fwdSwap));

  temp = fwd_sabrBetaVol * pow(fwdSwapValue, fwd_sabrBeta - 1);

  temp = ((fwdSwap.SwapStartDate - compSwap.SwapStartDate) / 365.25) *
         comp_sabrBetaVol * pow(compSwapValue, comp_sabrBeta) * swapCorrel /
         (1 + ((fwdSwap.SwapStartDate - compSwap.SwapStartDate) / 365.25) *
                  compSwapValue);

  temp = fwd_sabrBetaVol * pow(fwdSwapValue, fwd_sabrBeta - 1) *
         (fwdSwap.fix_cvgs[0] * fwdSwap.fix_NDates * fwdSwapValue *
              pow(1 + fwdSwap.fix_cvgs[0] * fwdSwapValue,
                  -fwdSwap.fix_NDates - 1) /
              (1 -
               pow(1 + fwdSwap.fix_cvgs[0] * fwdSwapValue, -fwdSwap.fix_NDates))

          - 1);

  return fwd_sabrBetaVol * pow(fwdSwapValue, fwd_sabrBeta - 1) *
         (((fwdSwap.SwapStartDate - compSwap.SwapStartDate) / 365.25) *
              comp_sabrBetaVol * pow(compSwapValue, comp_sabrBeta) *
              swapCorrel /
              (1 + ((fwdSwap.SwapStartDate - compSwap.SwapStartDate) / 365.25) *
                       compSwapValue) -
          fwd_sabrBetaVol * pow(fwdSwapValue, fwd_sabrBeta - 1) *
              (fwdSwap.fix_cvgs[0] * fwdSwap.fix_NDates * fwdSwapValue *
                   pow(1 + fwdSwap.fix_cvgs[0] * fwdSwapValue,
                       -fwdSwap.fix_NDates - 1) /
                   (1 - pow(1 + fwdSwap.fix_cvgs[0] * fwdSwapValue,
                            -fwdSwap.fix_NDates))

               - 1));
}

double CMSAdjustementWithoutDelay(double fact1, double fact2, long today,
                                  long eventDate, swap_for_lgm2f fwdSwap,
                                  double fwd_sabrBetaVol, double fwd_sabrBeta,
                                  TermStruct *l, char *yc) {
  double fwdSwapValue;

  fwdSwapValue =
      max(0.000001,
          SWAP_FOR_LGM2F_SwapRate_From_Factors_Value(fact1, fact2, &fwdSwap));

  return fwd_sabrBetaVol * pow(fwdSwapValue, fwd_sabrBeta - 1) *
         (-fwd_sabrBetaVol * pow(fwdSwapValue, fwd_sabrBeta - 1) *
          (fwdSwap.fix_cvgs[0] * fwdSwap.fix_NDates * fwdSwapValue *
               pow(1 + fwdSwap.fix_cvgs[0] * fwdSwapValue,
                   -fwdSwap.fix_NDates - 1) /
               (1 - pow(1 + fwdSwap.fix_cvgs[0] * fwdSwapValue,
                        -fwdSwap.fix_NDates))

           - 1));
}

Err Price_VolBondCoupon(long today,
                        //	The underlying
                        char *yc,      //	yc
                        TermStruct *l, // underlying
                        //	SABR Parameters
                        double sabrBetaVol, double sabrAlpha, double sabrBeta,
                        double sabrRho,
                        // Complementary Swap Sabr parameters
                        double swapRho, double CompSabrBeta,
                        double CompSabrBetaVol, int NGaussQuadrature, double *x,
                        double *w, volBondCoupon *volBondCpn,
                        double *VolBondCouponPrice) {
  long i, j;
  Err err;
  double fact1, fact2;
  double CMSAdj;
  double Payoff;
  swap_for_lgm2f swapForCMSAdj;
  my_lgm2F lgm2fForMeanVarCov;
  double GQRho;
  long Test;

  err = Fill_MY_LGM2F(1, &(lgm2fForMeanVarCov));
  if (err) {
    goto FREE_RETURN;
  }

  Test = volBondCpn->second_startDate - volBondCpn->second_fixingDate;

  if (Test > 2) {
    err = Fill_SWAP_FOR_LGM2F(
        today, volBondCpn->second_fixingDate, volBondCpn->second_startDate,
        volBondCpn->first_SwapFreq, volBondCpn->first_SwapBasis,
        volBondCpn->first_SwapRefRate, &(swapForCMSAdj));
    if (err) {
      goto FREE_RETURN;
    }
  }

  err = MY_LGM2F_Compute_Mean_Variance_Covariance(
      today, volBondCpn->first_fixingDate, volBondCpn->payment_Date,
      &(lgm2fForMeanVarCov), l, yc);
  if (err) {
    goto FREE_RETURN;
  }

  if ((lgm2fForMeanVarCov.var1 > 0) && (lgm2fForMeanVarCov.var2 > 0)) {
    GQRho = lgm2fForMeanVarCov.covar /
            sqrt(lgm2fForMeanVarCov.var1 * lgm2fForMeanVarCov.var2);
  } else {
    GQRho = 0.0;
  }

  err = SWAP_FOR_LGM2F_Preliminary_Computation(
      today, volBondCpn->first_fixingDate, &(volBondCpn->firstSwap), l, yc);
  if (err) {
    goto FREE_RETURN;
  }

  err = SWAP_FOR_LGM2F_Preliminary_Computation(
      today, volBondCpn->first_fixingDate, &(volBondCpn->secondSwap), l, yc);
  if (err) {
    goto FREE_RETURN;
  }

  if (Test > 2) {
    err = SWAP_FOR_LGM2F_Preliminary_Computation(
        today, volBondCpn->first_fixingDate, &swapForCMSAdj, l, yc);
    if (err) {
      goto FREE_RETURN;
    }
  }

  //--------------------------------------------------------------------------
  //-------------------------	Gaussian Quadrature
  //----------------------------
  //--------------------------------------------------------------------------

  CMSAdj = 0;
  *VolBondCouponPrice = 0;
  for (i = 0; i < NGaussQuadrature; ++i) {
    for (j = 0; j < NGaussQuadrature; ++j) {
      fact1 =
          lgm2fForMeanVarCov.esp1 + sqrt(lgm2fForMeanVarCov.var1) * x[i + 1];
      fact2 = lgm2fForMeanVarCov.esp2 +
              sqrt(lgm2fForMeanVarCov.var2) *
                  (GQRho * x[i + 1] + sqrt(1 - GQRho * GQRho) * x[j + 1]);

      if (Test > 2) {
        CMSAdj = CMSAdjustementWithDelay(
            fact1, fact2, today, volBondCpn->first_fixingDate, swapForCMSAdj,
            CompSabrBetaVol, CompSabrBeta, volBondCpn->secondSwap, sabrBetaVol,
            sabrBeta, swapRho, l, yc);
      } else {
        CMSAdj = CMSAdjustementWithoutDelay(
            fact1, fact2, today, volBondCpn->first_fixingDate,
            volBondCpn->secondSwap, sabrBetaVol, sabrBeta, l, yc);
      }

      err = PayOff_VolBondCoupon(
          today, fact1, fact2, volBondCpn,
          exp(CMSAdj *
              (volBondCpn->second_fixingDate - volBondCpn->first_fixingDate) /
              365.25),
          sabrAlpha, sabrBeta, sabrBetaVol, sabrRho, &Payoff, l, yc);
      if (err) {
        goto FREE_RETURN;
      }

      *VolBondCouponPrice += w[i + 1] * w[j + 1] * Payoff;
    }
  }

  *VolBondCouponPrice = (*VolBondCouponPrice) * volBondCpn->notional *
                        swp_f_df(today, volBondCpn->payment_Date, yc);

  //--------------------------------------------------------------------------
  //-----------------------	End Of Gaussian Quadrature
  //--------------------
  //--------------------------------------------------------------------------

FREE_RETURN:

  if (&(lgm2fForMeanVarCov)) {
    Free_MY_LGM2F(&(lgm2fForMeanVarCov));
  }

  if (Test > 2) {
    if (&(swapForCMSAdj)) {
      Free_SWAP_FOR_LGM2F(&(swapForCMSAdj));
    }
  }

  return err;
}

Err Price_VolBond(long today,
                  //	The underlying
                  char *yc,      //	yc
                  TermStruct *l, // underlying
                  //	SABR Parameters
                  double *sabrBetaVol, double *sabrAlpha, double *sabrBeta,
                  double *sabrRho,
                  // Complementary Swap Sabr parameters
                  double *swapRho, double *CompSabrBeta,
                  double *CompSabrBetaVol,
                  // Numerical Parameter
                  int NGaussQuadrature, volBondStruct *volBond,
                  double *VolBondCpnPrices, double *VolBondPrice) {
  long i;
  double price;
  double temp;
  Err err = NULL;
  double *x = NULL, *w = NULL;
  x = (double *)calloc(NGaussQuadrature + 1, sizeof(double));
  if (!x) {
    err = "Error in Price_VolBond : Gauss Quadra x allocation failed";
    goto FREE_RETURN;
  }
  w = (double *)calloc(NGaussQuadrature + 1, sizeof(double));
  if (!w) {
    err = "Error in Price_VolBond : Gauss Quadra w allocation failed";
    goto FREE_RETURN;
  }

  err = HermiteStandard(x, w, NGaussQuadrature);
  if (err) {
    err = "Error in Price_VolBond : Gauss Quadra failed";
    goto FREE_RETURN;
  }

  price = 0;
  for (i = 0; i < volBond->NbOfCoupons; ++i) {
    err = Price_VolBondCoupon(
        today, yc, l, sabrBetaVol[i], sabrAlpha[i], sabrBeta[i], sabrRho[i],
        swapRho[i], CompSabrBeta[i], CompSabrBetaVol[i], NGaussQuadrature, x, w,
        &(volBond->VolBondCpn[i]), &temp);
    if (err) {
      goto FREE_RETURN;
    }
    price += temp;
    VolBondCpnPrices[i] = temp;
  }

  *VolBondPrice = price;

FREE_RETURN:

  if (x) {
    free(x);
  }
  if (w) {
    free(w);
  }

  return err;
}

Err Price_VolBondAutoCal(long today,
                         //	The underlying
                         char *yc, //	yield curve
                         char *vc, //	vol curve
                         double defaultLambda, char *refrate, double alpha,
                         double gamma, double rho, // LGM2F parameters
                         //	SABR Parameters
                         double *sabrBetaVol, double *sabrAlpha,
                         double *sabrBeta, double *sabrRho,
                         // Complementary Swap Sabr parameters
                         double *swapRho, double *CompSabrBeta,
                         double *CompSabrBetaVol,
                         // Numerical Parameter
                         int NGaussQuadrature, double *lambdaVec,
                         double *sigmaVec, volBondStruct *volBond,
                         double *cpnPrices, double *VolBondPrice) {
  long i;
  double price;
  double temp;
  Err err = NULL;
  TermStruct *ts = NULL;
  SrtReceiverType pay_rec1 = 1;
  SrtReceiverType pay_rec2 = 1;
  SrtMdlDim mdl_dim = TWO_FAC;
  double lambda;
  double sig;

  double **sig_data;
  int sig_cols = 4;
  int num_sigs = 1;
  double **tau_data;
  int tau_cols = 3;
  int num_taus = 1;
  SrtMdlType mdl_type = LGM;
  double beta = 0;
  double omega = 0;
  double *x = NULL, *w = NULL;

  x = (double *)calloc(NGaussQuadrature + 1, sizeof(double));
  if (!x) {
    err = "Error in Price_VolBondCoupon : Gauss Quadra x allocation failed";
    goto FREE_RETURN;
  }
  w = (double *)calloc(NGaussQuadrature + 1, sizeof(double));
  if (!w) {
    err = "Error in Price_VolBondCoupon : Gauss Quadra w allocation failed";
    goto FREE_RETURN;
  }

  err = HermiteStandard(x, w, NGaussQuadrature);
  if (err) {
    err = "Error in Price_VolBondCoupon : Gauss Quadra failed";
    goto FREE_RETURN;
  }

  sig_data = dmatrix(0, 3, 0, 0);
  tau_data = dmatrix(0, 2, 0, 0);

  price = 0;
  sig = 0.0085 / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
  lambda = defaultLambda;
  for (i = 0; i < volBond->NbOfCoupons; ++i) {
    err = VolBondCalibrationOld(yc, vc,
                                (volBond->VolBondCpn[i]).first_SwapRefRate,
                                (volBond->VolBondCpn[i]).first_fixingDate,
                                (volBond->VolBondCpn[i]).first_endDate,
                                (volBond->VolBondCpn[i]).first_startDate,
                                (volBond->VolBondCpn[i]).first_SwapFreq,
                                (volBond->VolBondCpn[i]).first_SwapBasis, alpha,
                                gamma, rho, &lambda, &sig);
    if (err) {
      goto FREE_RETURN;
    }

    sig_data[0][0] = max(today + 1, (volBond->VolBondCpn[i]).first_fixingDate);
    sig_data[1][0] = sig;
    sig_data[2][0] = sig * alpha;
    sig_data[3][0] = rho;

    tau_data[0][0] = max(today + 1, (volBond->VolBondCpn[i]).first_fixingDate);
    tau_data[1][0] = 1 / lambda;
    tau_data[2][0] = 1 / (lambda + gamma);

    err = srt_f_init_IRM_TwoFac_TermStruct(
        &ts, today, sig_data, sig_cols, num_sigs, tau_data, tau_cols, num_taus,
        mdl_type, beta, alpha, gamma, rho, omega);
    if (err) {
      goto FREE_RETURN;
    }

    lambdaVec[i] = lambda;
    sigmaVec[i] = sig;

    err = Price_VolBondCoupon(
        today, yc, ts, sabrBetaVol[i], sabrAlpha[i], sabrBeta[i], sabrRho[i],
        swapRho[i], CompSabrBeta[i], CompSabrBetaVol[i], NGaussQuadrature, x, w,
        &(volBond->VolBondCpn[i]), &temp);
    if (err) {
      goto FREE_RETURN;
    }
    price += temp;
    cpnPrices[i] = temp;
  }

  *VolBondPrice = price;

FREE_RETURN:

  if (sig_data) {
    free_dmatrix(sig_data, 0, 3, 0, 0);
  }
  if (tau_data) {
    free_dmatrix(tau_data, 0, 2, 0, 0);
  }

  return err;
}

Err VolBond_Price(
    long today, char *yc, TermStruct *l,
    // VolBond Parameters//
    long NbCoupons, double *notional, char **first_SwapFreq,
    char **first_SwapBasis, char **first_SwapRefRate, char **second_SwapFreq,
    char **second_SwapBasis, char **second_SwapRefRate, long *first_fixingDate,
    long *first_startDate, long *first_endDate, long *second_fixingDate,
    long *second_startDate, long *second_endDate, long *payment_Date,
    double *alphaC, double *alphaP, double *betaC, double *betaP,
    double *strikeC, double *strikeP,
    // Main Swap Sabr Parameters//
    double *sabrBetaVol, double *sabrAlpha, double *sabrBeta, double *sabrRho,
    // Complementary Swap Sabr Parameters//
    double *swapRho, double *CompSabrBeta, double *CompSabrBetaVol,
    // Numerical Parameter//
    int NGaussQuadrature, double *cpnPrices, double *priceValue) {
  Err err = NULL;
  double price;
  volBondStruct volBond;

  err = Fill_VolBondStruct(
      today, NbCoupons, notional, first_SwapFreq, first_SwapBasis,
      first_SwapRefRate, second_SwapFreq, second_SwapBasis, second_SwapRefRate,
      first_fixingDate, first_startDate, first_endDate, second_fixingDate,
      second_startDate, second_endDate, payment_Date, alphaC, alphaP, betaC,
      betaP, strikeC, strikeP, &volBond);
  if (err) {
    goto FREE_RETURN;
  }

  err = Price_VolBond(today, yc, l, sabrBetaVol, sabrAlpha, sabrBeta, sabrRho,
                      swapRho, CompSabrBeta, CompSabrBetaVol, NGaussQuadrature,
                      &volBond, cpnPrices, &price);
  if (err) {
    goto FREE_RETURN;
  }

  *priceValue = price;

FREE_RETURN:

  Free_VolBondStruct(&volBond);

  return err;
}

Err VolBond_PriceAutoCal(
    long today, char *yc, char *vc, char *refrate, double defaultLambda,
    double lgm2f_alpha, double lgm2f_gamma, double lgm2f_rho,
    // VolBond Parameters//
    long NbCoupons, double *notional, char **first_SwapFreq,
    char **first_SwapBasis, char **first_SwapRefRate, char **second_SwapFreq,
    char **second_SwapBasis, char **second_SwapRefRate, long *first_fixingDate,
    long *first_startDate, long *first_endDate, long *second_fixingDate,
    long *second_startDate, long *second_endDate, long *payment_Date,
    double *alphaC, double *alphaP, double *betaC, double *betaP,
    double *strikeC, double *strikeP,
    // Main Swap Sabr Parameters//
    double *sabrBetaVol, double *sabrAlpha, double *sabrBeta, double *sabrRho,
    // Complementary Swap Sabr Parameters//
    double *swapRho, double *CompSabrBeta, double *CompSabrBetaVol,
    // Numerical Parameter//
    int NGaussQuadrature, double *lambdaVec, double *sigmaVec,
    double *cpnPrices, double *priceValue) {
  Err err = NULL;
  double price;
  volBondStruct volBond;

  err = Fill_VolBondStruct(
      today, NbCoupons, notional, first_SwapFreq, first_SwapBasis,
      first_SwapRefRate, second_SwapFreq, second_SwapBasis, second_SwapRefRate,
      first_fixingDate, first_startDate, first_endDate, second_fixingDate,
      second_startDate, second_endDate, payment_Date, alphaC, alphaP, betaC,
      betaP, strikeC, strikeP, &volBond);
  if (err) {
    goto FREE_RETURN;
  }

  err = Price_VolBondAutoCal(today, yc, vc, defaultLambda, refrate, lgm2f_alpha,
                             lgm2f_gamma, lgm2f_rho, sabrBetaVol, sabrAlpha,
                             sabrBeta, sabrRho, swapRho, CompSabrBeta,
                             CompSabrBetaVol, NGaussQuadrature, lambdaVec,
                             sigmaVec, &volBond, cpnPrices, &price);
  if (err) {
    goto FREE_RETURN;
  }

  *priceValue = price;

FREE_RETURN:

  Free_VolBondStruct(&volBond);

  return err;
}
