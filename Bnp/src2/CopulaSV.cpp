
/*******************************************************************************
 *
 * FUNCTION      : Copula SV
 *
 * PURPOSE       : Copula function in a double Heston Model
 *
 * DESCRIPTION   : uses a MonteCarlo type numerical integration with SOBOL
 *				  technique to reduce the variance of the result
 *
 * CALLS         : Sobol_init  , Sobol_Cube
 *
 *******************************************************************************/

#include "CopulaSV.h"
#include "CopulaGaussian.h"
#include "OPFNCTNS.h"
#include "SABRCalibRRBT.h"
#include "UTALLHDR.h"
#include "math.h"
#include "num_h_proba.h"
#include "num_h_sobol.h"
#include "opsabrgenericinterp.h"

void copula_sv_set_default_num_params(COPULASV_PARAMS sParams) {
  sParams->lNumPaths = 10000;
  sParams->dMinTime = 1.0 / 24.0;
  sParams->iMCType = 0;

  sParams->lNbPoints = (long)(pow(2.0, 7.0) + 0.5);
  sParams->iStudDegree = 4;
  sParams->iNConv = 0;
  sParams->eMCType = SOBOL;

  sParams->iAdaptGrid = 1;
  sParams->lAdaptNbPoints = 50;
  sParams->iUseOldGridMethod = 1;
  sParams->dNbStd = 5.0;
  sParams->dMaxError = 0.0005;

  sParams->iCalcImpliedRho = 0;

  sParams->iUseSABR = 0;
  sParams->dPi = 0.25;
  sParams->iForceLogBeta = 0;
  sParams->dCalibBetaStd = 1.0;

  sParams->iHasSavedGaussian = 0;
  sParams->iFreeSavedGaussian = 1;
  sParams->iNbFwd = 0;
  sParams->dSavedGaussian = NULL;

  sParams->iIntegSkipFirst = 0;

  /* For debug */
  sParams->cFileName[0] = '\0';
}

/*	The Approximation:
        X(T) = dConst + dCoef1_1 * U + dCoef2_1 * V + dCoef2_2 * V^2 + dCross *
   U * V U and V follow independant normal distributions */

Err copula_sv_get_polynomial_coefs(/* Model Parameters */
                                   double dTime, double dSigma, double dAlpha,
                                   double dRho, double dVar0, double dVarInf,
                                   double dLambda,

                                   /*	Outputs */
                                   double *dConst, double *dCoef1_1,
                                   double *dCoef2_1, double *dCross) {
  Err err = NULL;

  double dSigma0, dRho2, dRhoAlpha, dMeanLevel;
  double dSqTime;

  /* Constants */
  dSigma0 = sqrt(dVar0);
  dRho2 = sqrt(1.0 - dRho * dRho);
  dRhoAlpha = dRho * dAlpha;
  dMeanLevel = dLambda * dVarInf;
  dSqTime = sqrt(dTime);

  /* Outputs */
  *dConst =
      -0.5 * dVar0 * dSigma * dSigma * dTime +
      0.25 * dSigma * dSigma * (dMeanLevel - dLambda * dVar0) * dTime * dTime -
      0.25 * dRho * dAlpha * dSigma * dTime;

  *dCoef1_1 = dSigma *
              (dSigma0 + 1.0 / dSigma0 *
                             ((dMeanLevel - dLambda * dVar0) / 4.0 -
                              dAlpha * dAlpha / 32.0) *
                             dTime) *
              dSqTime;

  *dCoef2_1 = -dSigma * dAlpha *
              (dSigma * dSigma0 / 4.0 + dRho * dAlpha / dSigma0 / 32.0) *
              dTime * dSqTime;

  *dCross = 0.25 * dSigma * dAlpha * dTime;

  return err;
}

Err copula_sv_single_pricer(/* Product Parameters */
                            double dTime, int iNbStrike, double *dStrikes,

                            /* Model Parameters */
                            double dForward, double dShift, double dSigma,
                            double dAlpha, double dRho, double dVar0,
                            double dVarInf, double dLambda,

                            /* Numerical Params */
                            COPULASV_PARAMS sParams,

                            /* Outputs */
                            double *dValue) {
  Err err = NULL;

  long iNumPaths, rand;
  double step, prob;
  long seed = -123456789;
  long i, k;

  double dConst, dCoef1_1, dCoef2_1, dCross, dRho2;

  double X;

  double *dInitGauss = NULL, **dBrownian = NULL;

  /* Balantisam Simulation */
  iNumPaths = 2 * ((long)(sParams->lNumPaths / 2)) + 1;
  dInitGauss = calloc(iNumPaths, sizeof(double));
  dBrownian = dmatrix(0, 1, 0, iNumPaths - 1);

  if (!dInitGauss || !dBrownian) {
    err = "Memory allocation faillure in copula_sv_single_pricer_balantisam";
    goto FREE_RETURN;
  }

  /* Get the model new parameters */
  err = copula_sv_get_polynomial_coefs(dTime, dSigma, dAlpha, dRho, dVar0,
                                       dVarInf, dLambda, &dConst, &dCoef1_1,
                                       &dCoef2_1, &dCross);

  if (err)
    goto FREE_RETURN;

  dConst += log(dForward + dShift);
  dRho2 = sqrt(1.0 - dRho * dRho);

  /* Random numbers generation */
  if (sParams->iMCType == 0) {
    /* Gauss initialisation */
    iNumPaths -= 1;
    iNumPaths /= 2;
    step = 0.5 / (iNumPaths + 1);
    prob = step;

    /* Generation of the fractiles of the gaussian */
    for (i = 0; i < iNumPaths; i++) {
      dInitGauss[i] = inv_cumnorm_fast(prob);
      dInitGauss[iNumPaths + i + 1] = -dInitGauss[i];
      prob += step;
    }

    dInitGauss[iNumPaths] = 0.0;
    iNumPaths *= 2;
    iNumPaths += 1;

    /* Generate the random numbers */
    for (k = 0; k < 2; k++) {
      for (i = 0; i < iNumPaths - 1; i++) {
        /* rand = random_int(nbPaths-1-i  , &seed) + i; */
        rand = i + (int)((iNumPaths - i) * uniform(&seed));
        dBrownian[k][i] = dInitGauss[rand];
        dInitGauss[rand] = dInitGauss[i];
        dInitGauss[i] = dBrownian[k][i];
      }

      dBrownian[k][iNumPaths - 1] = dInitGauss[iNumPaths - 1];
    }
  }

  /* MC simulation */
  memset(dValue, 0, iNbStrike * sizeof(double));

  for (i = 0; i < iNumPaths; i++) {
    /* Correlation */
    dBrownian[1][i] = dRho * dBrownian[0][i] + dRho2 * dBrownian[1][i];

    /* Reconstruction */
    X = dConst + dCoef1_1 * dBrownian[0][i] +
        dBrownian[1][i] * (dCoef2_1 + dCross * dBrownian[0][i]);

    X = exp(X) - dShift;

    for (k = 0; k < iNbStrike; k++) {
      dValue[k] += max(X - dStrikes[k], 0.0);
    }
  }

  for (k = 0; k < iNbStrike; k++) {
    dValue[k] /= iNumPaths;
  }

FREE_RETURN:

  if (dInitGauss)
    free(dInitGauss);
  if (dBrownian)
    free_dmatrix(dBrownian, 0, 1, 0, iNumPaths - 1);

  return err;
}

/* Transformation of SABR into approx-Heston */
Err copula_sv_from_sabr_to_heston(double dMaturity, double dFwd, double dSigma,
                                  double dAlpha, double dBeta, double dRho,
                                  SrtDiffusionType eTypeInput,

                                  double *dSigmaH, double *dAlphaH,
                                  double *dShiftH, double *dRhoH,
                                  double *dVar0H, double *dVarInfH,
                                  double *dLambdaH) {
  Err err = NULL;

  /* First get the Sigma Beta */
  srt_f_optsarbvol(dFwd, dFwd, dMaturity, dSigma, dAlpha, dBeta, dRho,
                   eTypeInput, SRT_BETAVOL, &dSigma);

  /* Then Get the shift */
  *dShiftH = dFwd * (1.0 - dBeta) / dBeta;

  /* The Heston Vol */
  *dSigmaH = dSigma * exp(dBeta * log(dFwd)) / (dFwd + (*dShiftH));

  *dAlphaH = 2.0 * dAlpha;
  *dRhoH = dRho;
  *dVar0H = 1.0;
  *dVarInfH = 1.0;
  *dLambdaH = 0.001;

  return err;
}

double copula_sv_get_SABR_get_cumul(double dMaturity, double dForward,
                                    double dStrike, double dSigmaBeta,
                                    double dAlpha, double dBeta, double dRho,
                                    double dShift) {
  double price1, price2, res;

  price1 = srt_f_optblkschbetastochquick(
      dForward, dStrike - dShift, dMaturity, dSigmaBeta, dAlpha, dBeta, dRho,
      1.0, SRT_BETAVOL, SRT_LOGNORMAL, SRT_PUT, PREMIUM);

  price2 = srt_f_optblkschbetastochquick(
      dForward, dStrike + dShift, dMaturity, dSigmaBeta, dAlpha, dBeta, dRho,
      1.0, SRT_BETAVOL, SRT_LOGNORMAL, SRT_PUT, PREMIUM);

  res = (price2 - price1) / (2.0 * dShift);

  return res;
}

Err copula_sv_get_SABR_cumulative(double dMaturity, double dForward,
                                  double dSigmaBeta, double dAlpha,
                                  double dBeta, double dRho,

                                  int iAdaptGrid, long lAdaptPoints,
                                  int iUsePoints,

                                  long lNumPoints, double *dPoints,
                                  double *dCumulative) {
  Err err = NULL;

  double *dFirstGuessPoints = NULL, *dFirstGuessCumul = NULL;

  long lFirstGuessNbPoints;
  double dNbStdDev, dLogStd, dLogCvx;
  double dGridCoef, dShift;
  double dNextCumul, dCurrentCumul;
  double dProba, dProbaStep, dTemp;

  int i;

  if (iAdaptGrid && !iUsePoints) {
    lFirstGuessNbPoints = lAdaptPoints;
    dFirstGuessPoints = calloc(lFirstGuessNbPoints, sizeof(double));
    dFirstGuessCumul = calloc(lFirstGuessNbPoints, sizeof(double));

    if (!dFirstGuessCumul || !dFirstGuessCumul) {
      err = "Memory allocation faillure in copula_sv_get_SABR_cumulative";
      goto FREE_RETURN;
    }
  } else {
    lFirstGuessNbPoints = lNumPoints;
    dFirstGuessPoints = dPoints;
    dFirstGuessCumul = dCumulative;
  }

  /* First Guess */

  srt_f_optsarbvol(dForward, dForward, dMaturity, dSigmaBeta, dAlpha, dBeta,
                   dRho, SRT_BETAVOL, SRT_LOGNORMAL, &dLogStd);

  dLogStd *= sqrt(dMaturity);
  dLogCvx = -0.5 * dLogStd * dLogStd;

  if (!iUsePoints) {
    dFirstGuessPoints[0] = dForward * exp(dLogCvx - 5.0 * dLogStd);
    dGridCoef = exp(2.0 * 5.0 * dLogStd / lFirstGuessNbPoints);
  }

  for (i = 1; i < lFirstGuessNbPoints; i++) {
    if (!iUsePoints) {
      dFirstGuessPoints[i] = dFirstGuessPoints[i - 1] * dGridCoef;
    }

    dShift = (dFirstGuessPoints[i] - dFirstGuessPoints[i - 1]) / 2.0;

    dFirstGuessCumul[i] =
        copula_sv_get_SABR_get_cumul(dMaturity, dForward, dFirstGuessPoints[i],
                                     dSigmaBeta, dAlpha, dBeta, dRho, dShift);
  }

  dNbStdDev = max(min(5.0 * sqrt(dMaturity), 10.0), 5.0);
  dFirstGuessPoints[0] = dForward * exp(dLogCvx - dNbStdDev * dLogStd);
  dFirstGuessCumul[0] = 1.0E-10;

  /* Correct SABR problems... */
  i = (long)(lFirstGuessNbPoints / 2.0);

  while (i > 0 && dFirstGuessCumul[i] > dFirstGuessCumul[i - 1]) {
    i--;
  }

  i--;

  if (i >= 0) {
    dShift = (dFirstGuessPoints[i + 2] - dFirstGuessPoints[i + 1]) / 2.0;

    dNextCumul = copula_sv_get_SABR_get_cumul(
        dMaturity, dForward, dFirstGuessPoints[i + 1], dSigmaBeta, 0.0, dBeta,
        0.0, dShift);
    for (i = i; i >= 0; i--) {
      /* correction through beta distribution */
      dShift = (dFirstGuessPoints[i + 1] - dFirstGuessPoints[i]) / 2.0;

      dCurrentCumul = copula_sv_get_SABR_get_cumul(
          dMaturity, dForward, dFirstGuessPoints[i], dSigmaBeta, 0.0, dBeta,
          0.0, dShift);

      dFirstGuessCumul[i] =
          dFirstGuessCumul[i + 1] * dCurrentCumul / dNextCumul;

      dNextCumul = dCurrentCumul;
    }
  }

  if (iAdaptGrid && !iUsePoints) {
    /* Derive now an adapted grid */
    dProba = dFirstGuessCumul[0] * 1.000001;
    dProbaStep =
        (dFirstGuessCumul[lFirstGuessNbPoints - 1] - dFirstGuessCumul[0]) /
        (lNumPoints - 1.0);

    for (i = 0; i < lNumPoints; i++) {
      dPoints[i] = interp(dFirstGuessCumul, dFirstGuessPoints,
                          lFirstGuessNbPoints, dProba, 0, &dTemp);

      dProba += dProbaStep;
    }

    /* recall */
    err = copula_sv_get_SABR_cumulative(dMaturity, dForward, dSigmaBeta, dAlpha,
                                        dBeta, dRho, 0, 0, 1, lNumPoints,
                                        dPoints, dCumulative);
  }

FREE_RETURN:

  if (iAdaptGrid) {
    if (dFirstGuessPoints)
      free(dFirstGuessPoints);
    if (dFirstGuessCumul)
      free(dFirstGuessCumul);
  }

  return err;
}

Err copula_sv_basket_SABR(int iNbFwd, double *dFwds, double *dWeights,
                          int iNbStrikes, double *dStrike, double *dVols,
                          double *dAlpha, double *dBeta, double *dRho,
                          double **dCorrelation, double dMaturity,
                          SrtCallPutType eCallPut, SrtDiffusionType eTypeInput,
                          COPULASV_PARAMS sParams, double *dPremium) {
  Err err = NULL;
  double **dResult = NULL, **dX = NULL, **dCumulative = NULL;

  double *dShiftSV = NULL, *dSigmaSV = NULL, *dAlphaSV = NULL, *dRhoSV = NULL,
         *dVar0SV = NULL, *dLambdaSV = NULL;

  long *lNbPoints = NULL;

  double dSigmaBeta, dSigmaLog, dNewAlpha, dNewBeta, dNewRho;
  double dPi, dForward1, dSigmaBeta1, dForward2, dSigmaBeta2, dCalibRes;
  double dBasket, dFactor, dCashFlow;
  long i, j;
  FILE *stream;

  dResult = dmatrix(0, sParams->lNumPaths - 1, 0, iNbFwd - 1);
  dX = calloc(iNbFwd, sizeof(double *));
  dCumulative = calloc(iNbFwd, sizeof(double *));
  lNbPoints = calloc(iNbFwd, sizeof(long *));

  dShiftSV = calloc(iNbFwd, sizeof(double));
  dSigmaSV = calloc(iNbFwd, sizeof(double));
  dAlphaSV = calloc(iNbFwd, sizeof(double));
  dRhoSV = calloc(iNbFwd, sizeof(double));
  dVar0SV = calloc(iNbFwd, sizeof(double));
  dLambdaSV = calloc(iNbFwd, sizeof(double));

  if (!dResult || !dX || !dCumulative || !dShiftSV || !dSigmaSV || !dAlphaSV ||
      !dRhoSV || !dVar0SV || !dLambdaSV) {
    err = "Memory allocation faillure in copula_gaussian_basket_SABR";
    goto FREE_RETURN;
  }

  for (i = 0; i < iNbFwd; i++) {
    dX[i] = calloc(sParams->lNbPoints, sizeof(double));
    dCumulative[i] = calloc(sParams->lNbPoints, sizeof(double));
    lNbPoints[i] = sParams->lNbPoints;

    if (!dX[i] || !dCumulative[i]) {
      err = "Memory allocation faillure in copula_gaussian_basket_SABR";
      goto FREE_RETURN;
    }
  }

  /* Constant Initialisation */
  if (eCallPut == SRT_CALL) {
    dFactor = 1.0;
  } else {
    dFactor = -1.0;
  }

  for (j = 0; j < iNbFwd; j++) {
    /* Get Sigma Beta */
    srt_f_optsarbvol(dFwds[j], dFwds[j], dMaturity, dVols[j], dAlpha[j],
                     dBeta[j], dRho[j], eTypeInput, SRT_BETAVOL, &dSigmaBeta);

    srt_f_optsarbvol(dFwds[j], dFwds[j], dMaturity, dSigmaBeta, dAlpha[j],
                     dBeta[j], dRho[j], SRT_BETAVOL, SRT_LOGNORMAL, &dSigmaLog);

    if (sParams->iForceLogBeta) {
      dNewBeta = 1.0;

      /* Recalibrate the Smile with Beta = 1 */
      err = transform_sabr_beta(
          dMaturity, dFwds[j], dSigmaBeta, dAlpha[j], dBeta[j], dRho[j],
          SRT_BETAVOL, &dSigmaBeta, &dNewAlpha, dNewBeta, &dNewRho,
          sParams->dCalibBetaStd, 1, 0.0001, 10, &dCalibRes);

      if (err)
        goto FREE_RETURN;
    } else {
      dNewAlpha = dAlpha[j];
      dNewBeta = dBeta[j];
      dNewRho = dRho[j];
    }

    if (sParams->iUseSABR) {
      err = copula_gaussian_get_SABR_cumulative(
          dMaturity, dFwds[j], dSigmaBeta, dNewAlpha, dNewBeta, dNewRho,
          sParams->iAdaptGrid, sParams->lAdaptNbPoints, 0, sParams->lNbPoints,
          dX[j], dCumulative[j]);

      if (err)
        goto FREE_RETURN;
    } else {
      dPi = sParams->dPi;

      err = BMMCalibOnSabrStates(dFwds[j], dMaturity, dSigmaLog, dNewBeta,
                                 dNewAlpha, dNewRho, sParams->dPi, 1.0,
                                 &dForward1, &dForward2, &dSigmaBeta1,
                                 &dSigmaBeta2, &dPi, SRT_LOGNORMAL, &dCalibRes);

      if (err)
        goto FREE_RETURN;

      err = copula_gaussian_get_BMM_linterp_cumulative(
          dMaturity, dFwds[j], dForward1, dSigmaBeta1, dForward2, dSigmaBeta2,
          dNewBeta, dPi, sParams->dMaxError, sParams->dNbStd, sParams->dNbStd,
          &(lNbPoints[j]), &(dX[j]), &(dCumulative[j]));

      if (err)
        goto FREE_RETURN;
    }
  }

  if (strlen(sParams->cFileName) > 0) {
    stream = fopen(sParams->cFileName, "w+");
    if (!stream)
      goto FREE_RETURN;

    for (j = 0; j < iNbFwd; j++) {
      for (i = 0; i < lNbPoints[j]; i++) {
        fprintf(stream, "%f	%f\n", dX[j][i], dCumulative[j][i]);
      }

      fprintf(stream, "\n");
    }

    fclose(stream);
  }

  /* Set the variables and the transformation */
  for (i = 0; i < iNbFwd; i++) {
    err = copula_sv_from_sabr_to_heston(
        dMaturity, dFwds[i], dVols[i], dAlpha[i], dBeta[i], dRho[i], eTypeInput,
        &(dSigmaSV[i]), &(dAlphaSV[i]), &(dShiftSV[i]), &(dRhoSV[i]),
        &(dVar0SV[i]), &(dVar0SV[i]), &(dLambdaSV[i]));

    if (err)
      goto FREE_RETURN;
  }

  err = copula_sv_numer(iNbFwd, dX, dCumulative, lNbPoints, dMaturity, dFwds,
                        dShiftSV, dSigmaSV, dAlphaSV, dRhoSV, dVar0SV, dVar0SV,
                        dLambdaSV, dCorrelation, sParams, dResult);

  if (err)
    goto FREE_RETURN;

  /* Compute the Price and the normal vols via a Monte Carlo Methods */
  memset(dPremium, 0, iNbStrikes * sizeof(double));

  if (sParams->iCalcImpliedRho && iNbStrikes >= 5) {
    for (i = sParams->iIntegSkipFirst; i < sParams->lNumPaths; i++) {
      dPremium[0] += dResult[i][0] / (i - sParams->iIntegSkipFirst + 1.0);
      dPremium[1] +=
          dResult[i][0] * dResult[i][0] / (i - sParams->iIntegSkipFirst + 1.0);
      dPremium[2] += dResult[i][1] / (i - sParams->iIntegSkipFirst + 1.0);
      dPremium[3] +=
          dResult[i][1] * dResult[i][1] / (i - sParams->iIntegSkipFirst + 1.0);
      dPremium[4] +=
          dResult[i][0] * dResult[i][1] / (i - sParams->iIntegSkipFirst + 1.0);
    }
  } else {
    for (i = sParams->iIntegSkipFirst; i < sParams->lNumPaths; i++) {
      dBasket = 0.0;

      for (j = 0; j < iNbFwd; j++) {
        dBasket += dWeights[j] * dResult[i][j];
      }

      for (j = 0; j < iNbStrikes; j++) {
        dCashFlow = DMAX(dFactor * (dBasket - dStrike[j]), 0.0);

        dPremium[j] +=
            (dCashFlow - dPremium[j]) / (i - sParams->iIntegSkipFirst + 1.0);
      }
    }
  }

FREE_RETURN:

  if (dResult)
    free_dmatrix(dResult, 0, sParams->lNumPaths - 1, 0, iNbFwd - 1);

  if (dX) {
    for (i = 0; i < iNbFwd; i++) {
      if (dX[i])
        free(dX[i]);
    }

    free(dX);
  }

  if (dCumulative) {
    for (i = 0; i < iNbFwd; i++) {
      if (dCumulative[i])
        free(dCumulative[i]);
    }

    free(dCumulative);
  }

  if (lNbPoints)
    free(lNbPoints);

  if (dShiftSV)
    free(dShiftSV);
  if (dSigmaSV)
    free(dSigmaSV);
  if (dAlphaSV)
    free(dAlphaSV);
  if (dRhoSV)
    free(dRhoSV);
  if (dVar0SV)
    free(dVar0SV);
  if (dLambdaSV)
    free(dLambdaSV);

  return err;
}

Err copula_sv_numer(/* Marginales Distributions */
                    int iNbFwd, double **xa, double **ya, long *lNbPoints,

                    /* Copula Parameters */
                    double dTime, double *dForward, double *dShift,
                    double *dSigma, double *dAlpha, double *dRho, double *dVar0,
                    double *dVarInf, double *dLambda, double **dCorrMatrix,

                    /* Numerical Parameters */
                    COPULASV_PARAMS sParams,

                    /* Result */
                    double **dResult)

{
  long i, idum = -8935807;
  int j, k, eff_d;
  double **GaussSample = NULL, **XtSample = NULL, *dBrownian = NULL,
         **dSortMatrix = NULL, **dSqCorrMatrix = NULL;

  double *dConst = NULL, *dCoef1_1 = NULL, *dCoef2_1 = NULL, *dCross = NULL;

  double IntermRes, ProbaStep, q;

  long nNumPoints;
  double dImpliedRho;

  Err err = NULL;

  /* Memory allocation */
  eff_d = 2 * iNbFwd;
  nNumPoints = 2 * sParams->lNbPoints + 1;
  GaussSample = dmatrix(0, sParams->lNumPaths, 0, eff_d - 1);
  XtSample = dmatrix(0, iNbFwd - 1, 0, sParams->lNumPaths);
  dSqCorrMatrix = dmatrix(0, eff_d - 1, 0, eff_d - 1);
  dBrownian = calloc(eff_d, sizeof(double));
  dSortMatrix = dmatrix(0, 1, 0, sParams->lNumPaths);

  dConst = calloc(iNbFwd, sizeof(double));
  dCoef1_1 = calloc(iNbFwd, sizeof(double));
  dCoef2_1 = calloc(iNbFwd, sizeof(double));
  dCross = calloc(iNbFwd, sizeof(double));

  if (!GaussSample || !XtSample || !dSqCorrMatrix || !dBrownian ||
      !dSortMatrix || !dConst || !dCoef1_1 || !dCoef2_1 || !dCross) {
    err = "Memory allocation faillure in copula_sv_numer";
    goto FREE_RETURN;
  }

  /* Get the coeffs */
  for (j = 0; j < iNbFwd; j++) {
    err = copula_sv_get_polynomial_coefs(
        dTime, dSigma[j], dAlpha[j], dRho[j], dVar0[j], dVarInf[j], dLambda[j],
        &(dConst[j]), &(dCoef1_1[j]), &(dCoef2_1[j]), &(dCross[j]));

    if (err)
      goto FREE_RETURN;
  }

  ProbaStep = 1.0 / sParams->lNumPaths;

  /* performs Chol decomposition of the correlation matrix corr_mtx */
  err = choldc(eff_d, dCorrMatrix, dSqCorrMatrix);

  if (err)
    goto FREE_RETURN;

  switch (sParams->eMCType) {
  case ABS: // uses ABS cube
    break;

  case SOBOLBM: // uses SOBOL with Box-Muller algorithm
    break;

  case RANDOM_GAUSS:
    break;

  default: // by default uses SOBOL with Normal Inversion

    GetSobolMatrix(sParams->lNumPaths, eff_d, GaussSample);

    //////////  Step3: Generates the matrix of r.v. compatible with the Copula
    ///and the marginals
    for (i = 0; i < sParams->lNumPaths; i++) {
      /* Get the uncorrelated gaussians */
      for (j = 0; j < eff_d; j++) {
        GaussSample[i][j] = inv_cumnorm_fast(GaussSample[i][j]);
      }

      /* Correlate the gaussians */
      for (j = 0; j < eff_d; j++) {
        dBrownian[j] = 0.0;

        for (k = 0; k < eff_d; k++) {
          dBrownian[j] += dSqCorrMatrix[j][k] * GaussSample[i][k];
        }
      }

      /* Get the true variables */
      for (j = 0; j < iNbFwd; j++) {
        XtSample[j][i] =
            dConst[j] + dCoef1_1[j] * dBrownian[j] +
            dBrownian[j + iNbFwd] * (dCoef2_1[j] + dCross[j] * dBrownian[j]);
      }
    }

    for (j = 0; j < iNbFwd; j++) {
      /* Calculates the empirical distribution */
      for (i = 0; i < sParams->lNumPaths; i++) {
        dSortMatrix[0][i] = i;
      }

      memcpy(dSortMatrix[1], XtSample[j], sParams->lNumPaths * sizeof(double));

      num_f_sort_matrix(sParams->lNumPaths, 2, 1, dSortMatrix);

      IntermRes = ProbaStep;

      for (i = 0; i < sParams->lNumPaths; i++) {
        dResult[(long)(dSortMatrix[0][i])][j] =
            interp_columns(ya, xa, lNbPoints[j], IntermRes, 0, &q, j);

        IntermRes += ProbaStep;
      }
    }

    break;
  }

  /* Extra informations */
  if (sParams->iCalcImpliedRho) {
    err = copula_sv_get_implied_linear_correl(
        iNbFwd, 0, 1, dCoef1_1, dCoef2_1, dCross, dCorrMatrix, &dImpliedRho);

    if (err)
      goto FREE_RETURN;
  }

FREE_RETURN:

  if (GaussSample)
    free_dmatrix(GaussSample, 0, sParams->lNumPaths, 0, eff_d - 1);
  if (XtSample)
    free_dmatrix(XtSample, 0, iNbFwd - 1, 0, sParams->lNumPaths);
  if (dSqCorrMatrix)
    free_dmatrix(dSqCorrMatrix, 0, eff_d - 1, 0, eff_d - 1);
  if (dBrownian)
    free(dBrownian);
  if (dSortMatrix)
    free_dmatrix(dSortMatrix, 0, 1, 0, sParams->lNumPaths);

  if (dConst)
    free(dConst);
  if (dCoef1_1)
    free(dCoef1_1);
  if (dCoef2_1)
    free(dCoef2_1);
  if (dCross)
    free(dCross);

  return err;
}

Err copula_sv_quadrat(/* Forward */
                      int iNbFwd, double *dWeights,

                      /* Fractile */
                      int iNbFract, double *dFract,

                      /* Copula Parameters */
                      double dTime, double *dForward, double *dShift,
                      double *dSigma, double *dAlpha, double *dRho,
                      double *dVar0, double *dVarInf, double *dLambda,
                      double **dCorrMatrix,

                      /* Numerical Parameters */
                      COPULASV_PARAMS sParams,

                      /* Result */
                      double *dResult)

{
  long i, idum = -8935807;
  int j, k, eff_d;
  double **GaussSample = NULL, **XtSample = NULL, *dBrownian = NULL,
         **dSortMatrix = NULL, **dSqCorrMatrix = NULL;

  double *dConst = NULL, *dCoef1_1 = NULL, *dCoef2_1 = NULL, *dCross = NULL;

  double IntermRes, ProbaStep, dCashFlow;

  Err err = NULL;

  /* Memory allocation */
  eff_d = 2 * iNbFwd;
  GaussSample = dmatrix(0, sParams->lNumPaths, 0, eff_d - 1);
  XtSample = dmatrix(0, iNbFwd - 1, 0, sParams->lNumPaths);
  dSqCorrMatrix = dmatrix(0, eff_d - 1, 0, eff_d - 1);
  dBrownian = calloc(eff_d, sizeof(double));
  dSortMatrix = dmatrix(0, 1, 0, sParams->lNumPaths);

  dConst = calloc(iNbFwd, sizeof(double));
  dCoef1_1 = calloc(iNbFwd, sizeof(double));
  dCoef2_1 = calloc(iNbFwd, sizeof(double));
  dCross = calloc(iNbFwd, sizeof(double));

  if (!GaussSample || !XtSample || !dSqCorrMatrix || !dBrownian ||
      !dSortMatrix || !dConst || !dCoef1_1 || !dCoef2_1 || !dCross) {
    err = "Memory allocation faillure in copula_sv_numer";
    goto FREE_RETURN;
  }

  /* Get the coeffs */
  for (j = 0; j < iNbFwd; j++) {
    err = copula_sv_get_polynomial_coefs(
        dTime, dSigma[j], dAlpha[j], dRho[j], dVar0[j], dVarInf[j], dLambda[j],
        &(dConst[j]), &(dCoef1_1[j]), &(dCoef2_1[j]), &(dCross[j]));

    if (err)
      goto FREE_RETURN;
  }

  ProbaStep = 1.0 / sParams->lNumPaths;

  /* performs Chol decomposition of the correlation matrix corr_mtx */
  err = choldc(eff_d, dCorrMatrix, dSqCorrMatrix);

  if (err)
    goto FREE_RETURN;

  switch (sParams->eMCType) {
  case ABS: // uses ABS cube
    break;

  case SOBOLBM: // uses SOBOL with Box-Muller algorithm
    break;

  case RANDOM_GAUSS:
    break;

  default: // by default uses SOBOL with Normal Inversion

    GetSobolMatrix(sParams->lNumPaths, eff_d, GaussSample);

    //////////  Step3: Generates the matrix of r.v. compatible with the Copula
    ///and the marginals
    for (i = 0; i < sParams->lNumPaths; i++) {
      /* Get the uncorrelated gaussians */
      for (j = 0; j < eff_d; j++) {
        GaussSample[i][j] = inv_cumnorm_fast(GaussSample[i][j]);
      }

      /* Correlate the gaussians */
      for (j = 0; j < eff_d; j++) {
        dBrownian[j] = 0.0;

        for (k = 0; k < eff_d; k++) {
          dBrownian[j] += dSqCorrMatrix[j][k] * GaussSample[i][k];
        }
      }

      /* Get the true variables */
      for (j = 0; j < iNbFwd; j++) {
        XtSample[j][i] =
            dCoef1_1[j] * dBrownian[j] +
            dBrownian[j + iNbFwd] * (dCoef2_1[j] + dCross[j] * dBrownian[j]);
      }
    }

    /* Calculates the distributions */
    for (j = 0; j < iNbFwd; j++) {
      /* Calculates the empirical distribution */
      for (i = 0; i < sParams->lNumPaths; i++) {
        dSortMatrix[0][i] = i;
      }

      memcpy(dSortMatrix[1], XtSample[j], sParams->lNumPaths * sizeof(double));

      num_f_sort_matrix(sParams->lNumPaths, 2, 1, dSortMatrix);

      IntermRes = ProbaStep;

      for (i = 0; i < sParams->lNumPaths; i++) {
        XtSample[j][(long)(dSortMatrix[0][i])] = IntermRes;

        IntermRes += ProbaStep;
      }
    }

    /* Calculates the quadrat */

    memset(dResult, 0, iNbFract * sizeof(double));

    for (i = 3; i < sParams->lNumPaths; i++) {
      for (j = 0; j < iNbFract; j++) {
        dCashFlow = 1.0;

        for (k = 0; k < iNbFwd; k++) {
          if (dWeights[k] * XtSample[k][i] < dWeights[k] * dFract[j]) {
            dCashFlow = 0.0;
            break;
          }
        }

        dResult[j] += (dCashFlow - dResult[j]) / (i - 2);
      }
    }

    break;
  }

FREE_RETURN:

  if (GaussSample)
    free_dmatrix(GaussSample, 0, sParams->lNumPaths, 0, eff_d - 1);
  if (XtSample)
    free_dmatrix(XtSample, 0, iNbFwd - 1, 0, sParams->lNumPaths);
  if (dSqCorrMatrix)
    free_dmatrix(dSqCorrMatrix, 0, eff_d - 1, 0, eff_d - 1);
  if (dBrownian)
    free(dBrownian);
  if (dSortMatrix)
    free_dmatrix(dSortMatrix, 0, 1, 0, sParams->lNumPaths);

  if (dConst)
    free(dConst);
  if (dCoef1_1)
    free(dCoef1_1);
  if (dCoef2_1)
    free(dCoef2_1);
  if (dCross)
    free(dCross);

  return err;
}

/* Pricing of (dWeigths[0]*dFwds[0] + dMargin) * max(dWeights[1->n-1] *
 * dFwds[1->n-1] - K  , 0) */
Err copula_sv_GearedOption_SABR(int nFwds, double *dFwds, double *dWeights,
                                double dMargin, int iNbStrikes, double *dStrike,
                                double *dVols, double *dAlpha, double *dBeta,
                                double *dRho, double **dCorrelation,
                                double dMaturity, SrtCallPutType eCallPut,
                                SrtDiffusionType eTypeInput,
                                COPULASV_PARAMS sParams, double *dPremium) {
  Err err = NULL;

  double *pdFwdsAtMaturity = NULL, **ppdPoints = NULL,
         **ppdCumProbaValue = NULL, **dCorrCube = NULL, *dAvg = NULL;

  double *dShiftSV = NULL, *dSigmaSV = NULL, *dAlphaSV = NULL, *dRhoSV = NULL,
         *dVar0SV = NULL, *dLambdaSV = NULL;

  long lNbPoints[10];

  double dApproxStdev, dCoeffAdjustStrike, dCoeff, dAdjStrike, dTemp1, dTemp2,
      dGearFwd, dBasket;
  double dCashFlow, dSqMat, cp;
  double n_StdDev, PrMin, RaccRatio;
  double dSigmaBeta;

  long i, j, k, kfail, n_eff, nNumPoints;

  SrtCallPutType CallPut;
  SrtDiffusionType log_norm;

  /* Constant Initialisation */
  cp = (eCallPut == SRT_CALL) ? 1 : -1;
  dSqMat = sqrt(dMaturity);
  CallPut = SRT_CALL;
  log_norm = SRT_LOGNORMAL;
  nNumPoints = sParams->lNbPoints;

  for (i = 0; i < nFwds; i++) {
    lNbPoints[i] = 2 * nNumPoints;
  }

  /* Memory allocation */
  pdFwdsAtMaturity = calloc(nFwds, sizeof(double));
  dAvg = calloc(nFwds, sizeof(double));

  dShiftSV = calloc(nFwds, sizeof(double));
  dSigmaSV = calloc(nFwds, sizeof(double));
  dAlphaSV = calloc(nFwds, sizeof(double));
  dRhoSV = calloc(nFwds, sizeof(double));
  dVar0SV = calloc(nFwds, sizeof(double));
  dLambdaSV = calloc(nFwds, sizeof(double));

  ppdPoints = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);
  ppdCumProbaValue = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);
  dCorrCube = dmatrix(0, sParams->lNumPaths, 0, nFwds - 1);

  if (!pdFwdsAtMaturity || !dAvg || !dShiftSV || !dSigmaSV || !dAlphaSV ||
      !dRhoSV || !dVar0SV || !dLambdaSV || !ppdPoints || !ppdPoints ||
      !ppdCumProbaValue || !dCorrCube) {
    err = "Memory allocation faillure in copula_sv_basket_SABR";
    goto FREE_RETURN;
  }

  /* Define the points */
  /* Compute the cumulative proba = 1 + d_K Call(K) */

  for (i = 0; i < nFwds; i++) {
    /* Sigma Beta */
    srt_f_optsarbvol(dFwds[i], dFwds[i], dMaturity, dVols[i], dAlpha[i],
                     dBeta[i], dRho[i], eTypeInput, SRT_BETAVOL, &dSigmaBeta);

    err = copula_sv_get_SABR_cumulative(
        dMaturity, dFwds[i], dSigmaBeta, dAlpha[i], dBeta[i], dRho[i],
        sParams->iAdaptGrid, sParams->lAdaptNbPoints, 0, 2 * nNumPoints + 1,
        ppdPoints[i], ppdCumProbaValue[i]);

    if (err)
      goto FREE_RETURN;

    if (sParams->iUseOldGridMethod) {
      // set the Approximation for the Stdev //
      dApproxStdev = dSigmaBeta * pow(dFwds[i], dBeta[i] - 1.0) * dSqMat;

      // For the shift of the Strike
      dCoeffAdjustStrike = exp(dApproxStdev / nNumPoints);
      dCoeff =
          dCoeffAdjustStrike / (dCoeffAdjustStrike * dCoeffAdjustStrike - 1.0);

      // set the value at the boundary n_StdDev number of standard deviations
      n_StdDev = 5.0 * dSqMat;
      ppdPoints[i][0] = dFwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev -
                                       n_StdDev * dApproxStdev);
      ppdPoints[i][2 * nNumPoints] =
          dFwds[i] *
          exp(-0.5 * dApproxStdev * dApproxStdev + n_StdDev * dApproxStdev);

      ppdCumProbaValue[i][0] = 0.0;

      for (j = 1; j <= 2 * nNumPoints; j++) {
        ppdPoints[i][j] =
            dFwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev +
                           (j - nNumPoints) * 4.0 / nNumPoints * dApproxStdev);

        // Compute the Cum Proba associated to this points P(S < pdPoints)
        dAdjStrike = ppdPoints[i][j] * dCoeffAdjustStrike;

        dTemp1 = srt_f_optblkschbetastochquick(
            dFwds[i], dAdjStrike, dMaturity, dSigmaBeta, dAlpha[i], dBeta[i],
            dRho[i], 1.0, SRT_BETAVOL, log_norm, CallPut, PREMIUM);

        dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

        dTemp2 = srt_f_optblkschbetastochquick(
            dFwds[i], dAdjStrike, dMaturity, dSigmaBeta, dAlpha[i], dBeta[i],
            dRho[i], 1.0, SRT_BETAVOL, log_norm, CallPut, PREMIUM);

        ppdCumProbaValue[i][j] =
            1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][j];
      }

      //// find the point where eventually SABR function fails

      kfail = 0;

      if (dFwds[i] >= 0.036) {
        n_eff = 0;
      } else if (dFwds[i] >= 0.02 && dFwds[i] < 0.036) {
        n_eff = nNumPoints;
      } else if (dFwds[i] < 0.02) {
        n_eff = 2 * nNumPoints - 5;
      }

      for (j = 0; j < n_eff; j++) {
        if (ppdCumProbaValue[i][j + 1] <= ppdCumProbaValue[i][j]) {
          kfail = j + 1;
          PrMin = ppdCumProbaValue[i][j + 1];
        }
      }

      //// then use a lower vovol

      switch (kfail) {
      case 0:

        break;

      default:

        for (k = 1; k <= kfail + 1; k++) {
          dAdjStrike = ppdPoints[i][k] * dCoeffAdjustStrike;

          dTemp1 = srt_f_optblkschbetastochquick(
              dFwds[i], dAdjStrike, dMaturity, dSigmaBeta, 0.0, dBeta[i], 0.0,
              1.0, SRT_BETAVOL, log_norm, CallPut, PREMIUM);

          dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

          dTemp2 = srt_f_optblkschbetastochquick(
              dFwds[i], dAdjStrike, dMaturity, dSigmaBeta, 0.0, dBeta[i], 0.0,
              1.0, SRT_BETAVOL, log_norm, CallPut, PREMIUM);

          ppdCumProbaValue[i][k] =
              1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][k];
        }

        RaccRatio = PrMin / ppdCumProbaValue[i][kfail + 1];

        for (k = 1; k <= kfail + 1; k++) {
          ppdCumProbaValue[i][k] *= RaccRatio;
        }

        break;
      }
    }

    for (j = 1; j <= 2 * nNumPoints; j++) {
      ppdCumProbaValue[i][j] /= ppdCumProbaValue[i][2 * nNumPoints];
    }
  }

  /* Set the variables and the transformation */
  for (i = 0; i < nFwds; i++) {
    err = copula_sv_from_sabr_to_heston(
        dMaturity, dFwds[i], dVols[i], dAlpha[i], dBeta[i], dRho[i], eTypeInput,
        &(dSigmaSV[i]), &(dAlphaSV[i]), &(dShiftSV[i]), &(dRhoSV[i]),
        &(dVar0SV[i]), &(dVar0SV[i]), &(dLambdaSV[i]));

    if (err)
      goto FREE_RETURN;
  }

  err = copula_sv_numer(nFwds, ppdPoints, ppdCumProbaValue, &lNbPoints[0],
                        dMaturity, dFwds, dShiftSV, dSigmaSV, dAlphaSV, dRhoSV,
                        dVar0SV, dVar0SV, dLambdaSV, dCorrelation, sParams,
                        dCorrCube);

  if (err)
    goto FREE_RETURN;

  /* Compute the Price and the normal vols via a Monte Carlo Methods */

  /* Initialisation */
  memset(dPremium, 0, iNbStrikes * sizeof(double));
  memset(dAvg, 0, nFwds * sizeof(double));

  for (i = 3; i < sParams->lNumPaths; i++) {
    dBasket = 0.0;

    pdFwdsAtMaturity[0] = dCorrCube[i][0];
    dAvg[0] += (pdFwdsAtMaturity[0] - dAvg[0]) / (i - 2);

    for (j = 1; j < nFwds; j++) {
      pdFwdsAtMaturity[j] = dCorrCube[i][j];

      dAvg[j] += (pdFwdsAtMaturity[j] - dAvg[j]) / (i - 2);
      dBasket += dWeights[j] * pdFwdsAtMaturity[j] * dFwds[j] / dAvg[j];
    }

    dGearFwd = dWeights[0] * pdFwdsAtMaturity[0] * dFwds[0] / dAvg[0] + dMargin;

    for (j = 0; j < iNbStrikes; j++) {
      dCashFlow = cp * (dBasket - dStrike[j]);
      dPremium[j] += (dGearFwd * DMAX(dCashFlow, 0.0) - dPremium[j]) / (i - 2);
    }
  }

FREE_RETURN:

  if (pdFwdsAtMaturity)
    free(pdFwdsAtMaturity);
  if (dAvg)
    free(dAvg);

  if (dShiftSV)
    free(dShiftSV);
  if (dSigmaSV)
    free(dSigmaSV);
  if (dAlphaSV)
    free(dAlphaSV);
  if (dRhoSV)
    free(dRhoSV);
  if (dVar0SV)
    free(dVar0SV);
  if (dLambdaSV)
    free(dLambdaSV);

  if (ppdPoints)
    free_dmatrix(ppdPoints, 0, nFwds - 1, 0, 2 * nNumPoints);
  if (ppdCumProbaValue)
    free_dmatrix(ppdCumProbaValue, 0, nFwds - 1, 0, 2 * nNumPoints);
  if (dCorrCube)
    free_dmatrix(dCorrCube, 0, sParams->lNumPaths, 0, nFwds - 1);

  return err;
}

/* Monte Carlo to check results */
Err copula_sv_basket_SABR_MC(int nFwds, double *dFwds, double *dWeights,
                             int iNbStrikes, double *dStrike, double *dVols,
                             double *dAlpha, double *dBeta, double *dRho,
                             double **dCorrelation, double dMaturity,
                             SrtCallPutType eCallPut,
                             SrtDiffusionType eTypeInput,
                             COPULASV_PARAMS sParams, double *dPremium) {
  Err err = NULL;

  int i, j, k, l, eff_d;
  long lNumSteps;
  double dt, sqdt;

  double **dCorrMatrix = NULL, **dSqCorrMatrix = NULL, *dXtSample = NULL,
         *dIndBrownian = NULL, *dCorBrownian = NULL, *dCoefVol = NULL,
         *dStartVol = NULL, *dCoefVoVol = NULL;

  double cp, dBasket;

  long seed = -123456789;

  /* Memory allocation */
  eff_d = 2 * nFwds;
  dCorrMatrix = dmatrix(0, eff_d - 1, 0, eff_d - 1);
  dSqCorrMatrix = dmatrix(0, eff_d - 1, 0, eff_d - 1);
  dXtSample = calloc(eff_d, sizeof(double));
  dIndBrownian = calloc(eff_d, sizeof(double));
  dCorBrownian = calloc(eff_d, sizeof(double));
  dCoefVol = calloc(nFwds, sizeof(double));
  dStartVol = calloc(nFwds, sizeof(double));
  dCoefVoVol = calloc(nFwds, sizeof(double));

  if (!dSqCorrMatrix || !dXtSample || !dIndBrownian || !dCorBrownian ||
      !dCoefVol || !dStartVol || !dCoefVoVol) {
    err = "Memory allocation faillure in copula_sv_basket_SABR_MC";
    goto FREE_RETURN;
  }

  /* Fill the Corr Matrix */
  if (1) {
    for (i = 0; i < nFwds; i++) {
      for (j = 0; j < nFwds; j++) {
        dCorrMatrix[i][j] = dCorrelation[i][j];
      }
    }
  } else {
    for (i = 0; i < nFwds; i++) {
      dCorrMatrix[i][i] = 1.0;

      /* Rates / Rates */
      for (j = i + 1; j < nFwds; j++) {
        dCorrMatrix[i][j] = dCorrelation[i][j];
        dCorrMatrix[j][i] = dCorrelation[i][j];
      }

      for (j = 0; j < nFwds; j++) {
        /* Rates Vol */
        dCorrMatrix[i][nFwds + j] = dRho[i];
        dCorrMatrix[nFwds + j][i] = dRho[i];

        /* Vol Vol */
        dCorrMatrix[nFwds + i][nFwds + j] = 0.999;
        dCorrMatrix[nFwds + j][nFwds + i] = 0.999;
      }
    }
  }

  /* performs Chol decomposition of the correlation matrix corr_mtx */
  err = choldc(eff_d, dCorrMatrix, dSqCorrMatrix);

  if (err)
    goto FREE_RETURN;

  lNumSteps = (long)(dMaturity / sParams->dMinTime);
  dt = dMaturity / (lNumSteps * 1.0);
  sqdt = sqrt(dt);

  for (i = 0; i < nFwds; i++) {
    srt_f_optsarbvol(dFwds[i], dFwds[i], dMaturity, dVols[i], dAlpha[i],
                     dBeta[i], dRho[i], eTypeInput, SRT_BETAVOL,
                     &(dStartVol[i]));

    dCoefVol[i] = sqdt;
    dCoefVoVol[i] = dAlpha[i] * sqdt;
  }

  cp = (eCallPut == SRT_CALL) ? 1 : -1;

  memset(dPremium, 0, iNbStrikes * sizeof(double));

  for (i = 0; i < sParams->lNumPaths; i++) {
    /* Initialisation */
    for (k = 0; k < nFwds; k++) {
      dXtSample[k] = dFwds[k];
      dXtSample[nFwds + k] = dStartVol[k];
    }

    for (j = 0; j < lNumSteps; j++) {
      /* Get the brownians */
      for (k = 0; k < eff_d; k++) {
        dIndBrownian[k] = inv_cumnorm_fast(uniform(&seed));
      }

      for (k = 0; k < eff_d; k++) {
        dCorBrownian[k] = 0.0;

        for (l = 0; l < eff_d; l++) {
          dCorBrownian[k] += dSqCorrMatrix[k][l] * dIndBrownian[l];
        }
      }

      /* Move the Fwds and Vols */
      for (k = 0; k < nFwds; k++) {
        dXtSample[k] =
            max(dXtSample[k] + dXtSample[nFwds + k] *
                                   exp(dBeta[k] * log(dXtSample[k])) * sqdt *
                                   dCorBrownian[k],
                0.0001);
        dXtSample[nFwds + k] =
            max(dXtSample[nFwds + k] *
                    (1.0 + dCoefVoVol[k] * dCorBrownian[nFwds + k]),
                0.0001);
      }
    }

    /* Payoff */
    dBasket = 0.0;

    for (k = 0; k < nFwds; k++) {
      dBasket += dWeights[k] * dXtSample[k];
    }

    for (k = 0; k < iNbStrikes; k++) {
      dPremium[k] +=
          DMAX(cp * (dBasket - dStrike[k]), 0.0) / sParams->lNumPaths;
    }
  }

FREE_RETURN:

  if (dCorrMatrix)
    free_dmatrix(dCorrMatrix, 0, eff_d - 1, 0, eff_d - 1);
  if (dSqCorrMatrix)
    free_dmatrix(dSqCorrMatrix, 0, eff_d - 1, 0, eff_d - 1);
  if (dXtSample)
    free(dXtSample);
  if (dIndBrownian)
    free(dIndBrownian);
  if (dCorBrownian)
    free(dCorBrownian);
  if (dCoefVol)
    free(dCoefVol);
  if (dStartVol)
    free(dStartVol);
  if (dCoefVoVol)
    free(dCoefVoVol);

  return err;
}

Err copula_sv_quadrant_dependence_SABR(
    int nFwds, double *dFwds, double *dWeights, int iNbFract, double *dFract,
    double *dVols, double *dAlpha, double *dBeta, double *dRho,
    double **dCorrelation, double dMaturity, SrtDiffusionType eTypeInput,
    COPULASV_PARAMS sParams, double *dResult) {
  Err err = NULL;

  double *dShiftSV = NULL, *dSigmaSV = NULL, *dAlphaSV = NULL, *dRhoSV = NULL,
         *dVar0SV = NULL, *dLambdaSV = NULL;

  long i;

  /* Memory allocation */
  dShiftSV = calloc(nFwds, sizeof(double));
  dSigmaSV = calloc(nFwds, sizeof(double));
  dAlphaSV = calloc(nFwds, sizeof(double));
  dRhoSV = calloc(nFwds, sizeof(double));
  dVar0SV = calloc(nFwds, sizeof(double));
  dLambdaSV = calloc(nFwds, sizeof(double));

  if (!dShiftSV || !dSigmaSV || !dAlphaSV || !dRhoSV || !dVar0SV ||
      !dLambdaSV) {
    err = "Memory allocation faillure in copula_sv_quadrant_dependence_SABR";
    goto FREE_RETURN;
  }

  /* Set the variables and the transformation */
  for (i = 0; i < nFwds; i++) {
    err = copula_sv_from_sabr_to_heston(
        dMaturity, dFwds[i], dVols[i], dAlpha[i], dBeta[i], dRho[i], eTypeInput,
        &(dSigmaSV[i]), &(dAlphaSV[i]), &(dShiftSV[i]), &(dRhoSV[i]),
        &(dVar0SV[i]), &(dVar0SV[i]), &(dLambdaSV[i]));

    if (err)
      goto FREE_RETURN;
  }

  err = copula_sv_quadrat(nFwds, dWeights, iNbFract, dFract, dMaturity, dFwds,
                          dShiftSV, dSigmaSV, dAlphaSV, dRhoSV, dVar0SV,
                          dVar0SV, dLambdaSV, dCorrelation, sParams, dResult);

  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (dShiftSV)
    free(dShiftSV);
  if (dSigmaSV)
    free(dSigmaSV);
  if (dAlphaSV)
    free(dAlphaSV);
  if (dRhoSV)
    free(dRhoSV);
  if (dVar0SV)
    free(dVar0SV);
  if (dLambdaSV)
    free(dLambdaSV);

  return err;
}

Err copula_student_quadrant_dependence(int iNbFwd, double *dWeights,
                                       int iNbFract, double *dFract,
                                       double **dCorrMatrix,
                                       COPULASV_PARAMS sParams,
                                       double *dResult) {
  long i, idum = -8935807;
  int j, k, eff_d, total_d;
  double **GaussSample = NULL, **XtSample = NULL, *dBrownian = NULL,
         **dSortMatrix = NULL, **dSqCorrMatrix = NULL;

  double chi2, IntermRes, dCashFlow;

  Err err = NULL;

  /* Memory allocation */
  total_d = iNbFwd + sParams->iStudDegree;
  eff_d = iNbFwd;

  GaussSample = dmatrix(0, sParams->lNumPaths, 0, total_d - 1);
  XtSample = dmatrix(0, eff_d - 1, 0, sParams->lNumPaths);
  dSqCorrMatrix = dmatrix(0, eff_d - 1, 0, eff_d - 1);
  dBrownian = calloc(total_d, sizeof(double));

  if (!GaussSample || !XtSample || !dSqCorrMatrix || !dBrownian) {
    err = "Memory allocation faillure in copula_student_quadrant_dependence";
    goto FREE_RETURN;
  }

  /* performs Chol decomposition of the correlation matrix corr_mtx */
  err = choldc(eff_d, dCorrMatrix, dSqCorrMatrix);

  if (err)
    goto FREE_RETURN;

  switch (sParams->eMCType) {
  case ABS: // uses ABS cube
    break;

  case SOBOLBM: // uses SOBOL with Box-Muller algorithm
    break;

  case RANDOM_GAUSS:
    break;

  default: // by default uses SOBOL with Normal Inversion

    GetSobolMatrix(sParams->lNumPaths, total_d, GaussSample);

    //////////  Step3: Generates the matrix of r.v. compatible with the Copula
    ///and the marginals
    for (i = 3; i < sParams->lNumPaths; i++) {
      /* Get the uncorrelated gaussians */
      for (j = 0; j < total_d; j++) {
        GaussSample[i][j] = inv_cumnorm_fast(GaussSample[i][j]);
      }

      /* Correlate the gaussians */
      for (j = 0; j < eff_d; j++) {
        dBrownian[j] = 0.0;

        for (k = 0; k < eff_d; k++) {
          dBrownian[j] += dSqCorrMatrix[j][k] * GaussSample[i][k];
        }
      }

      /* Compute the draw of chi2 */
      chi2 = 0.0;

      for (j = eff_d; j < total_d; j++) {
        chi2 += GaussSample[i][j] * GaussSample[i][j];
      }

      for (j = 0; j < eff_d; j++) {
        if (sParams->iStudDegree > 0) {
          if (i > 0) {
            IntermRes = sqrt(sParams->iStudDegree / chi2) * dBrownian[j];
          } else {
            IntermRes = 0.0;
          }
        } else {
          IntermRes = dBrownian[j];
        }

        XtSample[j][i] = Student_Dis(sParams->iStudDegree, IntermRes);
      }
    }

    /* Calculates the quadrat */
    memset(dResult, 0, iNbFract * sizeof(double));

    for (i = 3; i < sParams->lNumPaths; i++) {
      for (j = 0; j < iNbFract; j++) {
        dCashFlow = 1.0;

        for (k = 0; k < iNbFwd; k++) {
          if (dWeights[k] * XtSample[k][i] < dWeights[k] * dFract[j]) {
            dCashFlow = 0.0;
            break;
          }
        }

        dResult[j] += (dCashFlow - dResult[j]) / (i - 2);
      }
    }

    break;
  }

FREE_RETURN:

  if (GaussSample)
    free_dmatrix(GaussSample, 0, sParams->lNumPaths, 0, total_d - 1);
  if (XtSample)
    free_dmatrix(XtSample, 0, eff_d - 1, 0, sParams->lNumPaths);
  if (dSqCorrMatrix)
    free_dmatrix(dSqCorrMatrix, 0, eff_d - 1, 0, eff_d - 1);
  if (dBrownian)
    free(dBrownian);

  return err;
}

Err copula_sv_get_implied_linear_correl(int iNbFwd, int iIndex1, int iIndex2,
                                        double *dCoef1, double *dCoef2,
                                        double *dCross, double **dCorrelation,
                                        double *dRho) {
  double **dNewCorrel = NULL, **dCholesky = NULL;

  int i, j;
  double dVar1, dVar2, dCovar;

  Err err = NULL;

  /* Memory allocation */
  dNewCorrel = dmatrix(0, 3, 0, 3);
  dCholesky = dmatrix(0, 3, 0, 3);

  if (!dNewCorrel || !dCholesky) {
    err = "Memory allocation faillure in copula_sv_get_implied_linear_correl";
    goto FREE_RETURN;
  }

  /* Fill the new correl */
  dNewCorrel[0][1] = dCorrelation[iIndex2 + iNbFwd][iIndex1 + iNbFwd];
  dNewCorrel[0][2] = dCorrelation[iIndex2 + iNbFwd][iIndex2];
  dNewCorrel[0][3] = dCorrelation[iIndex2 + iNbFwd][iIndex1];
  dNewCorrel[1][2] = dCorrelation[iIndex1 + iNbFwd][iIndex2];
  dNewCorrel[1][3] = dCorrelation[iIndex1 + iNbFwd][iIndex1];
  dNewCorrel[2][3] = dCorrelation[iIndex1][iIndex2];

  for (i = 0; i < 4; i++) {
    dNewCorrel[i][i] = 1.0;

    for (j = i; j < 4; j++) {
      dNewCorrel[j][i] = dNewCorrel[i][j];
    }
  }

  /* Get the Cholesky */
  err = choldc(4, dNewCorrel, dCholesky);

  if (err)
    goto FREE_RETURN;

  /* Calculates the variances and covariances */
  dVar1 =
      pow(dCoef1[iIndex1] * dCholesky[3][0] + dCoef2[iIndex1] * dCholesky[1][0],
          2.0) +
      pow(dCoef1[iIndex1] * dCholesky[3][1] + dCoef2[iIndex1] * dCholesky[1][1],
          2.0) +
      pow(dCoef1[iIndex1] * dCholesky[3][2], 2.0) +
      pow(dCoef1[iIndex1] * dCholesky[3][3], 2.0) +
      3.0 * pow(dCross[iIndex1] * dCholesky[3][0] * dCholesky[1][0], 2.0) +
      3.0 * pow(dCross[iIndex1] * dCholesky[3][1] * dCholesky[1][1], 2.0);

  dVar2 =
      pow(dCoef1[iIndex2] * dCholesky[2][0] + dCoef2[iIndex2] * dCholesky[0][0],
          2.0) +
      pow(dCoef1[iIndex2] * dCholesky[2][1], 2.0) +
      pow(dCoef1[iIndex2] * dCholesky[2][2], 2.0) +
      3.0 * pow(dCross[iIndex2] * dCholesky[2][0] * dCholesky[0][0], 2.0);

  dCovar =
      (dCoef1[iIndex1] * dCholesky[3][0] + dCoef2[iIndex1] * dCholesky[1][0]) *
          (dCoef1[iIndex2] * dCholesky[2][0] +
           dCoef2[iIndex2] * dCholesky[0][0]) +
      (dCoef1[iIndex1] * dCholesky[3][1] + dCoef2[iIndex1] * dCholesky[1][1]) *
          (dCoef1[iIndex2] * dCholesky[2][1]) +
      (dCoef1[iIndex1] * dCholesky[3][2]) *
          (dCoef1[iIndex2] * dCholesky[2][2]) +
      3.0 * dCross[iIndex1] * dCholesky[3][0] * dCholesky[1][0] *
          dCross[iIndex2] * dCholesky[2][0] * dCholesky[0][0];

  if (dVar1 < 0.0 || dVar2 < 0.0) {
    err = "Correlation Matrix not definite positive in "
          "copula_sv_get_implied_linear_correl";
    goto FREE_RETURN;
  }

  *dRho = dCovar / sqrt(dVar1 * dVar2);

FREE_RETURN:

  if (dNewCorrel)
    free_dmatrix(dNewCorrel, 0, 3, 0, 3);
  if (dCholesky)
    free_dmatrix(dCholesky, 0, 3, 0, 3);

  return err;
}