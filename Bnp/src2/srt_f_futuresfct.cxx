/* =======================================================================================

   FILENAME :       srt_f_futuresfct.cxx

   AUTHOR:          C. Drozo

   DATE:            August-September 2000
                    This version: 21/09/2000

   PURPOSE:         Various functions used to compute the convexity adjustments
                    between forward and futures contracts
   =======================================================================================
 */

// Header files & other declarations & definitions

#include "math.h"
#include "num_h_interp.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_futuresfct.h"
#include "srtaccess.h"
#include "swp_h_all.h"

#define FREE_ARG char *
#define NR_END 1
static void free_dvector1(double **v, long nl) {
  free((FREE_ARG)(*v + nl - NR_END));
  *v = NULL;
}
#undef NR_END
#undef FREE_ARG

/*=================================================================
    Interpolate the volatilities through the futures options one
    and complete the vectors with the caplet ones for longer expiries
-------------------------------------------------------------------*/
Err eBuildVols(long lToday) {
  int i;
  int iMax = nbTenors2;
  double *pdX;
  double *pdY;
  double q;

  dVolStruct2 = dvector(0, iNumFut - 1);
  pdX = dvector(0, iMax - 1);
  pdY = dvector(0, iMax - 1);
  for (i = 0; i < iMax; ++i) {
    pdX[i] = (ppdOptFutVol[i][0] - lToday) / 365.0;
    pdY[i] = ppdOptFutVol[i][1];
  }

  for (i = 0; i < iNumFut; ++i) {
    if (dDates[i] >= pdX[0] && dDates[i] <= pdX[iMax - 1]) {
      dVolStruct2[i] = interp(pdX, pdY, iMax, dDates[i], 1, &q);
    } else {
      dVolStruct2[i] = dVolStruct[i];
    }
  }
  free_dvector1(&pdX, 0);
  free_dvector1(&pdY, 0);
  return NULL;
}

/*=================================================================
    Compute the correlation between two LIBORs        , depending on their
    relative time to maturity. Exponential decreasing law
    beta can be fitted
-------------------------------------------------------------------*/
double dCorrel(double dT, double beta) {
  double dRhoMin = beta; // correlation to fit
  double dTmax = 10.0;
  double e = exp(1.);
  //    return dRhoMin+(1.0-dRhoMin)*exp(-dT/dTmax);
  //    return 1.+(dRhoMin-1.0)*dT/dTmax;
  return 1. / (1. - e) *
         (1. - dRhoMin * e + (dRhoMin - 1.0) * exp(1. - dT / dTmax));
}

/*=================================================================
    Compute volatility integral of the type:
    Int (vol(u        ,TP)vol(u        ,Tk)du) between 0 and Tk
    Use of sliding stationary volatility structure or of
    forward volatility stored in **ppdFwdVol
-------------------------------------------------------------------*/
double dCrossIntegral(int indexP, int indexK, double *dPeriod) {
  // fixed index K
  int i, j;
  double dInteg = 0.0;
  double dTmpVol;
  double *pdVolk;
  double *pdVolp;

  // sliding stationary volatility structure
  if (!(ppdFwdVol) && !(ppdOptFutVol)) {
    for (i = 0; i < indexK; ++i) {
      dInteg +=
          dVolStruct[indexK - i - 1] * dVolStruct[indexP - i - 1] * dPeriod[i];
    }
  } else if (ppdFwdVol) {
    // forward volatility reconstruction
    pdVolk = dvector(0, indexK - 1);
    for (j = 0; j < indexK; ++j) {
      dTmpVol = 0.;
      for (i = j; i < indexK; ++i)
        dTmpVol += pow(ppdFwdVol[i][indexK], 2) * dPeriod[i];
      pdVolk[j] = sqrt(dTmpVol * 365.0 / (dDates[indexK] - dDates[j]));
    }

    pdVolp = dvector(0, indexP - 1);
    for (j = 0; j < indexP; ++j) {
      dTmpVol = 0.;
      for (i = j; i < indexP; ++i)
        dTmpVol += pow(ppdFwdVol[i][indexP], 2) * dPeriod[i];
      pdVolp[j] = sqrt(dTmpVol * 365.0 / (dDates[indexP] - dDates[j]));
    }

    for (i = 0; i < indexK; ++i)
      dInteg += pdVolk[i] * pdVolp[i] * dPeriod[i];
    free_dvector(pdVolp, 0, indexP - 1);
    pdVolp = NULL;
    free_dvector(pdVolk, 0, indexK - 1);
    pdVolk = NULL;
  } else if (ppdOptFutVol) {
    for (i = 0; i < indexK; ++i) {
      dInteg += dVolStruct2[indexK - i - 1] * dVolStruct2[indexP - i - 1] *
                dPeriod[i];
    }
  }

  return dInteg;
}

/*=================================================================
    Compute the drift induced by the change of probability with
    the market model
-------------------------------------------------------------------*/
double dDrift(int iNumFut, int indexK, double *dPeriod, double *dFwd,
              double beta) {
  // fixed index K
  int p;
  double drift = 0.0;
  double dLevel;

  for (p = indexK + 1; p <= iNumFut; ++p) {
    dLevel = dPeriod[p] * dFwd[p - 1];
    drift += dLevel / (1. + dLevel) * dCrossIntegral(p, indexK, dPeriod);
  }
  return drift;
}

/*=================================================================
    Compute the variance under PTj+1 of the LIBOR L(Tk        ,Tk)
-------------------------------------------------------------------*/
double dVar(int indexK, double *dPeriod) {
  int p;
  double dInteg = 0.;

  if (!ppdOptFutVol) {
    for (p = 0; p < indexK; ++p) {
      dInteg +=
          dVolStruct[indexK - p - 1] * dVolStruct[indexK - p - 1] * dPeriod[p];
    }
  } else {
    for (p = 0; p < indexK; ++p) {
      dInteg += dVolStruct2[indexK - p - 1] * dVolStruct2[indexK - p - 1] *
                dPeriod[p];
    }
  }
  return exp(dInteg) - 1.;
}
//=======================================================================
//      Forward-Futures convexity adjustment computation
//      There are 4 models:
//          -BGM:           future_Rutkowski
//          -HJM:           future_HJM
//          -Ho & Lee:      future_HL
//          -Hull & White:  future_Bloom
/*=======================================================================
    Forward-Futures adjustments computed according to a market model
    See Marek Rutkowski article
    2nd order approximation
------------------------------------------------------------------------*/
Err future_Rutkowski(long index, double *dPeriod, double *dParams,
                     double *pdFuturesRate) {
  double tmpProd = 1., tmpDen = 1., tmpSum = 0., dVarj, dFutRate;
  double beta;
  int i, k, p;
  /*beta introduced as a parameter to fit the minimal correlation between 2
    Libors depending on the time to maturity*/
  if (!dParams) {
    beta = 1.2;
  } else {
    beta = dParams[0];
  }

  tmpProd = tmpDen = 1.0;

  for (k = 1; k < index; ++k) {
    tmpDen *= (1.0 + dPeriod[k] * dL0T[k - 1]);
  }

  dVarj = dVar(index, dPeriod);
  for (k = 1; k < index; ++k) {
    tmpProd += dPeriod[k] * dL0T[k - 1] *
               (exp(-dDrift(index, k, dPeriod, dL0T, beta)) +
                dCorrel(fabs(dDates[index] - dDates[k]), beta) *
                    sqrt(dVarj * dVar(k, dPeriod)));
  }
  for (i = 1; i < index; ++i)
    for (p = 1; p < index; ++p) {
      tmpSum += dL0T[i - 1] * dL0T[k - 1] * dPeriod[k] * dPeriod[i] *
                (exp(-dDrift(index, i, dPeriod, dL0T, beta) -
                     dDrift(index, p, dPeriod, dL0T, beta)) +
                 dCorrel(fabs(dDates[i] - dDates[p]), beta) *
                     sqrt(dVar(i, dPeriod) * dVar(p, dPeriod)));
    }
  for (k = 1; k < index; ++k) {
    tmpSum -= dPeriod[k] * dPeriod[k] * dL0T[k - 1] * dL0T[k - 1] *
              exp(-2 * dDrift(index, k, dPeriod, dL0T, beta) +
                  dCrossIntegral(k, k, dPeriod));
  }

  dFutRate = tmpProd + 0.5 * tmpSum;
  dFutRate *= dL0T[index - 1] / tmpDen;
  *pdFuturesRate = dFutRate;

  return NULL;
}

/*=======================================================================
    Forward-Futures adjustments computed according to HJM model
    modelisation of covariances and correlations
    2nd order approximation
------------------------------------------------------------------------*/
Err future_HJM(long index, double *dPeriod,
               double *dParams, // first: short term rate        , following:
                                // initial correlations
               double *pdFuturesRate) {
  double dFutRate;
  double dTau = 10.0;   // max time of dissipation
  double dShift = 0.08; // shift to get the short term rate form the IFR
  double dRhoI = 0.98;  // default value for the initial correlation law
  int i;
  double *pdFinVol;
  double *pdShRate;
  double *pdVol;
  double dCumTime;
  double dRho;
  double tmpSum1;
  double dBetaAdj, dBetaFut;
  double dG;
  double dShort;

  dShort = dL0T[index - 1] - dParams[0] / 100;
  if (sizeof(dParams) / sizeof(double) > 1) {
    dRho = dParams[1];
  } else {
    dRho = max(0.9, pow(dRhoI, iNumFut - 1));
  }

  pdFinVol = (double *)calloc(iNumFut, sizeof(double));
  pdShRate = (double *)calloc(iNumFut, sizeof(double));
  pdVol = (double *)calloc(iNumFut, sizeof(double));

  if (ppdOptFutVol) {
    for (i = 0; i < iNumFut; ++i)
      pdVol[i] = dVolStruct2[i];
  } else {
    for (i = 0; i < iNumFut; ++i)
      pdVol[i] = dVolStruct[i];
  }

  dCumTime = dDates[index] - dDates[0];

  pdShRate[0] = dL0T[0] - dShift / 100;
  pdFinVol[0] = pdVol[0] * dL0T[0];
  for (i = 1; i < index; ++i) {
    pdShRate[i] = (pdShRate[i - 1] * i + dFwd[i] - dShift / 100) / (i + 1);
    pdFinVol[i] =
        (pdFinVol[i - 1] * i + pdVol[i] * dFwd[i]) / (i + 1); // short term vol
  }

  for (i = 0; i < index; ++i)
    pdFinVol[i] /= pdShRate[i];

  tmpSum1 = dCrossIntegral(index, index, dPeriod);

  dBetaAdj = (1. + dL0T[index - 1] * dCoverage[index - 1]) *
             exp((pdVol[index - 1] * dFwd[index - 1]) *
                 (pdVol[index - 1] * dFwd[index - 1]) * dCumTime *
                 dCoverage[index - 1] * dCoverage[index - 1]);
  dG = 8. * dCumTime / dTau / 105 - 4. / 15;
  dG *= dCumTime / dTau;
  dG += 2. / 3;
  dG *= dCumTime * sqrt(dCumTime);
  dBetaFut = dBetaAdj * exp(dRho * (pdVol[index - 1] * dFwd[index - 1]) *
                            (pdFinVol[index - 1] * dShort) * sqrt(dCumTime) *
                            dG * dCoverage[index - 1]);
  dFutRate = (dBetaFut - 1.) / dCoverage[index - 1];

  free(pdFinVol);
  free(pdShRate);
  free(pdVol);

  *pdFuturesRate = dFutRate;

  return NULL;
}

/*=======================================================================
    Forward-Futures adjustments computed according to Hull and White model
    adjustments used by Bloomberg
------------------------------------------------------------------------*/
Err future_Bloom(long index, double *dPeriod, double *dParams,
                 double *pdFuturesRate) {
  double dZ;
  double dCoeff;
  double dExp1, dExp2, dExp3;
  double dFutRate;
  double dMeanRev;
  double dVolSr;

  if (!dParams) {
    dMeanRev = 0.03;
    dVolSr = 0.015;
  } else {
    dMeanRev = dParams[0];
    dVolSr = dParams[1];
  }

  dCoeff = pow(dVolSr, 2) / (4. * pow(dMeanRev, 3));
  dExp1 = 1.0 - exp(-dMeanRev * dCoverage[index - 1]);
  dExp2 = 1.0 - exp(-2.0 * dMeanRev * dDates[index]);
  dExp3 = 1.0 - exp(-dMeanRev * dDates[index]);

  dZ = dExp2 * dExp1;
  dZ += 2. * pow(dExp3, 2);
  dZ *= dCoeff * dExp1 / dCoverage[index - 1];

  dFutRate =
      (dB0T[index - 1] / dB0T[index] * exp(dZ) - 1.) / dCoverage[index - 1];

  *pdFuturesRate = dFutRate;

  return NULL;
}

/*=======================================================================
    Forward-Futures adjustments computed according to Ho and Lee model
------------------------------------------------------------------------*/
Err future_HL(long index, double *dPeriod, double *pdFuturesRate) {
  double dCumTime;
  double dVol;
  double dInteg, dAdjConv;
  double dFutRate;

  dCumTime = dDates[iNumFut] - dDates[0];

  dVol = ppdOptFutVol ? dVolStruct2[iNumFut - 1] : dVolStruct[iNumFut - 1];
  dInteg = 0.5 * dVol * dVol * dFwd[iNumFut - 1] * dFwd[iNumFut - 1] *
           pow(dCumTime, 2);
  dAdjConv = exp(dVol * dVol * dL0T[iNumFut - 1] * dL0T[iNumFut - 1] *
                 dCumTime * dCoverage[iNumFut - 1] * dCoverage[iNumFut - 1]);
  dAdjConv *= exp(dInteg * dCoverage[iNumFut - 1]);
  dFutRate =
      ((1. + dL0T[iNumFut - 1] * dCoverage[iNumFut - 1]) * dAdjConv - 1.) /
      dCoverage[iNumFut - 1];

  *pdFuturesRate = dFutRate;

  return NULL;
}

/*======================================================================
    Function used to optimize the parameters for the market model
------------------------------------------------------------------------*/
Err eFutFunc_Rut(double dIndex, double *pdParams, double *pdConvexity,
                 double *pdDerivative, int nParams) {
  Err err;
  double *dPeriod;
  double *dParams;
  double dFuturesRate;
  double dConv;
  double dEps;
  double dFp, dFm, dDer;
  int i;

  dParams = pdParams;

  dEps = 1.e-4;

  dPeriod = (double *)calloc((long)dIndex + 1, sizeof(double));
  for (i = 0; i <= (long)dIndex; i++)
    dPeriod[i] = dDates[i + 1] - dDates[i];

  err = future_Rutkowski((long)dIndex, dPeriod, dParams, &dFuturesRate);
  if (err) {
    free(dPeriod);
    return err;
  }
  dConv = 10000 * (dFuturesRate - dL0T[(long)dIndex - 1]);
  *pdConvexity = dConv;

  dParams[0] += dEps;
  err = future_Rutkowski((long)dIndex, dPeriod, dParams, &dFp);
  if (err) {
    free(dPeriod);
    return err;
  }
  dParams[0] -= 2. * dEps;
  err = future_Rutkowski((long)dIndex, dPeriod, dParams, &dFm);
  if (err) {
    free(dPeriod);
    return err;
  }

  dDer = 10000 * (dFp - dFm) / (2 * dEps);
  pdDerivative[1] = dDer;

  free(dPeriod);

  return NULL;
}

/*======================================================================
    Function used to optimized the parameters for the HJM model
------------------------------------------------------------------------*/
Err eFutFunc_HJM(double dIndex, double *pdParams, double *pdConvexity,
                 double *pdDerivative, int nParams) {
  Err err;
  double *dPeriod;
  double *dParams;
  double dFuturesRate;
  double dConv;
  double dEps;
  double dFp, dFm, dDer;
  int i;

  dParams = pdParams; // dParams[0]: short term rate        , dParams[1]: correl

  dEps = 1.e-4;

  dPeriod = (double *)calloc((long)dIndex + 1, sizeof(double));
  for (i = 0; i <= (long)dIndex; i++)
    dPeriod[i] = dDates[i + 1] - dDates[i];

  err = future_HJM((long)dIndex, dPeriod, dParams, &dFuturesRate);
  if (err) {
    free(dPeriod);
    return err;
  }
  dConv = 10000 * (dFuturesRate - dL0T[(long)dIndex - 1]);
  *pdConvexity = dConv;

  dParams[0] += dEps;
  err = future_HJM((long)dIndex, dPeriod, dParams, &dFp);
  if (err) {
    free(dPeriod);
    return err;
  }
  dParams[0] -= 2 * dEps;
  err = future_HJM((long)dIndex, dPeriod, dParams, &dFm);
  if (err) {
    free(dPeriod);
    return err;
  }
  dDer = 10000 * (dFp - dFm) / (2 * dEps);
  pdDerivative[1] = dDer;

  if (sizeof(dParams) / sizeof(double) > 1) {
    dParams[0] += dEps;
    dParams[1] += dEps;
    err = future_HJM((long)dIndex, dPeriod, dParams, &dFp);
    if (err) {
      free(dPeriod);
      return err;
    }
    dParams[1] -= 2 * dEps;
    err = future_HJM((long)dIndex, dPeriod, dParams, &dFm);
    if (err) {
      free(dPeriod);
      return err;
    }
    dDer = 10000 * (dFp - dFm) / (2 * dEps);
    pdDerivative[2] = dDer;
  }

  free(dPeriod);
  free(dParams);

  return NULL;
}

/*======================================================================
    Function used to optimized the parameters for the Hull & White model
------------------------------------------------------------------------*/
Err eFutFunc_Bloom(double dIndex, double *pdParams, double *pdConvexity,
                   double *pdDerivative, int nParams) {
  Err err;
  double *dPeriod;
  double *dParams;
  double dFuturesRate;
  double dConv;
  double dEps;
  double dFp, dFm, dDer;
  int i;

  dParams = pdParams;
  dEps = 1.e-4;

  dPeriod = (double *)calloc((long)dIndex + 1, sizeof(double));
  for (i = 0; i <= (long)dIndex; i++)
    dPeriod[i] = dDates[i + 1] - dDates[i];

  err = future_Bloom((long)dIndex, dPeriod, dParams, &dFuturesRate);
  if (err)
    return err;
  dConv = 10000 * (dFuturesRate - dL0T[(long)dIndex - 1]);
  *pdConvexity = dConv;

  dParams[0] += dEps;
  err = future_Bloom((long)dIndex, dPeriod, dParams, &dFp);
  if (err)
    return err;
  dParams[0] -= 2 * dEps;
  err = future_Bloom((long)dIndex, dPeriod, dParams, &dFm);
  if (err)
    return err;
  dDer = 10000 * (dFp - dFm) / (2 * dEps);
  pdDerivative[1] = dDer;

  dParams[0] += dEps;
  dParams[1] += dEps;
  err = future_Bloom((long)dIndex, dPeriod, dParams, &dFp);
  if (err)
    return err;
  dParams[1] -= 2 * dEps;
  err = future_Bloom((long)dIndex, dPeriod, dParams, &dFm);
  if (err)
    return err;
  dDer = 10000 * (dFp - dFm) / (2 * dEps);
  pdDerivative[2] = dDer;

  free(dPeriod);
  free(dParams);

  return NULL;
}
//========================END OF FILE======================================