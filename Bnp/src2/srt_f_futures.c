/* =======================================================================================

   FILENAME :       srt_f_futuresfct.c

   AUTHOR:          C. Drozo

   DATE:            August-September 2000
                    This version: 29/09/2000

   PURPOSE:         Convexity adjustments between forward and futures contracts
                    Using various models
   =======================================================================================
 */

// Header files & other declarations & definitions

#include "num_h_levenberg.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_futures.h"
#include "srt_h_futuresfct.h"
#include "srtaccess.h"
#include "swp_h_all.h"

double *dDates;
double *dVolStruct;
double *dVolStruct2;
double **ppdFwdVol;
double **ppdOptFutVol;
double *dB0T;
double *dL0T;
double *dFwd;
double *dCoverage;
int iNumFut;
int nbTenors2;

Err srt_f_futures(long lToday, int iType, long *plFutDates, int iNumFut,
                  double *dParams, double *pdFuturesPrice, double *pdConvexity,
                  double *pdRate) {
  double dFuturesRate;
  double dFutPrice, dConv;
  int j, k;
  double *dPeriod;
  Err err = NULL;

  dDates = dvector(0, iNumFut + 1);
  dPeriod = dvector(0, iNumFut);
  ;
  for (j = 0; j <= iNumFut + 1; j++)
    dDates[j] = ((double)plFutDates[j] - (double)lToday) / 365.0;
  for (k = 0; k <= iNumFut; k++)
    dPeriod[k] = dDates[k + 1] - dDates[k];

  if (ppdOptFutVol)
    err = eBuildVols(lToday);
  if (err) {
    free_dvector(dDates, 0, iNumFut + 1);
    dDates = NULL;
    free_dvector(dPeriod, 0, iNumFut);
    dPeriod = NULL;
    if (ppdOptFutVol)
      free_dvector(dVolStruct2, 0, iNumFut - 1);
    dVolStruct2 = NULL;
    return err;
  }

  switch (iType) {
    //=============================================================================
    /* third approx: HJM modelisation of the correlation rs  ,L(T  ,T1) */
  case 0:
    err = future_HJM(iNumFut, dPeriod, dParams, &dFuturesRate);
    if (err) {
      free_dvector(dDates, 0, iNumFut + 1);
      dDates = NULL;
      free_dvector(dPeriod, 0, iNumFut);
      dPeriod = NULL;
      if (ppdOptFutVol)
        free_dvector(dVolStruct2, 0, iNumFut - 1);
      dVolStruct2 = NULL;
      return err;
    }
    break;

    //=============================================================================
    /* fourth approx: HJM-Ho & Lee approx */
  case 1:
    err = future_HL(iNumFut, dPeriod, &dFuturesRate);
    if (err) {
      free_dvector(dDates, 0, iNumFut + 1);
      dDates = NULL;
      free_dvector(dPeriod, 0, iNumFut);
      dPeriod = NULL;
      if (ppdOptFutVol)
        free_dvector(dVolStruct2, 0, iNumFut - 1);
      dVolStruct2 = NULL;
      return err;
    }
    break;

    //=============================================================================
    /* fifth approx: modified Rutkowski                              */
  case 2:
    err = future_Rutkowski(iNumFut, dPeriod, dParams, &dFuturesRate);
    if (err) {
      free_dvector(dDates, 0, iNumFut + 1);
      dDates = NULL;
      free_dvector(dPeriod, 0, iNumFut);
      dPeriod = NULL;
      if (ppdOptFutVol)
        free_dvector(dVolStruct2, 0, iNumFut - 1);
      dVolStruct2 = NULL;
      return err;
    }
    break;

    //=============================================================================
    /* sixth approx: Bloomberg 1                              */
  case 3:
    err = future_Bloom(iNumFut, dPeriod, dParams, &dFuturesRate);
    if (err) {
      free_dvector(dDates, 0, iNumFut + 1);
      dDates = NULL;
      free_dvector(dPeriod, 0, iNumFut);
      dPeriod = NULL;
      if (ppdOptFutVol)
        free_dvector(dVolStruct2, 0, iNumFut - 1);
      dVolStruct2 = NULL;
      return err;
    }
    break;

    //=============================================================================
  default:
    return "Unknown adjustment type";
  }

  /* futures rate and other stuff */

  dFutPrice = 100 * (1.0 - dFuturesRate);
  dConv = 100 * (dFuturesRate - dL0T[iNumFut - 1]);

  *pdFuturesPrice = dFutPrice;
  *pdConvexity = max(0, dConv);
  *pdRate = dFuturesRate;

  free_dvector(dDates, 0, iNumFut + 1);
  dDates = NULL;
  free_dvector(dPeriod, 0, iNumFut);
  dPeriod = NULL;
  return NULL;
}

//================================================================================
Err srt_f_calibfutures(long lToday, long *lFutDates, int iNumFut, int iType,
                       double *alpha, double *beta, double *convexities,
                       double *error) {
  // voir pour lFutDates[0]=Today....
  Err err;
  extern double *dDates;
  double *dIndex;
  double *dWeights;
  double *dParams;
  double *dConv;
  int i, nParams, nIter;
  char buffer[20];

  dWeights = dvector(1, iNumFut);
  for (i = 1; i <= iNumFut; i++)
    dWeights[i] = 1.0; // convexities[i-1];//2*pow(2*i/iNumFut-0.5
                       // ,2);//pow(0.98  ,iNumFut-i);

  dIndex = dvector(1, iNumFut);
  for (i = 1; i <= iNumFut; i++)
    dIndex[i] = (double)i;

  dConv = dvector(1, iNumFut);
  for (i = 1; i <= iNumFut; i++)
    dConv[i] = convexities[i - 1];

  nIter = 10;
  switch (iType) {
  case 0:
    nParams = 2;
    dParams = dvector(1, nParams);
    dParams[1] = *alpha;
    dParams[2] = *beta;
    err = levenberg_marquardt(dIndex, dConv, dWeights, iNumFut, dParams,
                              nParams, nIter, eFutFunc_HJM, error);
    *alpha = dParams[1];
    *beta = dParams[2];
    break;
  case 2:
    nParams = 1;
    dParams = dvector(1, nParams);
    dParams[1] = *beta;
    err = levenberg_marquardt(dIndex, dConv, dWeights, iNumFut, dParams,
                              nParams, nIter, eFutFunc_Rut, error);
    *beta = dParams[1];
    break;
  case 3:
    nParams = 2;
    dParams = dvector(1, nParams);
    dParams[1] = *alpha;
    dParams[2] = *beta;
    err = levenberg_marquardt(dIndex, dConv, dWeights, iNumFut, dParams,
                              nParams, nIter, eFutFunc_Bloom, error);
    *alpha = dParams[1];
    *beta = dParams[2];
    break;
  default:
    sprintf(buffer, "No optimization for this model...\n");
    smessage(buffer);
  }

  free_dvector(dDates, 0, iNumFut + 1);
  free_dvector(dWeights, 1, iNumFut);
  free_dvector(dParams, 1, nParams);
  return NULL;
}

//=============================END OF FILE==================================