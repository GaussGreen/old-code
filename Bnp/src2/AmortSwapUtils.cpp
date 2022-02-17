/* ===================================================================================
   FILENAME:      swp_f_amortswaption.c

   PURPOSE:       Computes price of a european option on an amortized swap
                                in a co-initial Swap market model context
   ===================================================================================
 */

#pragma warning(disable : 4786) // Disable long name warnings

#include "math.h"
#include "num_h_allhdr.h"
#include "swp_h_all.h"
//#include "swp_h_amortswaption.h"
#include "opHeston.h"
#include "opfnctns.h"
#include "swp_h_vol.h"

//------------------------------------------------------------------------
//-----------------Convert an amortizing swap in a bond ------------------
//------------------------------------------------------------------------
Err ConvertAmortSwapWithMarginsInBond(char *cYCname, char *cRefRname,
                                      long StartDate, long EndDate,
                                      SrtCompounding srtFreq,
                                      SrtBasisCode srtBasis, long lNFixNot,
                                      double *dFixNotionals, double *dFixRates,
                                      long lNFloatNot, double *dFloatNotionals,
                                      double *dMargins, long *lNCpns,
                                      Date **lCpnDates, double **dCpns) {
  Err err = NULL;
  int i, k;
  SrtCurvePtr pCurve = lookup_curve(cYCname);
  SwapDP Swap;
  Date lToday;

  long iNFixPayDates, iNFixDates;
  long *lFixPayDates = NULL, *lFixStartDates = NULL, *lFixEndDates = NULL;
  double *dFixCoverages = NULL;

  long iNFloatPayDates, iNFloatDates;
  long *lFloatFixingDates = NULL, *lFloatPayDates = NULL,
       *lFloatStartDates = NULL, *lFloatEndDates = NULL;
  double *dFloatCoverages = NULL;
  double *dFloatSpreads = NULL;

  int FixFloatMult;
  int shortstub, floatshortstub;

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  err = swp_f_setSwapDP(StartDate, EndDate, srtFreq, srtBasis, &Swap);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FixedLegDatesAndCoverages(
      &Swap, lToday, &lFixPayDates, &iNFixPayDates, &lFixStartDates,
      &lFixEndDates, &dFixCoverages, &iNFixDates);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &Swap, lToday, cRefRname, &lFloatPayDates, &iNFloatPayDates,
      &lFloatFixingDates, &lFloatStartDates, &lFloatEndDates, &dFloatCoverages,
      &dFloatSpreads, &iNFloatDates);
  if (err)
    goto FREE_RETURN;

  if ((iNFloatDates != lNFloatNot) || (iNFixDates != lNFixNot)) {
    err = "Amortized Notionals : Wrong Dimensions";
    goto FREE_RETURN;
  }

  shortstub = 0;
  if ((iNFixDates > 1) &&
      ((int)(12 * (lFixEndDates[0] - lFixStartDates[0]) / 365.0 + 0.5) !=
       (int)(12 * (lFixEndDates[1] - lFixStartDates[1]) / 365.0 + 0.5))) {
    shortstub = 1; // stub;
  }

  floatshortstub = 0;
  if (shortstub) {
    for (i = 0; i < iNFloatPayDates; ++i) {
      if (lFloatEndDates[i] == lFixEndDates[0]) {
        floatshortstub = i;
        i = iNFloatPayDates;
      }
    }
  }

  FixFloatMult = (int)((double)(iNFloatDates - floatshortstub - shortstub) /
                           (iNFixDates - shortstub) +
                       0.5);

  *lCpnDates = (Date *)calloc(iNFixPayDates, sizeof(Date));
  *dCpns = (double *)calloc(iNFixPayDates, sizeof(double));

  (*lCpnDates)[0] = lFixPayDates[0];
  (*dCpns)[0] = -dFloatNotionals[0];
  for (i = 1; i <= shortstub; ++i) {
    (*lCpnDates)[i] = lFixPayDates[i];
    (*dCpns)[i] =
        dFixCoverages[i - 1] * dFixRates[i - 1] * dFixNotionals[i - 1];
    for (k = 1; k < floatshortstub + shortstub + 1; ++k) {
      (*dCpns)[i] += -(dFloatNotionals[k] - dFloatNotionals[k - 1] +
                       dFloatNotionals[k - 1] * dFloatCoverages[k - 1] *
                           (dMargins[k - 1] + dFloatSpreads[k - 1])) *
                     (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                      swp_f_df(lToday, lFixPayDates[i], cYCname));
    }
  }

  for (i = shortstub + 1; i < iNFixPayDates - 1; ++i) {
    (*lCpnDates)[i] = lFixPayDates[i];
    (*dCpns)[i] =
        dFixCoverages[i - 1] * dFixRates[i - 1] * dFixNotionals[i - 1];
    for (k = (i - 1 - shortstub) * FixFloatMult + shortstub + floatshortstub +
             1;
         k < (i - shortstub) * FixFloatMult + shortstub + floatshortstub + 1;
         ++k) {
      (*dCpns)[i] += -(dFloatNotionals[k] - dFloatNotionals[k - 1] +
                       dFloatNotionals[k - 1] * dFloatCoverages[k - 1] *
                           (dMargins[k - 1] + dFloatSpreads[k - 1])) *
                     (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                      swp_f_df(lToday, lFixPayDates[i], cYCname));
    }
  }

  (*lCpnDates)[iNFixPayDates - 1] = lFixPayDates[iNFixPayDates - 1];
  (*dCpns)[iNFixPayDates - 1] = dFixCoverages[iNFixPayDates - 2] *
                                dFixRates[iNFixPayDates - 2] *
                                dFixNotionals[iNFixPayDates - 2];
  for (k = (iNFixPayDates - 2 - shortstub) * FixFloatMult + shortstub +
           floatshortstub + 1;
       k < (iNFixPayDates - 1 - shortstub) * FixFloatMult + shortstub +
               floatshortstub;
       ++k) {
    (*dCpns)[iNFixPayDates - 1] +=
        -(dFloatNotionals[k] - dFloatNotionals[k - 1] +
          dFloatNotionals[k - 1] * dFloatCoverages[k - 1] *
              (dMargins[k - 1] + dFloatSpreads[k - 1])) *
        (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
         swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));
  }
  (*dCpns)[iNFixPayDates - 1] +=
      -(-dFloatNotionals[iNFloatPayDates - 2] +
        dFloatNotionals[iNFloatPayDates - 2] *
            dFloatCoverages[iNFloatPayDates - 2] *
            (dMargins[iNFloatPayDates - 2] +
             dFloatSpreads[iNFloatPayDates - 2])) *
      (swp_f_df(lToday, lFloatPayDates[iNFloatPayDates - 1], cYCname) /
       swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));

  *lNCpns = iNFixPayDates;

FREE_RETURN:

  if (lFixPayDates)
    free(lFixPayDates);
  if (lFixStartDates)
    free(lFixStartDates);
  if (lFixEndDates)
    free(lFixEndDates);
  if (dFixCoverages)
    free(dFixCoverages);

  if (lFloatPayDates)
    free(lFloatPayDates);
  if (lFloatFixingDates)
    free(lFloatFixingDates);
  if (lFloatStartDates)
    free(lFloatStartDates);
  if (lFloatEndDates)
    free(lFloatEndDates);
  if (dFloatCoverages)
    free(dFloatCoverages);
  if (dFloatSpreads)
    free(dFloatSpreads);

  return err;
}

Err ConvertAmortSwapWithMarginsInBond2(
    char *cYCname, char *cVCname, char *cRefRname, long StartDate, long EndDate,
    SrtCompounding srtFreq, SrtBasisCode srtBasis, long lNFixNot,
    double *dFixNotionals, double *dFixRates, long lNFloatNot,
    double *dFloatNotionals, double *dMargins,

    long *lNCpns, Date **lCpnDates, double **dCpns,

    int *iNFixPayDates, Date **lFixPayDates, double **dFixCoverages,

    int *iNFloatPayDates, Date **lFloatPayDates, double **dFloatCoverages,
    double **dFloatSpreads,

    int *shortstub, int *floatshortstub) {
  Err err = NULL;
  int i, k;
  SrtCurvePtr pCurve = lookup_curve(cYCname);
  SwapDP Swap;
  Date lToday;

  int iNFloatDates, iNFixDates;
  Date *lFixStartDates = NULL;
  Date *lFixEndDates = NULL;
  Date *lFloatStartDates = NULL;
  Date *lFloatEndDates = NULL;
  Date *lFloatFixingDates = NULL;

  int FixFloatMult;

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  err = swp_f_setSwapDP(StartDate, EndDate, srtFreq, srtBasis, &Swap);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FixedLegDatesAndCoverages(
      &Swap, lToday, lFixPayDates, iNFixPayDates, &lFixStartDates,
      &lFixEndDates, dFixCoverages, &iNFixDates);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &Swap, lToday, cRefRname, lFloatPayDates, iNFloatPayDates,
      &lFloatFixingDates, &lFloatStartDates, &lFloatEndDates, dFloatCoverages,
      dFloatSpreads, &iNFloatDates);
  if (err)
    goto FREE_RETURN;

  if ((iNFloatDates != lNFloatNot) || (iNFixDates != lNFixNot)) {
    err = "Amortized Notionals : Wrong Dimensions";
    goto FREE_RETURN;
  }

  *shortstub = 0;
  if ((iNFixDates > 1) &&
      ((int)(12 * (lFixEndDates[0] - lFixStartDates[0]) / 365.0 + 0.5) !=
       (int)(12 * (lFixEndDates[1] - lFixStartDates[1]) / 365.0 + 0.5))) {
    *shortstub = 1; // stub;
  }

  *floatshortstub = 0;
  if (*shortstub) {
    for (i = 0; i < *iNFloatPayDates; ++i) {
      if (lFloatEndDates[i] == lFixEndDates[0]) {
        *floatshortstub = i;
        i = *iNFloatPayDates;
      }
    }
  }

  FixFloatMult = (int)((double)(iNFloatDates - *floatshortstub - *shortstub) /
                           (iNFixDates - *shortstub) +
                       0.5);

  *lCpnDates = (Date *)calloc(*iNFixPayDates, sizeof(Date));
  *dCpns = (double *)calloc(*iNFixPayDates, sizeof(double));

  (*lCpnDates)[0] = (*lFixPayDates)[0];
  (*dCpns)[0] = -dFloatNotionals[0];
  for (i = 1; i <= *shortstub; ++i) {
    (*lCpnDates)[i] = (*lFixPayDates)[i];
    (*dCpns)[i] =
        (*dFixCoverages)[i - 1] * dFixRates[i - 1] * dFixNotionals[i - 1];
    for (k = 1; k < *floatshortstub + *shortstub + 1; ++k) {
      (*dCpns)[i] += -(dFloatNotionals[k] - dFloatNotionals[k - 1] +
                       dFloatNotionals[k - 1] * (*dFloatCoverages)[k - 1] *
                           (dMargins[k - 1] + (*dFloatSpreads)[k - 1])) *
                     (swp_f_df(lToday, (*lFloatPayDates)[k], cYCname) /
                      swp_f_df(lToday, (*lFixPayDates)[i], cYCname));
    }
  }

  for (i = *shortstub + 1; i < *iNFixPayDates - 1; ++i) {
    (*lCpnDates)[i] = (*lFixPayDates)[i];
    (*dCpns)[i] =
        (*dFixCoverages)[i - 1] * dFixRates[i - 1] * dFixNotionals[i - 1];
    for (k = (i - 1 - *shortstub) * FixFloatMult + *shortstub +
             *floatshortstub + 1;
         k < (i - *shortstub) * FixFloatMult + *shortstub + *floatshortstub + 1;
         ++k) {
      (*dCpns)[i] += -(dFloatNotionals[k] - dFloatNotionals[k - 1] +
                       dFloatNotionals[k - 1] * (*dFloatCoverages)[k - 1] *
                           (dMargins[k - 1] + (*dFloatSpreads)[k - 1])) *
                     (swp_f_df(lToday, (*lFloatPayDates)[k], cYCname) /
                      swp_f_df(lToday, (*lFixPayDates)[i], cYCname));
    }
  }

  (*lCpnDates)[*iNFixPayDates - 1] = (*lFixPayDates)[*iNFixPayDates - 1];
  (*dCpns)[*iNFixPayDates - 1] = (*dFixCoverages)[*iNFixPayDates - 2] *
                                 dFixRates[*iNFixPayDates - 2] *
                                 dFixNotionals[*iNFixPayDates - 2];
  for (k = (*iNFixPayDates - 2 - *shortstub) * FixFloatMult + *shortstub +
           *floatshortstub + 1;
       k < (*iNFixPayDates - 1 - *shortstub) * FixFloatMult + *shortstub +
               *floatshortstub;
       ++k) {
    (*dCpns)[*iNFixPayDates - 1] +=
        -(dFloatNotionals[k] - dFloatNotionals[k - 1] +
          dFloatNotionals[k - 1] * (*dFloatCoverages)[k - 1] *
              (dMargins[k - 1] + (*dFloatSpreads)[k - 1])) *
        (swp_f_df(lToday, (*lFloatPayDates)[k], cYCname) /
         swp_f_df(lToday, (*lFixPayDates)[*iNFixPayDates - 1], cYCname));
  }
  (*dCpns)[*iNFixPayDates - 1] +=
      -(-dFloatNotionals[*iNFloatPayDates - 2] +
        dFloatNotionals[*iNFloatPayDates - 2] *
            (*dFloatCoverages)[*iNFloatPayDates - 2] *
            (dMargins[*iNFloatPayDates - 2] +
             (*dFloatSpreads)[*iNFloatPayDates - 2])) *
      (swp_f_df(lToday, (*lFloatPayDates)[*iNFloatPayDates - 1], cYCname) /
       swp_f_df(lToday, (*lFixPayDates)[*iNFixPayDates - 1], cYCname));

  *lNCpns = *iNFixPayDates;

FREE_RETURN:

  if (lFixStartDates)
    free(lFixStartDates);
  if (lFixEndDates)
    free(lFixEndDates);

  if (lFloatFixingDates)
    free(lFloatFixingDates);
  if (lFloatStartDates)
    free(lFloatStartDates);
  if (lFloatEndDates)
    free(lFloatEndDates);

  return err;
}

//------------------------------------------------------------------------
//-----------------Convert an amortizing swap in a bond For MAD-----------
//------------------------------------------------------------------------
Err ConvertAmortSwapWithMarginsInBondForMAD(
    char *cYCname, char *cVCname, char *cRefRname,

    SrtCompounding srtFreq, SrtBasisCode srtBasis,

    long lNFixDates, long *lFixStartDates, long *lFixEndDates,
    double *dFixCoverages, double *dFixNotionals, double *dFixRates,

    long lNFloatDates, long *lFloatStartDates, long *lFloatEndDates,
    double *dFloatCoverages, double *dFloatNotionals, double *dFloatSpreads,
    double *dMargins,

    long *lNCpns, Date **lCpnDates, double **dCpns,

    int *shortstub, int *floatshortstub) {
  Err err = NULL;
  int i, k;
  SrtCurvePtr pCurve = lookup_curve(cYCname);
  Date lToday;

  int FixFloatMult;

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  *shortstub = 0;
  if ((lNFixDates > 1) &&
      ((int)(12 * (lFixEndDates[0] - lFixStartDates[0]) / 365.0 + 0.5) !=
       (int)(12 * (lFixEndDates[1] - lFixStartDates[1]) / 365.0 + 0.5))) {
    *shortstub = 1; // stub;
  }

  *floatshortstub = 0;
  if (*shortstub) {
    for (i = 0; i < lNFloatDates; ++i) {
      if (lFloatEndDates[i] == lFixEndDates[0]) {
        *floatshortstub = i;
        i = lNFloatDates;
      }
    }
  }

  FixFloatMult = (int)((double)(lNFloatDates - *floatshortstub - *shortstub) /
                           (lNFixDates - *shortstub) +
                       0.5);

  *lCpnDates = (Date *)calloc(lNFixDates + 1, sizeof(Date));
  *dCpns = (double *)calloc(lNFixDates + 1, sizeof(double));

  (*lCpnDates)[0] = lFixStartDates[0];
  (*dCpns)[0] = -dFloatNotionals[0];
  for (i = 1; i <= *shortstub; ++i) {
    (*lCpnDates)[i] = lFixStartDates[i];
    (*dCpns)[i] =
        dFixCoverages[i - 1] * dFixRates[i - 1] * dFixNotionals[i - 1];
    for (k = 1; k < *floatshortstub + *shortstub + 1; ++k) {
      (*dCpns)[i] += -(dFloatNotionals[k] - dFloatNotionals[k - 1] +
                       dFloatNotionals[k - 1] * dFloatCoverages[k - 1] *
                           (dMargins[k - 1] + dFloatSpreads[k - 1])) *
                     (swp_f_df(lToday, lFloatStartDates[k], cYCname) /
                      swp_f_df(lToday, lFixStartDates[i], cYCname));
    }
  }

  for (i = *shortstub + 1; i < lNFixDates; ++i) {
    (*lCpnDates)[i] = lFixStartDates[i];
    (*dCpns)[i] =
        dFixCoverages[i - 1] * dFixRates[i - 1] * dFixNotionals[i - 1];
    for (k = (i - 1 - *shortstub) * FixFloatMult + *shortstub +
             *floatshortstub + 1;
         k < (i - *shortstub) * FixFloatMult + *shortstub + *floatshortstub + 1;
         ++k) {
      (*dCpns)[i] += -(dFloatNotionals[k] - dFloatNotionals[k - 1] +
                       dFloatNotionals[k - 1] * dFloatCoverages[k - 1] *
                           (dMargins[k - 1] + dFloatSpreads[k - 1])) *
                     (swp_f_df(lToday, lFloatStartDates[k], cYCname) /
                      swp_f_df(lToday, lFixStartDates[i], cYCname));
    }
  }

  (*lCpnDates)[lNFixDates] = lFixEndDates[lNFixDates - 1];
  (*dCpns)[lNFixDates] = dFixCoverages[lNFixDates - 1] *
                         dFixRates[lNFixDates - 1] *
                         dFixNotionals[lNFixDates - 1];
  for (k = (lNFixDates - 1 - *shortstub) * FixFloatMult + *shortstub +
           *floatshortstub + 1;
       k <
       (lNFixDates - *shortstub) * FixFloatMult + *shortstub + *floatshortstub;
       ++k) {
    (*dCpns)[lNFixDates] +=
        -(dFloatNotionals[k] - dFloatNotionals[k - 1] +
          dFloatNotionals[k - 1] * dFloatCoverages[k - 1] *
              (dMargins[k - 1] + dFloatSpreads[k - 1])) *
        (swp_f_df(lToday, lFloatEndDates[k], cYCname) /
         swp_f_df(lToday, lFixEndDates[lNFixDates - 1], cYCname));
  }
  (*dCpns)[lNFixDates] +=
      -(-dFloatNotionals[lNFloatDates - 1] +
        dFloatNotionals[lNFloatDates - 1] * dFloatCoverages[lNFloatDates - 1] *
            (dMargins[lNFloatDates - 1] + dFloatSpreads[lNFloatDates - 1])) *
      (swp_f_df(lToday, lFloatEndDates[lNFloatDates - 1], cYCname) /
       swp_f_df(lToday, lFixEndDates[lNFixDates - 1], cYCname));

  *lNCpns = lNFixDates + 1;

  // FREE_RETURN :

  return err;
}

Err ConvertAmortSwapWithMarginsInBondForMAD2(
    char *cYCname, char *cVCname, char *cRefRname,

    SrtCompounding srtFreq, SrtBasisCode srtBasis,

    long lNFixDates, long *lFixStartDates, long *lFixEndDates,
    double *dFixCoverages, double *dFixNotionals, double *dFixRates,

    long lNFloatDates, long *lFloatStartDates, long *lFloatEndDates,
    double *dFloatCoverages, double *dFloatNotionals, double *dFloatSpreads,
    double *dMargins,

    long *lNCpns, Date **lCpnDates, double **dCpns) {
  Err err = NULL;
  int i, k;
  SrtCurvePtr pCurve = lookup_curve(cYCname);
  Date lToday;

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  *lCpnDates = (Date *)calloc(lNFixDates + 1, sizeof(Date));
  *dCpns = (double *)calloc(lNFixDates + 1, sizeof(double));

  (*lCpnDates)[0] = lFixStartDates[0];
  (*dCpns)[0] = -dFloatNotionals[0];

  k = 0;
  for (i = 0; i < lNFixDates - 1; ++i) {
    (*lCpnDates)[i + 1] = lFixEndDates[i];
    (*dCpns)[i + 1] = dFixCoverages[i] * dFixRates[i] * dFixNotionals[i];
    while (lFloatEndDates[k] <= lFixEndDates[i] + 10) {
      (*dCpns)[i + 1] += -(dFloatNotionals[k + 1] - dFloatNotionals[k] +
                           dFloatNotionals[k] * dFloatCoverages[k] *
                               (dMargins[k] + dFloatSpreads[k])) *
                         (swp_f_df(lToday, lFloatEndDates[k], cYCname) /
                          swp_f_df(lToday, lFixEndDates[i], cYCname));
      k = k + 1;
    }
  }

  (*lCpnDates)[i + 1] = lFixEndDates[i];
  (*dCpns)[i + 1] = dFixCoverages[i] * dFixRates[i] * dFixNotionals[i];
  while (lFloatEndDates[k] < lFixEndDates[i] - 10) {
    (*dCpns)[i + 1] += -(dFloatNotionals[k + 1] - dFloatNotionals[k] +
                         dFloatNotionals[k] * dFloatCoverages[k] *
                             (dMargins[k] + dFloatSpreads[k])) *
                       (swp_f_df(lToday, lFloatEndDates[k], cYCname) /
                        swp_f_df(lToday, lFixEndDates[i], cYCname));
    k = k + 1;
  }

  (*dCpns)[i + 1] +=
      -(-dFloatNotionals[k] + dFloatNotionals[k] * dFloatCoverages[k] *
                                  (dMargins[k] + dFloatSpreads[k])) *
      (swp_f_df(lToday, lFloatEndDates[k], cYCname) /
       swp_f_df(lToday, lFixEndDates[i], cYCname));

  *lNCpns = lNFixDates + 1;

  // FREE_RETURN :

  return err;
}

//---------------------------------------------------------------------------
//----------- Convert a vanilla swap in a bond ------------------------------
//---------------------------------------------------------------------------
Err ConvertVanillaSwapInBond(char *cYCname, char *cRefRname, long StartDate,
                             long EndDate, SrtCompounding srtFreq,
                             SrtBasisCode srtBasis, double dFixRate,
                             long *lNCpns, Date **lCpnDates, double **dCpns) {
  Err err = NULL;
  int i, k;
  SrtCurvePtr pCurve = lookup_curve(cYCname);
  SwapDP Swap;
  Date lToday;

  long iNFixPayDates, iNFixDates;
  long *lFixPayDates = NULL, *lFixStartDates = NULL, *lFixEndDates = NULL;
  double *dFixCoverages = NULL;

  long iNFloatPayDates, iNFloatDates;
  long *lFloatFixingDates = NULL, *lFloatPayDates = NULL,
       *lFloatStartDates = NULL, *lFloatEndDates = NULL;
  double *dFloatCoverages = NULL;
  double *dFloatSpreads = NULL;

  int FixFloatMult;

  int shortstub, floatshortstub;

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  err = swp_f_setSwapDP(StartDate, EndDate, srtFreq, srtBasis, &Swap);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FixedLegDatesAndCoverages(
      &Swap, lToday, &lFixPayDates, &iNFixPayDates, &lFixStartDates,
      &lFixEndDates, &dFixCoverages, &iNFixDates);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &Swap, lToday, cRefRname, &lFloatPayDates, &iNFloatPayDates,
      &lFloatFixingDates, &lFloatStartDates, &lFloatEndDates, &dFloatCoverages,
      &dFloatSpreads, &iNFloatDates);
  if (err)
    goto FREE_RETURN;

  shortstub = 0;
  if ((iNFixDates > 1) &&
      ((int)(12 * (lFixEndDates[0] - lFixStartDates[0]) / 365.0 + 0.5) !=
       (int)(12 * (lFixEndDates[1] - lFixStartDates[1]) / 365.0 + 0.5))) {
    shortstub = 1; // stub;
  }

  floatshortstub = 0;
  if (shortstub) {
    for (i = 0; i < iNFloatPayDates; ++i) {
      if (lFloatEndDates[i] == lFixEndDates[0]) {
        floatshortstub = i;
        i = iNFloatPayDates;
      }
    }
  }

  FixFloatMult = (int)((double)(iNFloatDates - floatshortstub - shortstub) /
                           (iNFixDates - shortstub) +
                       0.5);

  *lCpnDates = (Date *)calloc(iNFixPayDates, sizeof(Date));
  *dCpns = (double *)calloc(iNFixPayDates, sizeof(double));

  (*lCpnDates)[0] = lFixPayDates[0];
  (*dCpns)[0] = -1;
  //	FixFloatMult = (int)((double)(iNFloatDates)/iNFixDates+0.5);
  /*
          for(i=1;i<iNFixPayDates-1;++i)
          {
                  (*lCpnDates)[i] = lFixPayDates[i];
                  (*dCpns)[i] = dFixCoverages[i-1]*dFixRate;
                  for(k=(i-1)*FixFloatMult+1;k<i*FixFloatMult+1;++k)
                  {
                          (*dCpns)[i] += - (dFloatCoverages[k-1]
                                                          * dFloatSpreads[k-1])
                                                          * (swp_f_df(lToday  ,
     lFloatPayDates[k]  , cYCname) /swp_f_df(lToday  , lFixPayDates[i]  ,
     cYCname));
                  }
          }

          (*lCpnDates)[iNFixPayDates-1] = lFixPayDates[iNFixPayDates-1];
          (*dCpns)[iNFixPayDates-1] = dFixCoverages[iNFixPayDates-2]*dFixRate;
          for(k=(iNFixPayDates-2)*FixFloatMult+1;k<(iNFixPayDates-1)*FixFloatMult;++k)
          {
                  (*dCpns)[iNFixPayDates-1] += - ( dFloatCoverages[k-1]
                                                                  * dFloatSpreads[k-1])
                                                                  * (swp_f_df(lToday  , lFloatPayDates[k]  , cYCname)
                                                                  /swp_f_df(lToday
     , lFixPayDates[iNFixPayDates-1]  , cYCname));
          }
          (*dCpns)[iNFixPayDates-1] += - ( - 1.0
                                                  +
     dFloatCoverages[iNFloatPayDates-2]
                                                    * dFloatSpreads[iNFloatPayDates-2])
                                          * (swp_f_df(lToday  ,
     lFloatPayDates[iNFloatPayDates-1]  , cYCname) /swp_f_df(lToday  ,
     lFixPayDates[iNFixPayDates-1]  , cYCname));
  */
  for (i = 1; i <= shortstub; ++i) {
    (*lCpnDates)[i] = lFixPayDates[i];
    (*dCpns)[i] = dFixCoverages[i - 1] * dFixRate;
    for (k = 1; k < floatshortstub + shortstub + 1; ++k) {
      (*dCpns)[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                     (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                      swp_f_df(lToday, lFixPayDates[i], cYCname));
    }
  }

  for (i = shortstub + 1; i < iNFixPayDates - 1; ++i) {
    (*lCpnDates)[i] = lFixPayDates[i];
    (*dCpns)[i] = dFixCoverages[i - 1] * dFixRate;
    for (k = (i - 1 - shortstub) * FixFloatMult + shortstub + floatshortstub +
             1;
         k < (i - shortstub) * FixFloatMult + shortstub + floatshortstub + 1;
         ++k) {
      (*dCpns)[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                     (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                      swp_f_df(lToday, lFixPayDates[i], cYCname));
    }
  }

  (*lCpnDates)[iNFixPayDates - 1] = lFixPayDates[iNFixPayDates - 1];
  (*dCpns)[iNFixPayDates - 1] = dFixCoverages[iNFixPayDates - 2] * dFixRate;

  for (k = (iNFixPayDates - 2 - shortstub) * FixFloatMult + shortstub +
           floatshortstub + 1;
       k < (iNFixPayDates - 1 - shortstub) * FixFloatMult + shortstub +
               floatshortstub;
       ++k) {
    (*dCpns)[iNFixPayDates - 1] +=
        -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
        (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
         swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));
  }
  (*dCpns)[iNFixPayDates - 1] +=
      -(-1.0 + dFloatCoverages[iNFloatPayDates - 2] *
                   dFloatSpreads[iNFloatPayDates - 2]) *
      (swp_f_df(lToday, lFloatPayDates[iNFloatPayDates - 1], cYCname) /
       swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));

  *lNCpns = iNFixPayDates;

FREE_RETURN:

  if (lFixPayDates)
    free(lFixPayDates);
  if (lFixStartDates)
    free(lFixStartDates);
  if (lFixEndDates)
    free(lFixEndDates);
  if (dFixCoverages)
    free(dFixCoverages);

  if (lFloatPayDates)
    free(lFloatPayDates);
  if (lFloatFixingDates)
    free(lFloatFixingDates);
  if (lFloatStartDates)
    free(lFloatStartDates);
  if (lFloatEndDates)
    free(lFloatEndDates);
  if (dFloatCoverages)
    free(dFloatCoverages);
  if (dFloatSpreads)
    free(dFloatSpreads);

  return err;
}

Err ConvertVanillaSwapInBond2(char *cYCname,

                              long lToday,

                              int shortstub, int floatshortstub,

                              int iNFixPayDates, long *lFixPayDates,
                              double *dFixCoverages,

                              long iNFloatPayDates, long *lFloatPayDates,
                              double *dFloatCoverages, double *dFloatSpreads,

                              double dFixRate,

                              long *lNCpns, Date **lCpnDates, double **dCpns) {
  Err err = NULL;
  int i, k;
  SrtCurvePtr pCurve = lookup_curve(cYCname);

  int FixFloatMult;

  pCurve = lookup_curve(cYCname);

  FixFloatMult =
      (int)((double)(iNFloatPayDates - 1 - floatshortstub - shortstub) /
                (iNFixPayDates - 1 - shortstub) +
            0.5);

  *lCpnDates = (Date *)calloc(iNFixPayDates, sizeof(Date));
  *dCpns = (double *)calloc(iNFixPayDates, sizeof(double));

  (*lCpnDates)[0] = lFixPayDates[0];
  (*dCpns)[0] = -1;
  for (i = 1; i <= shortstub; ++i) {
    (*lCpnDates)[i] = lFixPayDates[i];
    (*dCpns)[i] = dFixCoverages[i - 1] * dFixRate;
    for (k = 1; k < floatshortstub + shortstub + 1; ++k) {
      (*dCpns)[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                     (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                      swp_f_df(lToday, lFixPayDates[i], cYCname));
    }
  }

  for (i = shortstub + 1; i < iNFixPayDates - 1; ++i) {
    (*lCpnDates)[i] = lFixPayDates[i];
    (*dCpns)[i] = dFixCoverages[i - 1] * dFixRate;
    for (k = (i - 1 - shortstub) * FixFloatMult + shortstub + floatshortstub +
             1;
         k < (i - shortstub) * FixFloatMult + shortstub + floatshortstub + 1;
         ++k) {
      (*dCpns)[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                     (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                      swp_f_df(lToday, lFixPayDates[i], cYCname));
    }
  }

  (*lCpnDates)[iNFixPayDates - 1] = lFixPayDates[iNFixPayDates - 1];
  (*dCpns)[iNFixPayDates - 1] = dFixCoverages[iNFixPayDates - 2] * dFixRate;

  for (k = (iNFixPayDates - 2 - shortstub) * FixFloatMult + shortstub +
           floatshortstub + 1;
       k < (iNFixPayDates - 1 - shortstub) * FixFloatMult + shortstub +
               floatshortstub;
       ++k) {
    (*dCpns)[iNFixPayDates - 1] +=
        -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
        (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
         swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));
  }
  (*dCpns)[iNFixPayDates - 1] +=
      -(-1.0 + dFloatCoverages[iNFloatPayDates - 2] *
                   dFloatSpreads[iNFloatPayDates - 2]) *
      (swp_f_df(lToday, lFloatPayDates[iNFloatPayDates - 1], cYCname) /
       swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));

  *lNCpns = iNFixPayDates;

  return err;
}

Err ConvertVanillaSwapInBond3(char *cYCname,

                              long lToday,

                              int iNFixPayDates, long *lFixPayDates,
                              double *dFixCoverages,

                              long *lFloatPayDates, double *dFloatCoverages,
                              double *dFloatSpreads,

                              double dFixRate,

                              long *lNCpns, Date **lCpnDates, double **dCpns) {
  Err err = NULL;
  int i, k;
  SrtCurvePtr pCurve = lookup_curve(cYCname);

  pCurve = lookup_curve(cYCname);

  *lCpnDates = (Date *)calloc(iNFixPayDates, sizeof(Date));
  *dCpns = (double *)calloc(iNFixPayDates, sizeof(double));

  (*lCpnDates)[0] = lFixPayDates[0];
  (*dCpns)[0] = -1;

  k = 1;
  for (i = 1; i < iNFixPayDates - 1; ++i) {
    (*lCpnDates)[i] = lFixPayDates[i];
    (*dCpns)[i] = dFixCoverages[i - 1] * dFixRate;

    while (lFloatPayDates[k] <= lFixPayDates[i] + 10) {
      (*dCpns)[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                     (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                      swp_f_df(lToday, lFixPayDates[i], cYCname));
      k = k + 1;
    }
  }

  (*lCpnDates)[i] = lFixPayDates[i];
  (*dCpns)[i] = dFixCoverages[i - 1] * dFixRate;

  while (lFloatPayDates[k] < lFixPayDates[i] - 10) {
    (*dCpns)[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                   (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                    swp_f_df(lToday, lFixPayDates[i], cYCname));
    k = k + 1;
  }

  (*dCpns)[i] += -(-1.0 + dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                 (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                  swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));

  *lNCpns = iNFixPayDates;

  return err;
}

//------------------------------------------------------------------------
//-----------------Compute PV of a bond from the IRR ---------------------
//------------------------------------------------------------------------
Err PVBOND_IRR(char *cYCname, char *cVCname, long lToday, char *cFreq,
               char *cBasis, char *cRefRname, int iNPayDates, long *lPayDates,
               double *dCpns, int UseVol, double R, double *Pv) {
  Err err = NULL;
  int i;
  double PV;
  double fwd, vol, power;

  PV = 0.0;
  for (i = 0; i < iNPayDates; ++i) {
    if ((i > 0) && (UseVol == 1)) {
      err = swp_f_ForwardRate(lPayDates[0], lPayDates[i], cFreq, cBasis,
                              cYCname, cRefRname, &fwd);
      if (err) {
        smessage("Error computing the forward swap in PVBondIRR");
        err = "Error computing the forward swap in PVBondIRR";
        goto FREE_RETURN;
      }

      err = swp_f_vol(cVCname, lPayDates[0], lPayDates[i], fwd, &vol, &power);
      if (err) {
        smessage("Error computing the volatility in PVBondIRR");
        err = "Error computing the volatility in PVBondIRR";
        goto FREE_RETURN;
      }
    } else {
      vol = 1;
    }

    PV += dCpns[i] *
          (swp_f_df(lToday, lPayDates[i], cYCname) /
           swp_f_df(lToday, lPayDates[0], cYCname)) *
          exp(-vol * R * (lPayDates[i] - lPayDates[0]) / 365.0);
  }

  *Pv = PV;

FREE_RETURN:

  return err;
}

//------------------------------------------------------------------------
//-----------------Compute the IRR of a bond -----------------------------
//------------------------------------------------------------------------
Err ComputeBondIRR(char *cYCname, char *cVCname, long lToday, char *cFreq,
                   char *cBasis, char *cRefRname, int nPayDates,
                   long *lCpnDates, double *dCpns, double *dIRR, int UseVol) {
  double R, R2;
  double PV, PV2;
  int Compt;
  Err err = NULL;

  R = 0.05;
  //	PV = PVBOND_IRR(cYCname  , lToday  , nPayDates  , lCpnDates  , dCpns  ,
  //R);

  err = PVBOND_IRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                   nPayDates, lCpnDates, dCpns, UseVol, R, &PV);

  Compt = 0;
  while ((fabs(PV) > 1e-12) && (Compt < 15)) {
    R2 = R + 0.0001;
    //		PV2 = PVBOND_IRR(cYCname  , lToday  , nPayDates  , lCpnDates  , dCpns
    //, R2);
    err = PVBOND_IRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                     nPayDates, lCpnDates, dCpns, UseVol, R2, &PV2);

    if (fabs(PV2 - PV) > 1e-14) {
      R = R - PV * (R2 - R) / (PV2 - PV);
      // PV = PVBOND_IRR(cYCname  , lToday  , nPayDates  , lCpnDates  , dCpns  ,
      // R);
      err = PVBOND_IRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                       nPayDates, lCpnDates, dCpns, UseVol, R, &PV);
    } else {
      Compt = 15;
    }

    Compt = Compt + 1;
  }

  *dIRR = R;

  return err;
}

//------------------------------------------------------------------------
//-----------Compute the swap rate corresponding to a given IRR ----------
//------------------------------------------------------------------------
Err ComputeSwapRateFromIRR(char *cYCname, char *cVCname, long lToday,
                           char *cRefRname, long StartDate, long EndDate,
                           SrtCompounding srtFreq, SrtBasisCode srtBasis,
                           double IRR, int UseVol, double *dSwapRate) {
  SwapDP Swap;
  long iNFixPayDates, iNFixDates;
  long *lFixPayDates = NULL, *lFixStartDates = NULL, *lFixEndDates = NULL;
  double *dFixCoverages = NULL;

  long iNFloatPayDates, iNFloatDates;
  long *lFloatFixingDates = NULL, *lFloatPayDates = NULL,
       *lFloatStartDates = NULL, *lFloatEndDates = NULL;
  double *dFloatCoverages = NULL;
  double *dFloatSpreads = NULL;
  int i, k, FixFloatMult;
  double *dCpns_Fix = NULL;
  double *dCpns_Float = NULL;
  double PV_Float, PV_Fix;

  char *cBasis;
  char *cFreq;

  int shortstub, floatshortstub;

  Err err = NULL;

  err = translate_basis(&cBasis, srtBasis);
  err = translate_compounding(&cFreq, srtFreq);

  err = swp_f_setSwapDP(StartDate, EndDate, srtFreq, srtBasis, &Swap);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FixedLegDatesAndCoverages(
      &Swap, lToday, &lFixPayDates, &iNFixPayDates, &lFixStartDates,
      &lFixEndDates, &dFixCoverages, &iNFixDates);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &Swap, lToday, cRefRname, &lFloatPayDates, &iNFloatPayDates,
      &lFloatFixingDates, &lFloatStartDates, &lFloatEndDates, &dFloatCoverages,
      &dFloatSpreads, &iNFloatDates);
  if (err)
    goto FREE_RETURN;

  dCpns_Float = dvector(0, iNFixPayDates - 1);
  dCpns_Fix = dvector(0, iNFixPayDates - 1);

  shortstub = 0;
  if ((iNFixDates > 1) &&
      ((int)(12 * (lFixEndDates[0] - lFixStartDates[0]) / 365.0 + 0.5) !=
       (int)(12 * (lFixEndDates[1] - lFixStartDates[1]) / 365.0 + 0.5))) {
    shortstub = 1; // stub;
  }

  floatshortstub = 0;
  if (shortstub) {
    for (i = 0; i < iNFloatPayDates; ++i) {
      if (lFloatEndDates[i] == lFixEndDates[0]) {
        floatshortstub = i;
        i = iNFloatPayDates;
      }
    }
  }

  FixFloatMult = (int)((double)(iNFloatDates - floatshortstub - shortstub) /
                           (iNFixDates - shortstub) +
                       0.5);

  dCpns_Float[0] = -1.0;
  dCpns_Fix[0] = 0.;

  for (i = 1; i <= shortstub; ++i) {
    dCpns_Fix[i] = dFixCoverages[i - 1];
    for (k = 1; k < shortstub + floatshortstub + 1; ++k) {
      dCpns_Float[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                        (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                         swp_f_df(lToday, lFixPayDates[i], cYCname));
    }
  }

  for (i = shortstub + 1; i < iNFixPayDates - 1; ++i) {
    dCpns_Fix[i] = dFixCoverages[i - 1];
    for (k = (i - 1 - shortstub) * FixFloatMult + shortstub + floatshortstub +
             1;
         k < (i - shortstub) * FixFloatMult + shortstub + floatshortstub + 1;
         ++k) {
      dCpns_Float[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                        (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                         swp_f_df(lToday, lFixPayDates[i], cYCname));
    }
  }

  dCpns_Fix[iNFixPayDates - 1] = dFixCoverages[iNFixPayDates - 2];
  for (k = (iNFixPayDates - 2 - shortstub) * FixFloatMult + shortstub +
           floatshortstub + 1;
       k < (iNFixPayDates - 1 - shortstub) * FixFloatMult + shortstub +
               floatshortstub;
       ++k) {
    dCpns_Float[iNFixPayDates - 1] +=
        -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
        (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
         swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));
  }
  dCpns_Float[iNFixPayDates - 1] +=
      -(-1.0 + dFloatCoverages[iNFloatPayDates - 2] *
                   dFloatSpreads[iNFloatPayDates - 2]) *
      (swp_f_df(lToday, lFloatPayDates[iNFloatPayDates - 1], cYCname) /
       swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));

  /*	for(i=1;i<iNFixPayDates-1;++i)
          {
                  for(k=(i-1)*FixFloatMult+1;k<i*FixFloatMult+1;++k)
                  {
                          dCpns_Float[i] += - (dFloatCoverages[k-1]
                                                          * dFloatSpreads[k-1])
                                                          * (swp_f_df(lToday  ,
     lFloatPayDates[k]  , cYCname) /swp_f_df(lToday  , lFixPayDates[i]  ,
     cYCname));
                  }

                  dCpns_Fix[i] = dFixCoverages[i-1];
          }

          dCpns_Fix[iNFixPayDates-1] = dFixCoverages[iNFixPayDates-2];

          for(k=(iNFixPayDates-2)*FixFloatMult+1;k<(iNFixPayDates-1)*FixFloatMult;++k)
          {
                  dCpns_Float[iNFixPayDates-1] += - ( dFloatCoverages[k-1]
                                                                  * dFloatSpreads[k-1])
                                                                  * (swp_f_df(lToday  , lFloatPayDates[k]  , cYCname)
                                                                  /swp_f_df(lToday
     , lFixPayDates[iNFixPayDates-1]  , cYCname));
          }
          dCpns_Float[iNFixPayDates-1] += - ( - 1.0
                                                  +
     dFloatCoverages[iNFloatPayDates-2]
                                                    * dFloatSpreads[iNFloatPayDates-2])
                                          * (swp_f_df(lToday  ,
     lFloatPayDates[iNFloatPayDates-1]  , cYCname) /swp_f_df(lToday  ,
     lFixPayDates[iNFixPayDates-1]  , cYCname));

  */

  //	PV_Float = PVBOND_IRR(cYCname  , lToday  , iNFixPayDates  , lFixPayDates
  //, dCpns_Float  , IRR);
  err = PVBOND_IRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                   iNFixPayDates, lFixPayDates, dCpns_Float, UseVol, IRR,
                   &PV_Float);

  //	PV_Fix = PVBOND_IRR(cYCname  , lToday  , iNFixPayDates  , lFixPayDates
  //, dCpns_Fix  , IRR);
  err =
      PVBOND_IRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                 iNFixPayDates, lFixPayDates, dCpns_Fix, UseVol, IRR, &PV_Fix);

  //	*dSwapRate = -PVBOND_IRR(cYCname  , lToday  , iNFixPayDates  ,
  //lFixPayDates  , dCpns_Float  , IRR); 	*dSwapRate = *dSwapRate /
  //PVBOND_IRR(cYCname  , lToday  , iNFixPayDates  , lFixPayDates  , dCpns_Fix
  //, IRR);
  *dSwapRate = -PV_Float / PV_Fix;

FREE_RETURN:

  if (dCpns_Float)
    free_dvector(dCpns_Float, 0, iNFixPayDates - 1);
  if (dCpns_Fix)
    free_dvector(dCpns_Fix, 0, iNFixPayDates - 1);

  if (lFixPayDates)
    free(lFixPayDates);
  if (lFixStartDates)
    free(lFixStartDates);
  if (lFixEndDates)
    free(lFixEndDates);
  if (dFixCoverages)
    free(dFixCoverages);

  if (lFloatPayDates)
    free(lFloatPayDates);
  if (lFloatFixingDates)
    free(lFloatFixingDates);
  if (lFloatStartDates)
    free(lFloatStartDates);
  if (lFloatEndDates)
    free(lFloatEndDates);
  if (dFloatCoverages)
    free(dFloatCoverages);
  if (dFloatSpreads)
    free(dFloatSpreads);

  return err;
}

Err ComputeSwapRateFromIRR2(char *cYCname, char *cVCname, long lToday,
                            char *cFreq, char *cBasis, char *cRefRname,
                            int shortstub, int floatshortstub,
                            int iNFixPayDates, long *lFixPayDates,
                            double *dFixCoverages, long iNFloatPayDates,
                            long *lFloatPayDates, double *dFloatCoverages,
                            double *dFloatSpreads, double IRR, int UseVol,
                            double *dSwapRate) {
  int i, k, FixFloatMult;
  double *dCpns_Fix = NULL;
  double *dCpns_Float = NULL;
  double PV_Float, PV_Fix;

  Err err = NULL;

  dCpns_Float = dvector(0, iNFixPayDates - 1);
  dCpns_Fix = dvector(0, iNFixPayDates - 1);

  FixFloatMult =
      (int)((double)(iNFloatPayDates - 1 - floatshortstub - shortstub) /
                (iNFixPayDates - 1 - shortstub) +
            0.5);

  dCpns_Float[0] = -1.0;
  dCpns_Fix[0] = 0.;

  for (i = 1; i <= shortstub; ++i) {
    dCpns_Fix[i] = dFixCoverages[i - 1];
    for (k = 1; k < shortstub + floatshortstub + 1; ++k) {
      dCpns_Float[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                        (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                         swp_f_df(lToday, lFixPayDates[i], cYCname));
    }
  }

  for (i = shortstub + 1; i < iNFixPayDates - 1; ++i) {
    dCpns_Fix[i] = dFixCoverages[i - 1];
    for (k = (i - 1 - shortstub) * FixFloatMult + shortstub + floatshortstub +
             1;
         k < (i - shortstub) * FixFloatMult + shortstub + floatshortstub + 1;
         ++k) {
      dCpns_Float[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                        (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                         swp_f_df(lToday, lFixPayDates[i], cYCname));
    }
  }

  dCpns_Fix[iNFixPayDates - 1] = dFixCoverages[iNFixPayDates - 2];
  for (k = (iNFixPayDates - 2 - shortstub) * FixFloatMult + shortstub +
           floatshortstub + 1;
       k < (iNFixPayDates - 1 - shortstub) * FixFloatMult + shortstub +
               floatshortstub;
       ++k) {
    dCpns_Float[iNFixPayDates - 1] +=
        -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
        (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
         swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));
  }
  dCpns_Float[iNFixPayDates - 1] +=
      -(-1.0 + dFloatCoverages[iNFloatPayDates - 2] *
                   dFloatSpreads[iNFloatPayDates - 2]) *
      (swp_f_df(lToday, lFloatPayDates[iNFloatPayDates - 1], cYCname) /
       swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));

  //	PV_Float = PVBOND_IRR(cYCname  , lToday  , iNFixPayDates  , lFixPayDates
  //, dCpns_Float  , IRR);
  err = PVBOND_IRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                   iNFixPayDates, lFixPayDates, dCpns_Float, UseVol, IRR,
                   &PV_Float);

  //	PV_Fix = PVBOND_IRR(cYCname  , lToday  , iNFixPayDates  , lFixPayDates
  //, dCpns_Fix  , IRR);
  err =
      PVBOND_IRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                 iNFixPayDates, lFixPayDates, dCpns_Fix, UseVol, IRR, &PV_Fix);

  //	*dSwapRate = -PVBOND_IRR(cYCname  , lToday  , iNFixPayDates  ,
  //lFixPayDates  , dCpns_Float  , IRR); 	*dSwapRate = *dSwapRate /
  //PVBOND_IRR(cYCname  , lToday  , iNFixPayDates  , lFixPayDates  , dCpns_Fix
  //, IRR);
  *dSwapRate = -PV_Float / PV_Fix;

  if (dCpns_Float)
    free_dvector(dCpns_Float, 0, iNFixPayDates - 1);
  if (dCpns_Fix)
    free_dvector(dCpns_Fix, 0, iNFixPayDates - 1);

  return err;
}

Err ComputeSwapRateFromIRR3(char *cYCname, char *cVCname, long lToday,
                            char *cFreq, char *cBasis, char *cRefRname,

                            int iNFixPayDates, long *lFixPayDates,
                            double *dFixCoverages,

                            long *lFloatPayDates, double *dFloatCoverages,
                            double *dFloatSpreads, double IRR, int UseVol,
                            double *dSwapRate) {
  int i, k;
  double *dCpns_Fix = NULL;
  double *dCpns_Float = NULL;
  double PV_Float, PV_Fix;

  Err err = NULL;

  dCpns_Float = dvector(0, iNFixPayDates - 1);
  dCpns_Fix = dvector(0, iNFixPayDates - 1);

  dCpns_Float[0] = -1.0;
  dCpns_Fix[0] = 0.;

  k = 1;
  for (i = 1; i < iNFixPayDates - 1; ++i) {
    dCpns_Fix[i] = dFixCoverages[i - 1];
    while (lFloatPayDates[k] <= lFixPayDates[i] + 10) {
      dCpns_Float[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                        (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                         swp_f_df(lToday, lFixPayDates[i], cYCname));
      k = k + 1;
    }
  }

  dCpns_Fix[iNFixPayDates - 1] = dFixCoverages[iNFixPayDates - 2];
  while (lFloatPayDates[k] < lFixPayDates[i] - 10) {
    dCpns_Float[i] += -(dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
                      (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
                       swp_f_df(lToday, lFixPayDates[i], cYCname));
    k = k + 1;
  }

  dCpns_Float[iNFixPayDates - 1] +=
      -(-1.0 + dFloatCoverages[k - 1] * dFloatSpreads[k - 1]) *
      (swp_f_df(lToday, lFloatPayDates[k], cYCname) /
       swp_f_df(lToday, lFixPayDates[iNFixPayDates - 1], cYCname));

  err = PVBOND_IRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                   iNFixPayDates, lFixPayDates, dCpns_Float, UseVol, IRR,
                   &PV_Float);

  err =
      PVBOND_IRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                 iNFixPayDates, lFixPayDates, dCpns_Fix, UseVol, IRR, &PV_Fix);

  *dSwapRate = -PV_Float / PV_Fix;

  if (dCpns_Float)
    free_dvector(dCpns_Float, 0, iNFixPayDates - 1);
  if (dCpns_Fix)
    free_dvector(dCpns_Fix, 0, iNFixPayDates - 1);

  return err;
}

Err ComputeLevel(char *cYCname, long lToday, int iNFixPayDates,
                 long *lFixPayDates, double *dFixCoverages, double *dLvl) {
  int i;
  Err err = NULL;

  *dLvl = 0.0;
  for (i = 0; i < iNFixPayDates - 1; ++i) {
    *dLvl += dFixCoverages[i] * swp_f_df(lToday, lFixPayDates[i + 1], cYCname);
  }

  return err;
}

Err ComputeAmortSwapDiagonalIRRs(char *cYCname, char *cVCname, char *cRefRname,
                                 long StartDate, long EndDate,
                                 SrtCompounding srtFreq, SrtBasisCode srtBasis,
                                 long lNFixNot, double *dFixNotionals,
                                 double *dFixRates, double *dExerciseFees,
                                 long lNFloatNot, double *dFloatNotionals,
                                 double *dMargins, double *dIRRs, int UseVol) {
  Err err = NULL;
  int i;
  long lNCpns;
  Date *lCpnDates = NULL;
  double *dCpns = NULL;

  SrtCurvePtr pCurve = lookup_curve(cYCname);
  SwapDP Swap;
  Date lToday;

  long iNFixPayDates, iNFixDates;
  long *lFixPayDates = NULL, *lFixStartDates = NULL, *lFixEndDates = NULL;
  double *dFixCoverages = NULL;

  long iNFloatPayDates, iNFloatDates;
  long *lFloatFixingDates = NULL, *lFloatPayDates = NULL,
       *lFloatStartDates = NULL, *lFloatEndDates = NULL;
  double *dFloatCoverages = NULL;
  double *dFloatSpreads = NULL;

  int FixFloatMult;

  double dIRR;
  char *cBasis;
  char *cFreq;

  double temp;

  err = translate_basis(&cBasis, srtBasis);
  err = translate_compounding(&cFreq, srtFreq);

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  err = swp_f_setSwapDP(StartDate, EndDate, srtFreq, srtBasis, &Swap);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FixedLegDatesAndCoverages(
      &Swap, lToday, &lFixPayDates, &iNFixPayDates, &lFixStartDates,
      &lFixEndDates, &dFixCoverages, &iNFixDates);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &Swap, lToday, cRefRname, &lFloatPayDates, &iNFloatPayDates,
      &lFloatFixingDates, &lFloatStartDates, &lFloatEndDates, &dFloatCoverages,
      &dFloatSpreads, &iNFloatDates);
  if (err)
    goto FREE_RETURN;

  if ((iNFloatDates != lNFloatNot) || (iNFixDates != lNFixNot)) {
    err = "Amortized Notionals : Wrong Dimensions";
    goto FREE_RETURN;
  }

  FixFloatMult = (int)((double)(iNFloatDates) / iNFixDates + 0.5);

  err = ConvertAmortSwapWithMarginsInBond(
      cYCname, cRefRname, StartDate, EndDate, srtFreq, srtBasis, lNFixNot,
      dFixNotionals, dFixRates, lNFloatNot, dFloatNotionals, dMargins, &lNCpns,
      &lCpnDates, &dCpns);

  for (i = 0; i < iNFixDates; ++i) {
    temp = dCpns[i];
    dCpns[i] = -dFixNotionals[i] - dExerciseFees[i];
    err = ComputeBondIRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                         lNCpns - i, lCpnDates + i, dCpns + i, &dIRR, UseVol);
    if (err) {
      goto FREE_RETURN;
    }
    dIRRs[i] = dIRR;
    dCpns[i] = temp;
  }

FREE_RETURN:

  if (lCpnDates)
    free(lCpnDates);
  if (dCpns)
    free(dCpns);

  if (lFixPayDates)
    free(lFixPayDates);
  if (lFixStartDates)
    free(lFixStartDates);
  if (lFixEndDates)
    free(lFixEndDates);
  if (dFixCoverages)
    free(dFixCoverages);

  if (lFloatPayDates)
    free(lFloatPayDates);
  if (lFloatFixingDates)
    free(lFloatFixingDates);
  if (lFloatStartDates)
    free(lFloatStartDates);
  if (lFloatEndDates)
    free(lFloatEndDates);
  if (dFloatCoverages)
    free(dFloatCoverages);
  if (dFloatSpreads)
    free(dFloatSpreads);

  return err;
}

Err ComputeAmortSwapDiagonalIRRsNew(
    char *cYCname, char *cVCname,

    char *cRefRname, SrtCompounding srtFreq, SrtBasisCode srtBasis,

    long lNFix, long *lFixStartDates, long *lFixEndDates, double *dFixCoverages,
    double *dFixRates, double *dFixNotionals, double *dExerciseFees,

    long lNFloat, long *lFloatStartDates, long *lFloatEndDates,
    double *dFloatCoverages, double *dFloatMargins, double *dFloatSpreads,
    double *dFloatNotionals,

    double *dIRRs, int UseVol) {
  Err err = NULL;
  int i;
  long lNCpns;
  Date *lCpnDates = NULL;
  double *dCpns = NULL;

  SrtCurvePtr pCurve = lookup_curve(cYCname);
  Date lToday;

  int FixFloatMult, shortstub, floatshortstub;

  double dIRR;
  char *cBasis;
  char *cFreq;

  double temp;

  err = translate_basis(&cBasis, srtBasis);
  err = translate_compounding(&cFreq, srtFreq);

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  FixFloatMult =
      (int)(dFixCoverages[lNFix - 1] / dFloatCoverages[lNFloat - 1] + 0.5);

  err = ConvertAmortSwapWithMarginsInBondForMAD(
      cYCname, cVCname, cRefRname,

      srtFreq, srtBasis,

      lNFix, lFixStartDates, lFixEndDates, dFixCoverages, dFixNotionals,
      dFixRates,

      lNFloat, lFloatStartDates, lFloatEndDates, dFloatCoverages,
      dFloatNotionals, dFloatSpreads, dFloatMargins,

      &lNCpns, &lCpnDates, &dCpns,

      &shortstub, &floatshortstub);

  for (i = 0; i < lNFix; ++i) {
    temp = dCpns[i];
    dCpns[i] = -dFixNotionals[i] - dExerciseFees[i];
    err = ComputeBondIRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                         lNCpns - i, lCpnDates + i, dCpns + i, &dIRR, UseVol);
    if (err) {
      goto FREE_RETURN;
    }
    dIRRs[i] = dIRR;
    dCpns[i] = temp;
  }

FREE_RETURN:

  if (lCpnDates)
    free(lCpnDates);
  if (dCpns)
    free(dCpns);

  return err;
}

Err ComputeAmortSwapDiagonalIRRsNew2(
    char *cYCname, char *cVCname,

    char *cRefRname, SrtCompounding srtFreq, SrtBasisCode srtBasis,

    long lNFix, long *lFixStartDates, long *lFixEndDates, double *dFixCoverages,
    double *dFixRates, double *dFixNotionals, double *dExerciseFees,

    long lNFloat, long *lFloatStartDates, long *lFloatEndDates,
    double *dFloatCoverages, double *dFloatMargins, double *dFloatSpreads,
    double *dFloatNotionals,

    double *dIRRs, int UseVol) {
  Err err = NULL;
  int i;
  long lNCpns;
  Date *lCpnDates = NULL;
  double *dCpns = NULL;

  SrtCurvePtr pCurve = lookup_curve(cYCname);
  Date lToday;

  double dIRR;
  char *cBasis;
  char *cFreq;

  double temp;

  err = translate_basis(&cBasis, srtBasis);
  err = translate_compounding(&cFreq, srtFreq);

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  err = ConvertAmortSwapWithMarginsInBondForMAD2(
      cYCname, cVCname, cRefRname,

      srtFreq, srtBasis,

      lNFix, lFixStartDates, lFixEndDates, dFixCoverages, dFixNotionals,
      dFixRates,

      lNFloat, lFloatStartDates, lFloatEndDates, dFloatCoverages,
      dFloatNotionals, dFloatSpreads, dFloatMargins,

      &lNCpns, &lCpnDates, &dCpns);

  for (i = 0; i < lNFix; ++i) {
    temp = dCpns[i];
    dCpns[i] = -dFixNotionals[i] - dExerciseFees[i];
    err = ComputeBondIRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                         lNCpns - i, lCpnDates + i, dCpns + i, &dIRR, UseVol);
    if (err) {
      goto FREE_RETURN;
    }
    dIRRs[i] = dIRR;
    dCpns[i] = temp;
  }

FREE_RETURN:

  if (lCpnDates)
    free(lCpnDates);
  if (dCpns)
    free(dCpns);

  return err;
}

Err ComputeCoInitStrikes(char *cYCname, char *cVCname, char *cRefRname,
                         long StartDate, long EndDate, SrtCompounding srtFreq,
                         SrtBasisCode srtBasis, long lNFixNot,
                         double *dFixNotionals, double *dFixRates,
                         long lNFloatNot, double *dFloatNotionals,
                         double *dMargins, double *dSwapRates, int UseVol) {
  Err err = NULL;
  int i;
  long lNCpns;
  Date *lCpnDates = NULL;
  double *dCpns = NULL;

  SrtCurvePtr pCurve = lookup_curve(cYCname);
  SwapDP Swap;
  Date lToday;

  long iNFixPayDates, iNFixDates;
  long *lFixPayDates = NULL, *lFixStartDates = NULL, *lFixEndDates = NULL;
  double *dFixCoverages = NULL;

  long iNFloatPayDates, iNFloatDates;
  long *lFloatFixingDates = NULL, *lFloatPayDates = NULL,
       *lFloatStartDates = NULL, *lFloatEndDates = NULL;
  double *dFloatCoverages = NULL;
  double *dFloatSpreads = NULL;

  int FixFloatMult;

  char *cFreq, *cBasis;

  double dIRR;
  double dSwapRate;

  err = translate_basis(&cBasis, srtBasis);
  err = translate_compounding(&cFreq, srtFreq);

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  err = swp_f_setSwapDP(StartDate, EndDate, srtFreq, srtBasis, &Swap);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FixedLegDatesAndCoverages(
      &Swap, lToday, &lFixPayDates, &iNFixPayDates, &lFixStartDates,
      &lFixEndDates, &dFixCoverages, &iNFixDates);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &Swap, lToday, cRefRname, &lFloatPayDates, &iNFloatPayDates,
      &lFloatFixingDates, &lFloatStartDates, &lFloatEndDates, &dFloatCoverages,
      &dFloatSpreads, &iNFloatDates);
  if (err)
    goto FREE_RETURN;

  if ((iNFloatDates != lNFloatNot) || (iNFixDates != lNFixNot)) {
    err = "Amortized Notionals : Wrong Dimensions";
    goto FREE_RETURN;
  }

  FixFloatMult = (int)((double)(iNFloatDates) / iNFixDates + 0.5);

  err = ConvertAmortSwapWithMarginsInBond(
      cYCname, cRefRname, StartDate, EndDate, srtFreq, srtBasis, lNFixNot,
      dFixNotionals, dFixRates, lNFloatNot, dFloatNotionals, dMargins, &lNCpns,
      &lCpnDates, &dCpns);

  err = ComputeBondIRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                       lNCpns, lCpnDates, dCpns, &dIRR, UseVol);

  for (i = 0; i < iNFixDates; ++i) {
    ComputeSwapRateFromIRR(cYCname, cVCname, lToday, cRefRname, StartDate,
                           lFixEndDates[i], srtFreq, srtBasis, dIRR, UseVol,
                           &dSwapRate);
    if (err) {
      goto FREE_RETURN;
    }
    dSwapRates[i] = dSwapRate;
  }

FREE_RETURN:

  if (lCpnDates)
    free(lCpnDates);
  if (dCpns)
    free(dCpns);

  if (lFixPayDates)
    free(lFixPayDates);
  if (lFixStartDates)
    free(lFixStartDates);
  if (lFixEndDates)
    free(lFixEndDates);
  if (dFixCoverages)
    free(dFixCoverages);

  if (lFloatPayDates)
    free(lFloatPayDates);
  if (lFloatFixingDates)
    free(lFloatFixingDates);
  if (lFloatStartDates)
    free(lFloatStartDates);
  if (lFloatEndDates)
    free(lFloatEndDates);
  if (dFloatCoverages)
    free(dFloatCoverages);
  if (dFloatSpreads)
    free(dFloatSpreads);

  return err;
}

Err ComputeEquivalentAmounts(int nPayDates, long *lCpnDates, double *dCpns,
                             int NCoInitSwap, int *lNCoupons,
                             long **lCouponDates, double **dCoupons,
                             double *dAmounts) {
  int i, j;
  Err err = NULL;
  double sum = 0;

  if (nPayDates != NCoInitSwap + 1) {
    err = "Dimension Problem in ComputeEquivalentNotionals";
    return err;
  }

  dAmounts[NCoInitSwap - 1] =
      dCpns[NCoInitSwap] / dCoupons[NCoInitSwap - 1][NCoInitSwap];

  for (i = NCoInitSwap - 2; i >= 0; --i) {
    sum = 0;
    for (j = i + 1; j < NCoInitSwap; ++j) {
      sum += dAmounts[j] * dCoupons[j][i + 1];
    }
    dAmounts[i] = (dCpns[i + 1] - sum) / dCoupons[i][i + 1];
  }

  return err;
}

Err ConvertBondInCoInitalSwapsPortfolio(
    char *cYCname, char *cVCname, long StartDate, long EndDate, int nPayDates,
    long *lCpnDates, double *dCpns, char *cRefRname,
    SrtCompounding
        srtFreq, // Bond and coinit swaps should have the same frequency//
    SrtBasisCode srtBasis, double *dStrikes, double *dAmounts, int UseVol) {
  Err err = NULL;
  int i;
  SrtCurvePtr pCurve = lookup_curve(cYCname);
  SwapDP Swap;
  Date lToday;

  double dIRR;

  //	double tempoSwap;

  long iNFixPayDates, iNFixDates;
  long *lFixPayDates = NULL, *lFixStartDates = NULL, *lFixEndDates = NULL;
  double *dFixCoverages = NULL;

  long *lNCoupons;
  long **lCouponDates;
  double **dCoupons;

  long start, end;
  int shortstub, floatshortstub;
  char *cBasis, *cFreq;
  int nMonth;

  long iNFloatPayDates, iNFloatDates;
  long *lFloatFixingDates = NULL, *lFloatPayDates = NULL,
       *lFloatStartDates = NULL, *lFloatEndDates = NULL;
  double *dFloatCoverages = NULL;
  double *dFloatSpreads = NULL;

  err = translate_basis(&cBasis, srtBasis);
  err = translate_compounding(&cFreq, srtFreq);

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  err = swp_f_setSwapDP(StartDate, EndDate, srtFreq, srtBasis, &Swap);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FixedLegDatesAndCoverages(
      &Swap, lToday, &lFixPayDates, &iNFixPayDates, &lFixStartDates,
      &lFixEndDates, &dFixCoverages, &iNFixDates);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &Swap, lToday, cRefRname, &lFloatPayDates, &iNFloatPayDates,
      &lFloatFixingDates, &lFloatStartDates, &lFloatEndDates, &dFloatCoverages,
      &dFloatSpreads, &iNFloatDates);
  if (err)
    goto FREE_RETURN;

  shortstub = 0;
  if ((iNFixDates > 1) &&
      ((int)(12 * (lFixEndDates[0] - lFixStartDates[0]) / 365.0 + 0.5) !=
       (int)(12 * (lFixEndDates[1] - lFixStartDates[1]) / 365.0 + 0.5))) {
    shortstub = 1;
  }

  floatshortstub = 0;
  if (shortstub) {
    for (i = 0; i < iNFloatPayDates; ++i) {
      if (lFloatEndDates[i] == lFixEndDates[0]) {
        floatshortstub = i;
        i = iNFloatPayDates;
      }
    }
  }

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  lNCoupons = (long *)calloc(iNFixDates, sizeof(long));
  lCouponDates = (long **)calloc(iNFixDates, sizeof(long *));
  dCoupons = (double **)calloc(iNFixDates, sizeof(double *));

  err = ComputeBondIRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                       nPayDates, lCpnDates, dCpns, &dIRR, UseVol);

  if (cFreq == "ANNUAL") {
    nMonth = 12;
  } else if (cFreq == "SEMIANNUAL") {
    nMonth = 6;
  } else if (cFreq == "QUARTERLY") {
    nMonth = 3;
  } else if (cFreq == "MONTHLY") {
    nMonth = 1;
  }

  for (i = 0; i < iNFixDates; ++i) {
    start = lFixStartDates[0];

    end = add_unit(lFixStartDates[shortstub], nMonth * (i + 1 - shortstub),
                   SRT_MONTH, NO_BUSDAY_CONVENTION);
    //		end  = lFixEndDates[i];

    ComputeSwapRateFromIRR(cYCname, cVCname, lToday, cRefRname, start, end,
                           srtFreq, srtBasis, dIRR, UseVol, &dStrikes[i]);

    err = ConvertVanillaSwapInBond(cYCname, cRefRname, start, end, srtFreq,
                                   srtBasis, dStrikes[i], &lNCoupons[i],
                                   &lCouponDates[i], &dCoupons[i]);
  }

  err = ComputeEquivalentAmounts(nPayDates, lCpnDates, dCpns, iNFixDates,
                                 lNCoupons, lCouponDates, dCoupons, dAmounts);

FREE_RETURN:

  if (lFixPayDates)
    free(lFixPayDates);
  if (lFixStartDates)
    free(lFixStartDates);
  if (lFixEndDates)
    free(lFixEndDates);
  if (dFixCoverages)
    free(dFixCoverages);

  if (lFloatFixingDates)
    free(lFloatFixingDates);
  if (lFloatPayDates)
    free(lFloatPayDates);
  if (lFloatStartDates)
    free(lFloatStartDates);
  if (lFloatEndDates)
    free(lFloatEndDates);
  if (dFloatCoverages)
    free(dFloatCoverages);
  if (dFloatSpreads)
    free(dFloatSpreads);

  if (lNCoupons)
    free(lNCoupons);
  for (i = 0; i < iNFixDates; ++i) {
    if (dCoupons[i])
      free(dCoupons[i]);
    if (lCouponDates[i])
      free(lCouponDates[i]);
  }
  if (dCoupons)
    free(dCoupons);
  if (lCouponDates)
    free(lCouponDates);

  return err;
}

Err ConvertBondInCoInitalSwapsPortfolio2(
    char *cYCname, char *cVCname,

    char *cRefRname, SrtCompounding srtFreq, SrtBasisCode srtBasis,

    int shortstub, int floatshortstub,

    int nPayDates, long *lCpnDates, double *dCpns,

    int iNFixPayDates, long *lFixPayDates, double *dFixCoverages,

    long iNFloatPayDates, long *lFloatPayDates, double *dFloatCoverages,
    double *dFloatSpreads,

    double *dStrikes, double *dAmounts,

    double *dLvls, double *dFwdSwaps,

    int UseVol) {
  Err err = NULL;
  int i, FixFloatMult;
  SrtCurvePtr pCurve = lookup_curve(cYCname);
  Date lToday;

  char *cFreq, *cBasis;

  double dIRR;

  long *lNCoupons;
  long **lCouponDates;
  double **dCoupons;

  err = translate_basis(&cBasis, srtBasis);
  err = translate_compounding(&cFreq, srtFreq);

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  lNCoupons = (long *)calloc(iNFixPayDates - 1, sizeof(long));
  lCouponDates = (long **)calloc(iNFixPayDates - 1, sizeof(long *));
  dCoupons = (double **)calloc(iNFixPayDates - 1, sizeof(double *));

  if ((!lNCoupons) || (!lCouponDates) || (!dCoupons)) {
    err = "Allocation failed in ConvertBondInCoInitalSwapsPortfolio2";
    goto FREE_RETURN;
  }

  err = ComputeBondIRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                       nPayDates, lCpnDates, dCpns, &dIRR, UseVol);
  if (err) {
    goto FREE_RETURN;
  }

  FixFloatMult =
      (int)((double)(iNFloatPayDates - 1 - floatshortstub - shortstub) /
                (iNFixPayDates - 1 - shortstub) +
            0.5);

  for (i = 0; i < iNFixPayDates - 1; ++i) {
    ComputeSwapRateFromIRR2(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,

                            shortstub, floatshortstub,

                            i + 2, lFixPayDates, dFixCoverages,

                            (i + 1 - shortstub) * FixFloatMult +
                                floatshortstub + shortstub + 1,
                            lFloatPayDates, dFloatCoverages, dFloatSpreads,

                            dIRR, UseVol, &dStrikes[i]);

    ComputeSwapRateFromIRR2(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,

                            shortstub, floatshortstub,

                            i + 2, lFixPayDates, dFixCoverages,

                            (i + 1 - shortstub) * FixFloatMult +
                                floatshortstub + shortstub + 1,
                            lFloatPayDates, dFloatCoverages, dFloatSpreads,

                            0.0, UseVol, &dFwdSwaps[i]);

    ComputeLevel(cYCname, lToday, i + 2, lFixPayDates, dFixCoverages,
                 &dLvls[i]);

    err = ConvertVanillaSwapInBond2(
        cYCname, lToday, shortstub, floatshortstub,

        i + 2, lFixPayDates, dFixCoverages,

        (i + 1 - shortstub) * FixFloatMult + floatshortstub + shortstub + 1,
        lFloatPayDates, dFloatCoverages, dFloatSpreads, dStrikes[i],
        &lNCoupons[i], &lCouponDates[i], &dCoupons[i]);
    if (err) {
      goto FREE_RETURN;
    }
  }

  err = ComputeEquivalentAmounts(nPayDates, lCpnDates, dCpns, iNFixPayDates - 1,
                                 lNCoupons, lCouponDates, dCoupons, dAmounts);
  if (err) {
    goto FREE_RETURN;
  }

FREE_RETURN:

  if (lNCoupons)
    free(lNCoupons);
  for (i = 0; i < iNFixPayDates - 1; ++i) {
    if (dCoupons[i])
      free(dCoupons[i]);
    if (lCouponDates[i])
      free(lCouponDates[i]);
  }
  if (dCoupons)
    free(dCoupons);
  if (lCouponDates)
    free(lCouponDates);

  return err;
}

Err ConvertBondInCoInitalSwapsPortfolio3(
    char *cYCname, char *cVCname,

    char *cRefRname, SrtCompounding srtFreq, SrtBasisCode srtBasis,

    int nPayDates, long *lCpnDates, double *dCpns,

    int iNFixPayDates, long *lFixPayDates, double *dFixCoverages,

    //								 long iNFloatPayDates
    //,
    long *lFloatPayDates, double *dFloatCoverages, double *dFloatSpreads,

    double *dStrikes, double *dAmounts,

    double *dLvls, double *dFwdSwaps,

    int UseVol) {
  Err err = NULL;
  int i;
  SrtCurvePtr pCurve = lookup_curve(cYCname);
  Date lToday;

  char *cFreq, *cBasis;

  double dIRR;

  long *lNCoupons;
  long **lCouponDates;
  double **dCoupons;

  err = translate_basis(&cBasis, srtBasis);
  err = translate_compounding(&cFreq, srtFreq);

  pCurve = lookup_curve(cYCname);
  lToday = (Date)get_today_from_curve(pCurve);

  lNCoupons = (long *)calloc(iNFixPayDates - 1, sizeof(long));
  lCouponDates = (long **)calloc(iNFixPayDates - 1, sizeof(long *));
  dCoupons = (double **)calloc(iNFixPayDates - 1, sizeof(double *));

  if ((!lNCoupons) || (!lCouponDates) || (!dCoupons)) {
    err = "Allocation failed in ConvertBondInCoInitalSwapsPortfolio2";
    goto FREE_RETURN;
  }

  err = ComputeBondIRR(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,
                       nPayDates, lCpnDates, dCpns, &dIRR, UseVol);
  if (err) {
    goto FREE_RETURN;
  }

  for (i = 0; i < iNFixPayDates - 1; ++i) {
    ComputeSwapRateFromIRR3(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,

                            i + 2, lFixPayDates, dFixCoverages,

                            lFloatPayDates, dFloatCoverages, dFloatSpreads,

                            dIRR, UseVol, &dStrikes[i]);

    ComputeSwapRateFromIRR3(cYCname, cVCname, lToday, cFreq, cBasis, cRefRname,

                            i + 2, lFixPayDates, dFixCoverages,

                            lFloatPayDates, dFloatCoverages, dFloatSpreads,

                            0.0, UseVol, &dFwdSwaps[i]);

    ComputeLevel(cYCname, lToday, i + 2, lFixPayDates, dFixCoverages,
                 &dLvls[i]);

    err = ConvertVanillaSwapInBond3(cYCname, lToday,

                                    i + 2, lFixPayDates, dFixCoverages,

                                    lFloatPayDates, dFloatCoverages,
                                    dFloatSpreads, dStrikes[i], &lNCoupons[i],
                                    &lCouponDates[i], &dCoupons[i]);
    if (err) {
      goto FREE_RETURN;
    }
  }

  err = ComputeEquivalentAmounts(nPayDates, lCpnDates, dCpns, iNFixPayDates - 1,
                                 lNCoupons, lCouponDates, dCoupons, dAmounts);
  if (err) {
    goto FREE_RETURN;
  }

FREE_RETURN:

  if (lNCoupons)
    free(lNCoupons);
  for (i = 0; i < iNFixPayDates - 1; ++i) {
    if (dCoupons[i])
      free(dCoupons[i]);
    if (lCouponDates[i])
      free(lCouponDates[i]);
  }
  if (dCoupons)
    free(dCoupons);
  if (lCouponDates)
    free(lCouponDates);

  return err;
}

Err ConvertAmortSwapWithMarginsInCoInitSwapPortfolio(
    char *cYCname, char *cVCname, char *cRefRname, SrtCallPutType srtCallPut,
    long StartDate, long EndDate, SrtCompounding srtFreq, SrtBasisCode srtBasis,
    double exer_fee, long lNFixNot, double *dFixNotionals, double *dFixRates,
    long lNFloatNot, double *dFloatNotionals, double *dMargins,
    double *dStrikes, double *dAmounts, int UseVol) {
  Err err = NULL;
  long lNCpns;
  Date *lCpnDates = NULL;
  double *dCpns = NULL;
  double recpay;

  recpay = -1;
  if (srtCallPut == SRT_CALL)
    recpay = 1;

  err = ConvertAmortSwapWithMarginsInBond(
      cYCname, cRefRname, StartDate, EndDate, srtFreq, srtBasis, lNFixNot,
      dFixNotionals, dFixRates, lNFloatNot, dFloatNotionals, dMargins, &lNCpns,
      &lCpnDates, &dCpns);
  if (err)
    goto FREE_RETURN;

  dCpns[0] = dCpns[0] + recpay * exer_fee;

  err = ConvertBondInCoInitalSwapsPortfolio(
      cYCname, cVCname, StartDate, EndDate, lNCpns, lCpnDates, dCpns, cRefRname,
      srtFreq, srtBasis, dStrikes, dAmounts, UseVol);
  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (lCpnDates)
    free(lCpnDates);
  if (dCpns)
    free(dCpns);

  return err;
}

Err ConvertAmortSwapWithMarginsInCoInitSwapPortfolio2(
    char *cYCname, char *cVCname, char *cRefRname, SrtCallPutType srtCallPut,
    long StartDate, long EndDate, SrtCompounding srtFreq, SrtBasisCode srtBasis,
    double exer_fee, long lNFixNot, double *dFixNotionals, double *dFixRates,
    long lNFloatNot, double *dFloatNotionals, double *dMargins,
    double *dStrikes, double *dAmounts,

    double *dLvls, double *dFwdSwaps,

    int UseVol) {
  Err err = NULL;
  long lNCpns;
  Date *lCpnDates = NULL;
  double *dCpns = NULL;
  double recpay;

  int iNFixPayDates;
  Date *lFixPayDates = NULL;
  double *dFixCoverages = NULL;

  int iNFloatPayDates;
  Date *lFloatPayDates = NULL;
  double *dFloatCoverages = NULL;
  double *dFloatSpreads = NULL;
  int shortstub, floatshortstub;

  recpay = -1;
  if (srtCallPut == SRT_CALL)
    recpay = 1;

  err = ConvertAmortSwapWithMarginsInBond2(
      cYCname, cVCname, cRefRname, StartDate, EndDate, srtFreq, srtBasis,
      lNFixNot, dFixNotionals, dFixRates, lNFloatNot, dFloatNotionals, dMargins,
      &lNCpns, &lCpnDates, &dCpns,

      &iNFixPayDates, &lFixPayDates, &dFixCoverages,

      &iNFloatPayDates, &lFloatPayDates, &dFloatCoverages, &dFloatSpreads,

      &shortstub, &floatshortstub);
  if (err)
    goto FREE_RETURN;

  dCpns[0] = dCpns[0] + recpay * exer_fee;

  err = ConvertBondInCoInitalSwapsPortfolio2(
      cYCname, cVCname,

      cRefRname, srtFreq, srtBasis,

      shortstub, floatshortstub,

      lNCpns, lCpnDates, dCpns,

      iNFixPayDates, lFixPayDates, dFixCoverages,

      iNFloatPayDates, lFloatPayDates, dFloatCoverages, dFloatSpreads,

      dStrikes, dAmounts,

      dLvls, dFwdSwaps,

      UseVol);
  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (lCpnDates)
    free(lCpnDates);
  if (lFixPayDates)
    free(lFixPayDates);
  if (dCpns)
    free(dCpns);
  if (dFixCoverages)
    free(dFixCoverages);
  if (lFloatPayDates)
    free(lFloatPayDates);
  if (dFloatCoverages)
    free(dFloatCoverages);
  if (dFloatSpreads)
    free(dFloatSpreads);

  return err;
}

Err ConvertAmortSwapWithMarginsInCoInitSwapPortfolioForMAD(
    char *cYCname, char *cVCname, char *cRefRname,

    SrtCallPutType srtCallPut, double exer_fee,

    SrtCompounding srtFreq, SrtBasisCode srtBasis,

    long lNFixDates, long *lFixStartDates, long *lFixEndDates,
    double *dFixCoverages, double *dFixNotionals, double *dFixRates,

    long lNFloatDates, long *lFloatStartDates, long *lFloatEndDates,
    double *dFloatCoverages, double *dFloatNotionals, double *dFloatSpreads,
    double *dMargins,

    double *dStrikes, double *dAmounts,

    double *dLvls, double *dFwdSwaps,

    int UseVol) {
  Err err = NULL;
  int i;
  long lNCpns;
  Date *lCpnDates = NULL;
  double *dCpns = NULL;
  double recpay;
  Date *lFixDates = NULL;
  Date *lFloatDates = NULL;

  int shortstub, floatshortstub;

  lFixDates = (Date *)calloc(lNFixDates + 1, sizeof(Date));
  lFloatDates = (Date *)calloc(lNFloatDates + 1, sizeof(Date));

  lFixDates[0] = lFixStartDates[0];
  for (i = 0; i < lNFixDates; ++i) {
    lFixDates[i + 1] = lFixEndDates[i];
  }

  lFloatDates[0] = lFloatStartDates[0];
  for (i = 0; i < lNFloatDates; ++i) {
    lFloatDates[i + 1] = lFloatEndDates[i];
  }

  recpay = -1;
  if (srtCallPut == SRT_CALL)
    recpay = 1;

  err = ConvertAmortSwapWithMarginsInBondForMAD(
      cYCname, cVCname, cRefRname,

      srtFreq, srtBasis,

      lNFixDates, lFixStartDates, lFixEndDates, dFixCoverages, dFixNotionals,
      dFixRates,

      lNFloatDates, lFloatStartDates, lFloatEndDates, dFloatCoverages,
      dFloatNotionals, dFloatSpreads, dMargins,

      &lNCpns, &lCpnDates, &dCpns,

      &shortstub, &floatshortstub);
  if (err)
    goto FREE_RETURN;

  dCpns[0] = dCpns[0] + recpay * exer_fee;

  err = ConvertBondInCoInitalSwapsPortfolio2(
      cYCname, cVCname,

      cRefRname, srtFreq, srtBasis,

      shortstub, floatshortstub,

      lNCpns, lCpnDates, dCpns,

      lNFixDates + 1, lFixDates, dFixCoverages,

      lNFloatDates + 1, lFloatDates, dFloatCoverages, dFloatSpreads,

      dStrikes, dAmounts,

      dLvls, dFwdSwaps,

      UseVol);
  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (lCpnDates)
    free(lCpnDates);
  if (dCpns)
    free(dCpns);

  if (lFixDates)
    free(lFixDates);
  if (lFloatDates)
    free(lFloatDates);

  return err;
}

Err ConvertAmortSwapWithMarginsInCoInitSwapPortfolioForMAD2(
    char *cYCname, char *cVCname, char *cRefRname,

    SrtCallPutType srtCallPut, double exer_fee,

    SrtCompounding srtFreq, SrtBasisCode srtBasis,

    long lNFixDates, long *lFixStartDates, long *lFixEndDates,
    double *dFixCoverages, double *dFixNotionals, double *dFixRates,

    long lNFloatDates, long *lFloatStartDates, long *lFloatEndDates,
    double *dFloatCoverages, double *dFloatNotionals, double *dFloatSpreads,
    double *dMargins,

    double *dStrikes, double *dAmounts,

    double *dLvls, double *dFwdSwaps,

    int UseVol) {
  Err err = NULL;
  int i;
  long lNCpns;
  Date *lCpnDates = NULL;
  double *dCpns = NULL;
  double recpay;
  Date *lFixDates = NULL;
  Date *lFloatDates = NULL;

  lFixDates = (Date *)calloc(lNFixDates + 1, sizeof(Date));
  lFloatDates = (Date *)calloc(lNFloatDates + 1, sizeof(Date));

  lFixDates[0] = lFixStartDates[0];
  for (i = 0; i < lNFixDates; ++i) {
    lFixDates[i + 1] = lFixEndDates[i];
  }

  lFloatDates[0] = lFloatStartDates[0];
  for (i = 0; i < lNFloatDates; ++i) {
    lFloatDates[i + 1] = lFloatEndDates[i];
  }

  recpay = -1;
  if (srtCallPut == SRT_CALL)
    recpay = 1;

  err = ConvertAmortSwapWithMarginsInBondForMAD2(
      cYCname, cVCname, cRefRname,

      srtFreq, srtBasis,

      lNFixDates, lFixStartDates, lFixEndDates, dFixCoverages, dFixNotionals,
      dFixRates,

      lNFloatDates, lFloatStartDates, lFloatEndDates, dFloatCoverages,
      dFloatNotionals, dFloatSpreads, dMargins,

      &lNCpns, &lCpnDates, &dCpns);
  if (err)
    goto FREE_RETURN;

  dCpns[0] = dCpns[0] + recpay * exer_fee;

  err = ConvertBondInCoInitalSwapsPortfolio3(
      cYCname, cVCname,

      cRefRname, srtFreq, srtBasis,

      lNCpns, lCpnDates, dCpns,

      lNFixDates + 1, lFixDates, dFixCoverages,

      lFloatDates, dFloatCoverages, dFloatSpreads,

      dStrikes, dAmounts,

      dLvls, dFwdSwaps,

      UseVol);
  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (lCpnDates)
    free(lCpnDates);
  if (dCpns)
    free(dCpns);

  if (lFixDates)
    free(lFixDates);
  if (lFloatDates)
    free(lFloatDates);

  return err;
}
