/* ===============================================================================

   FILENAME:       srt_f_grfnflexicap.cxx

   PURPOSE:        Pricing of a FlexiCap with Grfn for a given model
                   Builds a tableau and prices the deal

   ===============================================================================
 */
#include "grf_h_all.h"
#include "grf_h_public.h"
#include "srt_h_all.h"
#include "srt_h_grfnflexicap.h"
#include "srtaccess.h"

/* -------------------------------------------------------------------------- */
/* Cell in the first column of the flexi cap tableau : the FRA (a[0] ,a[1]) */
void GrfnFlexiCapFraCell(String szCashBasis, String szUndName,
                         String szRefRateCode, String szFraCell) {
  sprintf(szFraCell,
          "FRA(a[0        ,i]        ,a[1        ,i]        ,\"%s\"        "
          ",\"%s\"        ,\"%s\")",
          szCashBasis, szUndName, szRefRateCode);
}

void GrfnFlexiCapFloorCell(SrtReceiverType eCapFloor, double dStrike,
                           String szCapFloorCell) {
  if (eCapFloor == SRT_PAYER) {
    sprintf(szCapFloorCell,
            "a[2        ,i]*df(now        ,a[1        ,i])*MAX(0        ,c[0   "
            "     ,i]-%lf)",
            dStrike);
  } else if (eCapFloor == SRT_RECEIVER) {
    sprintf(szCapFloorCell,
            "a[2        ,i]*df(now        ,a[1        ,i])*MAX(0        "
            ",%lf-c[0        ,i])",
            dStrike);
  }
}

void GrfnFlexiDigitalCell(double dEps, String szDigitalCell) {
  sprintf(szDigitalCell,
          "(c[2        ,i]-c[3        ,i])/(a[2        ,i]*df(now        ,a[1  "
          "      ,i])*%lf)",
          2 * dEps);
}

void GrfnFlexiMultOpt(String szMultOptCell) {
  sprintf(szMultOptCell, "c[1        ,i]*c[4        ,i]");
}

static void GrfnNewFlexiCapNumExerciseCell(String szNumExerciseCell) {

  sprintf(szNumExerciseCell, "va:=va + c[4        ,i]");
}

static void GrfnFlexiCapFirstCashFlowCell(long lMaxNumExercise,
                                          String szCashFlowCell) {
  sprintf(szCashFlowCell,
          "if ( ( ( c[2        ,i] > 0 ) && (va<= %d) )         , c[5        "
          ",i]         , 0 ) ",
          lMaxNumExercise);
}

static void GrfnFlexiCapSecondCashFlowCell(long lMaxNumExercise,
                                           String szCashFlowCell) {
  sprintf(szCashFlowCell,
          "if ( ( ( va > %d ) && (va - c[4        ,i] < %d) )         , c[5  "
          "      ,i]*(%d-(va-c[4        ,i]))/c[4        ,i]        , 0 ) ",
          lMaxNumExercise, lMaxNumExercise, lMaxNumExercise);
}

static void NewGrfnFlexiCapCashFlowCell(String szCashFlowCell) {
  sprintf(szCashFlowCell, "c[7        ,i] + c[8        ,i]");
}

/* -------------------------------------------------------------------------- */
/* First cell (top left) of the flexi cap tableau: initialise va */
static void GrfnFlexiCapFirstCell(long lNumAlreadyExercised, String szCashBasis,
                                  String szUndName, String szRefRateCode,
                                  String szFirstCell) {
  char szTemp[GRFN_DEF_ARGBUFSZ];

  GrfnFlexiCapFraCell(szCashBasis, szUndName, szRefRateCode, szTemp);

  sprintf(szFirstCell, "(VA:= %d | %s)", lNumAlreadyExercised, szTemp);
}

/* -------------------------------------------------------------------------- */
/* Cell in the second column of the flexi cap tableau : the discounted floating
 * intrinsic */
static void GrfnFlexiCapDiscFloatCell(SrtReceiverType eCapFloor, double dStrike,
                                      String szDiscFloatCell) {
  if (eCapFloor == SRT_RECEIVER) {
    sprintf(szDiscFloatCell,
            "(%f - c[0        ,i]) * a[2        ,i] * df(now        , a[1      "
            "  ,i])",
            dStrike);
  } else if (eCapFloor == SRT_PAYER) {
    sprintf(szDiscFloatCell,
            "(c[0        ,i] - %f) * a[2        ,i] * df(now        , a[1      "
            "  ,i])",
            dStrike);
  }
}

/* -------------------------------------------------------------------------- */
/* Cell in the third column of the flexi cap tableau : the discounted floating
 * payment */
static void GrfnFlexiCapDiscPaymentCell(String szDiscPaymentCell) {
  sprintf(szDiscPaymentCell, "max ( c[1        ,i]         , 0 ) ");
}

/* -------------------------------------------------------------------------- */
/* Cell in the fourth column of the flexi cap tableau : the number of exercise
 */
static void GrfnFlexiCapNumExerciseCell(String szNumExerciseCell) {
  sprintf(szNumExerciseCell,
          "if ( c[1        ,i] > 0         , va:=va + 1        , 0 )");
}

/* -------------------------------------------------------------------------- */
/* Cell in the fifth and last column of the flexi cap tableau : the cashflow */
static void GrfnFlexiCapCashFlowCell(long lMaxNumExercise,
                                     String szCashFlowCell) {
  sprintf(szCashFlowCell,
          "if ( ( ( c[1        ,i] > 0 ) && (va<= %d) )         , c[2        "
          ",i]         , 0 ) ",
          lMaxNumExercise);
}

/* Function that build the full new grfn tableau for a Flexi cap */
Err NewGrfnFlexiCapMakeTableau(long lNumEventDates, double dStrike, double dEps,
                               long lMaxNumExercise, long lNumAlreadyExercised,
                               SrtReceiverType eCapFloor, String szCashBasis,
                               String szUndName, String szRefRateCode,
                               long *plNumRows, long *plNumCols,
                               GrfnCell ***pppsTableau) {
  Err err = NULL;
  int i;

  /* Sets the tableau dimensions : 10 columns */
  *plNumCols = 10;
  *plNumRows = lNumEventDates;

  /* Allocate space for the Tableau (including strings) */
  *pppsTableau = GrfnCellmatrix(*plNumRows, *plNumCols, GRFN_DEF_ARGBUFSZ);

  /* Fills in the tableau row by row */
  for (i = 0; i < *plNumRows; i++) {
    /* The first column : the FRA */
    if (i > 0) {
      GrfnFlexiCapFraCell(szCashBasis, szUndName, szRefRateCode,
                          (*pppsTableau)[i][0].sval);
    } else
        /* The very first cell is different */
        if (i == 0) {
      GrfnFlexiCapFirstCell(lNumAlreadyExercised, szCashBasis, szUndName,
                            szRefRateCode, (*pppsTableau)[0][0].sval);
    }
    (*pppsTableau)[i][0].type = GRFNSCELL;

    /* The second column: the caplet struck at K*/
    GrfnFlexiCapFloorCell(eCapFloor, dStrike, (*pppsTableau)[i][1].sval);
    (*pppsTableau)[i][1].type = GRFNSCELL;

    /* The third column: the caplet struck at K-eps*/
    GrfnFlexiCapFloorCell(eCapFloor, (dStrike - dEps),
                          (*pppsTableau)[i][2].sval);
    (*pppsTableau)[i][2].type = GRFNSCELL;

    /* The fourth column: the caplet struck at K+eps*/
    GrfnFlexiCapFloorCell(eCapFloor, (dStrike + dEps),
                          (*pppsTableau)[i][3].sval);
    (*pppsTableau)[i][3].type = GRFNSCELL;

    /*The fifth column: the digital */
    GrfnFlexiDigitalCell(dEps, (*pppsTableau)[i][4].sval);
    (*pppsTableau)[i][4].type = GRFNSCELL;

    GrfnFlexiMultOpt((*pppsTableau)[i][5].sval);
    (*pppsTableau)[i][5].type = GRFNSCELL;

    GrfnNewFlexiCapNumExerciseCell((*pppsTableau)[i][6].sval);
    (*pppsTableau)[i][6].type = GRFNSCELL;

    GrfnFlexiCapFirstCashFlowCell(lMaxNumExercise, (*pppsTableau)[i][7].sval);
    (*pppsTableau)[i][7].type = GRFNSCELL;

    GrfnFlexiCapSecondCashFlowCell(lMaxNumExercise, (*pppsTableau)[i][8].sval);
    (*pppsTableau)[i][8].type = GRFNSCELL;

    NewGrfnFlexiCapCashFlowCell((*pppsTableau)[i][9].sval);
    (*pppsTableau)[i][9].type = GRFNSCELL;

  } /* END loop on all tableau rows */

  return err;
}

/* ----------------------------------------------------------------------------------
 */

/* Function that build the full grfn tableau for a Flexi cap */
Err GrfnFlexiCapMakeTableau(long lNumEventDates, double dStrike,
                            long lMaxNumExercise, long lNumAlreadyExercised,
                            SrtReceiverType eCapFloor, String szCashBasis,
                            String szUndName, String szRefRateCode,
                            long *plNumRows, long *plNumCols,
                            GrfnCell ***pppsTableau) {
  Err err = NULL;
  int i;

  /* Sets the tableau dimensions : 5 columns */
  *plNumCols = 5;
  *plNumRows = lNumEventDates;

  /* Allocate space for the Tableau (including strings) */
  *pppsTableau = GrfnCellmatrix(*plNumRows, *plNumCols, GRFN_DEF_ARGBUFSZ);

  /* Fills in the tableau row by row */
  for (i = 0; i < *plNumRows; i++) {
    /* The first column : the FRA */
    if (i > 0) {
      GrfnFlexiCapFraCell(szCashBasis, szUndName, szRefRateCode,
                          (*pppsTableau)[i][0].sval);
    } else
        /* The very first cell is different */
        if (i == 0) {
      GrfnFlexiCapFirstCell(lNumAlreadyExercised, szCashBasis, szUndName,
                            szRefRateCode, (*pppsTableau)[0][0].sval);
    }
    (*pppsTableau)[i][0].type = GRFNSCELL;

    /* The second column: the discounted floating intrinsic */
    GrfnFlexiCapDiscFloatCell(eCapFloor, dStrike, (*pppsTableau)[i][1].sval);
    (*pppsTableau)[i][1].type = GRFNSCELL;

    /* The third column: the discounted floating payment */
    GrfnFlexiCapDiscPaymentCell((*pppsTableau)[i][2].sval);
    (*pppsTableau)[i][2].type = GRFNSCELL;

    /* The fourth column: the number of exercise */
    GrfnFlexiCapNumExerciseCell((*pppsTableau)[i][3].sval);
    (*pppsTableau)[i][3].type = GRFNSCELL;

    /* The fifth column: the cash flow  */
    GrfnFlexiCapCashFlowCell(lMaxNumExercise, (*pppsTableau)[i][4].sval);
    (*pppsTableau)[i][4].type = GRFNSCELL;

  } /* END loop on all tableau rows */

  return err;
}

/* ----------------------------------------------------------------------------
 */
#define FREE_GRFN_FLEXI_CAP_MEMORY                                             \
  {                                                                            \
    if (*pppsGrfnTableau)                                                      \
      grfn_free_GrfnCellmatrix((*pppsGrfnTableau), (*plNumRows),               \
                               (*plNumCols));                                  \
    *pppsGrfnTableau = NULL;                                                   \
    if (*pplAuxRangesLength)                                                   \
      srt_free((*pplAuxRangesLength));                                         \
    *pplAuxRangesLength = NULL;                                                \
    free_dmatrix((*pppdAuxRanges), 0, (*plNumAuxColumns) - 1, 0,               \
                 lNumDates - 1);                                               \
    *pppdAuxRanges = NULL;                                                     \
    free_dmatrix((*pppdLastPath), 0, (*plNumRows) - 1, 0, (*plNumCols) - 1);   \
    *pppdLastPath = NULL;                                                      \
  }

/* ---------------------------------------------------------------------------

                                 MAIN FUNCTION

        Prices a Flexi Cap in Grfn        , building the Tableau internally
   ---------------------------------------------------------------------------
 */

Err SrtGrfnFlexiCapPrice(
    long *plFixingDates, long *plStartDates, long *plPayDates,
    double *pdCoverages, long lNumDates, String szRefRateCode, double dStrike,
    double dEps, long lMaxNumExercise, long lNumAlreadyExercised,
    String szCapFloor, String szUndName, String *pszGrfnParamNames,
    String *pszGrfnParamValues, int iNumGrfnParams, double dFullCapRealPrice,

    /* Outputs from Grfn */
    double *pdFlexiCapPrice, double *pdGrfnFullCapFloorPrice, long *plNumRows,
    long *plNumCols, GrfnCell ***pppsGrfnTableau, double ***pppdLastPath,
    long *plNumAuxColumns, long **pplAuxRangesLength, double ***pppdAuxRanges) {
  Err err = NULL;
  SrtReceiverType eCapFloor;
  long lNumEventdates = 0;
  SrtGrfnParam sGrfnParams;
  SrtUndPtr sUndPtr;
  SrtIOStruct *psIOList;
  SrtBasisCode eBasis;
  SrtCompounding eFrequency;
  String szCashBasis;
  double *pdColumnPvs;
  long lNumPvs;
  double dGrfnFlexiCapPrice;
  double dGrfnFullCapPrice;
  int i;
  long lToday;

  /* Gets the Underlying from its name */
  sUndPtr = lookup_und(szUndName);

  /* Gets today from the Underlying */
  lToday = get_today_from_underlying(sUndPtr);

  /* Set and Overwrite defaults with user defined parameters */
  if (err = srt_f_set_GrfnParams(iNumGrfnParams, pszGrfnParamNames,
                                 pszGrfnParamValues, &sGrfnParams)) {
    return err;
  }

  /* Sets the Receiver/Payer type */
  err = interp_rec_pay(szCapFloor, &eCapFloor);
  if (err)
    return err;

  /* Gets the Reference Rate information : basis        , frequency */
  err = swp_f_get_ref_rate_details(szRefRateCode, &eBasis, &eFrequency);
  if (err)
    return err;

  /* Transforms the basis back into a string */
  translate_basis(&szCashBasis, eBasis);

  /* Create a I/O list for the prices */
  err = srt_f_IOstructcreate(&psIOList, "");
  if (err)
    return (err);

  /* Remove everything that is before Today */
  while (plFixingDates[0] < lToday + sGrfnParams.end_of_day_flg) {
    lNumDates--;
    plFixingDates++;
    plStartDates++;
    plPayDates++;
    pdCoverages++;
  }
  /* Construct a Grfn Tableau that describes the Flexi Cap  */
  (*pppsGrfnTableau) = NULL;

  if (dEps > 0) {
    err = NewGrfnFlexiCapMakeTableau(lNumDates, dStrike, dEps, lMaxNumExercise,
                                     lNumAlreadyExercised, eCapFloor,
                                     szCashBasis, szUndName, szRefRateCode,
                                     plNumRows, plNumCols, pppsGrfnTableau);
  }

  else {
    err = GrfnFlexiCapMakeTableau(lNumDates, dStrike, lMaxNumExercise,
                                  lNumAlreadyExercised, eCapFloor, szCashBasis,
                                  szUndName, szRefRateCode, plNumRows,
                                  plNumCols, pppsGrfnTableau);
  }

  if (err) {
    srt_f_IOstructfree(&psIOList);
    return err;
  }

  /* Sets the Auxiliary ranges ; A[0] = start        , A[1] = end        , A[2]
   * = coverage
   */
  (*plNumAuxColumns) = 3;
  (*pplAuxRangesLength) = (long *)srt_malloc((*plNumAuxColumns) * sizeof(long));
  (*pplAuxRangesLength)[0] = lNumDates;
  (*pplAuxRangesLength)[1] = lNumDates;
  (*pplAuxRangesLength)[2] = lNumDates;
  (*pppdAuxRanges) = dmatrix(0, (*plNumAuxColumns) - 1, 0, lNumDates - 1);
  for (i = 0; i < lNumDates; i++) {
    (*pppdAuxRanges)[0][i] = (double)plStartDates[i];
    (*pppdAuxRanges)[1][i] = (double)plPayDates[i];
    (*pppdAuxRanges)[2][i] = (double)pdCoverages[i];
  }

  /* Allocate space for the Grfn Cells (results of the last path) */
  (*pppdLastPath) = dmatrix(0, (*plNumRows) - 1, 0, (*plNumCols) - 1);

  /* Call GRFN to value the deal */
  err = srt_f_grfn(sUndPtr, &sGrfnParams, lNumDates, &plFixingDates, plNumRows,
                   plNumCols, pppsGrfnTableau, 0, NULL, (*plNumAuxColumns),
                   (*pplAuxRangesLength), (*pppdAuxRanges), psIOList,
                   (*pppdLastPath), 0);
  if (err) {
    FREE_GRFN_FLEXI_CAP_MEMORY;
    return err;
  }

  /* Gets all the columns PV from the I/O list */
  err = srt_f_IOstructgetcolpvs(*psIOList, &pdColumnPvs, &lNumPvs);
  if (err) {
    FREE_GRFN_FLEXI_CAP_MEMORY;
    return err;
  }

  /* The Flexi Cap PV is in the last column */
  dGrfnFlexiCapPrice = pdColumnPvs[(*plNumCols) - 1];

  /* The Full Cap PV is in the third column */
  if (dEps > 0)
    dGrfnFullCapPrice = pdColumnPvs[1];
  else
    dGrfnFullCapPrice = pdColumnPvs[2];

  *pdGrfnFullCapFloorPrice = dGrfnFullCapPrice;
  /* Use some sort of control variate if it can be done with the full cap */
  *pdFlexiCapPrice = dGrfnFlexiCapPrice;
  if (dFullCapRealPrice > 0) {
    *pdFlexiCapPrice =
        dGrfnFlexiCapPrice + (double)lMaxNumExercise / (double)lNumDates *
                                 (dFullCapRealPrice - dGrfnFullCapPrice);
  } else {
    *pdFlexiCapPrice = dGrfnFlexiCapPrice;
  }

  /* Returns the success message */
  return err;

} /* END Err SrtGrfnFlexiCapPrice(...) */

Err SrtGrfnFlexiCapCalibrate(
    long *plFixingDates, long *plStartDates, long lNumDates,
    String szRefRateCode, double dStrike, String szCapFloor, String szModelName,
    double dTau, SRT_Boolean bUseTwoFactor, double dAlpha, double dGamma,
    double dRho, String *pszGrfnParamNames, String *pszGrfnParamValues,
    int iNumGrfnParams, String szYieldCurveName,
    Err (*pfGetVol)(double dStart, double dEnd, double dStrike, double dForward,
                    double dSpread, double *pdBsVol),
    String szVolType,
    /* Outputs from this calibration */
    String szUndName, double **ppdCapletRealPrices, long *plNumCaplets) {
  Err err = NULL;
  int iNumSigmaRows;
  int iNumSigmaCols;
  double **ppdSigmaCurve;
  int iNumTauRows;
  int iNumTauCols;
  double **ppdTauCurve;
  char szRealModelName[64];
  SrtBasisCode eBasis;
  SrtCompounding eFrequency;
  String szBasis;
  String szFrequency;
  Date lToday;
  SrtCurvePtr sCurve;
  int iNumInstruments;
  long *plNfp;
  String *pszFrequency;
  String *pszBasis;
  double *pdStrike;
  double *pdBondStrike;
  String *pszProductType;
  String *pszRecPay;
  String *pszRefRate;
  double *pdOptionPrice;
  double *pdATMOptionPrice;
  double forward;
  SrtGrfnParam sGrfnParams;
  int i;
  char szPremium[32];

  /* Builds the Full Model Name */
  if (szModelName == NULL) {
    sprintf(szRealModelName, "LGM");
  } else {
    sprintf(szRealModelName, szModelName);
  }

  if (bUseTwoFactor == SRT_YES) {
    strcat(szRealModelName, "2F");
  }

  /* Creates a fake initial Term Structure with just one value for Sigma and Tau
   */
  iNumSigmaRows = 1;
  iNumSigmaCols = 1;
  ppdSigmaCurve = dmatrix(0, iNumSigmaCols - 1, 0, iNumSigmaRows - 1);
  ppdSigmaCurve[0][0] = 0.02;
  iNumTauRows = 1;
  iNumTauCols = 1;
  ppdTauCurve = dmatrix(0, iNumTauCols - 1, 0, iNumTauRows - 1);
  ppdTauCurve[0][0] = dTau;

  /* Sets a default name for the underlying */
  /* sprintf(szUndName        ,"FlexiCap%sUnd"        ,szModelName); */
  strcat(szUndName, szRealModelName);

  /* Creates an initial underlying with the input values for Tau        , Alpha
   * , Gamma and Rho */
  err = SrtInitIRUnd(szUndName, szYieldCurveName, szRealModelName,
                     iNumSigmaRows, iNumSigmaCols, ppdSigmaCurve, iNumTauRows,
                     iNumTauCols, ppdTauCurve, 0.0, dAlpha, dGamma, dRho, 0.0,
                     0.0, 0.0, /* vasicek parms */
                     0, 0, NULL);

  /* Free the memory allocated for the Sigma and Tau Cruves */
  free_dmatrix(ppdSigmaCurve, 0, iNumSigmaCols - 1, 0, iNumSigmaRows - 1);
  free_dmatrix(ppdTauCurve, 0, iNumTauCols - 1, 0, iNumTauRows - 1);

  if (err)
    return err;

  /* Set and Overwrite defaults with user defined parameters */
  if (err = srt_f_set_GrfnParams(iNumGrfnParams, pszGrfnParamNames,
                                 pszGrfnParamValues, &sGrfnParams)) {
    return err;
  }

  /* Gets the Reference Rate information : basis        , frequency */
  err = swp_f_get_ref_rate_details(szRefRateCode, &eBasis, &eFrequency);
  if (err)
    return err;

  /* Transforms the basis and frequency back into a string */
  translate_basis(&szBasis, eBasis);
  translate_compounding(&szFrequency, eFrequency);

  /* Gets the Yield Curve (as well as today) */
  sCurve = lookup_curve(szYieldCurveName);
  lToday = get_today_from_curve(sCurve);

  /* Select all the Fixing Dates on Today or after */
  iNumInstruments = lNumDates;
  while (plFixingDates[0] < lToday + sGrfnParams.end_of_day_flg) {
    if (iNumInstruments <= 1) {
      return serror("All caplets have expired in SrtGrfnFlexiCapCalibrate");
    } else {
      iNumInstruments--;
      plFixingDates++;
      plStartDates++;
    }
  }

  /* Allocate memory for the calibration instruments (caplets in this case) */
  plNfp = (long *)srt_malloc(iNumInstruments * sizeof(long));
  pszFrequency = svector_size(0, iNumInstruments - 1, 64);
  pszBasis = svector_size(0, iNumInstruments - 1, 64);
  pdStrike = (double *)srt_malloc(iNumInstruments * sizeof(double));
  pdBondStrike = (double *)srt_malloc(iNumInstruments * sizeof(double));
  pszRecPay = svector_size(0, iNumInstruments - 1, 64);
  pszRefRate = svector_size(0, iNumInstruments - 1, 64);
  pszProductType = svector_size(0, iNumInstruments - 1, 64);
  pdOptionPrice = (double *)srt_malloc(iNumInstruments * sizeof(double));
  pdATMOptionPrice = (double *)srt_malloc(iNumInstruments * sizeof(double));
  strcpy(szPremium, "PREMIUM");

  /* Fills in the Fields for the instruments */
  for (i = 0; i < iNumInstruments; i++) {
    forward = 0.0;
    plNfp[i] = 1;
    sprintf(pszFrequency[i], szFrequency);
    sprintf(pszBasis[i], szBasis);
    pdStrike[i] = dStrike;
    pdBondStrike[i] = 1.0;
    sprintf(pszRecPay[i], szCapFloor);
    sprintf(pszRefRate[i], szRefRateCode);
    sprintf(pszProductType[i], "CAPFLOOR");
    err = swp_f_CapFloor(plStartDates[i], plNfp[i], pdStrike[i], pfGetVol,
                         pszRecPay[i], pszRefRate[i], szYieldCurveName,
                         szPremium, szVolType, &pdOptionPrice[i]);
    err = swp_f_ForwardRate(plStartDates[i], plNfp[i], pszFrequency[i],
                            pszBasis[i], szYieldCurveName, pszRefRate[i],
                            &forward);
    err = swp_f_CapFloor(plStartDates[i], plNfp[i], forward, pfGetVol,
                         pszRecPay[i], pszRefRate[i], szYieldCurveName,
                         szPremium, szVolType, &pdATMOptionPrice[i]);
  }

  /* Calls the Calibration Routine */
  err = SrtBootstrap(iNumGrfnParams, pszGrfnParamNames, pszGrfnParamValues,
                     &iNumInstruments, plStartDates, plNfp, pszFrequency,
                     pszBasis, pdStrike, pdBondStrike, pszProductType,
                     pszRecPay, pszRefRate, pdOptionPrice, pdATMOptionPrice,
                     dTau, dAlpha, dGamma, dRho, szUndName);

  /* Transfers the Caplets Prices information */
  *plNumCaplets = iNumInstruments;
  (*ppdCapletRealPrices) = dvector(0, *plNumCaplets - 1);
  for (i = 0; i < *plNumCaplets; i++)
    (*ppdCapletRealPrices)[i] = pdOptionPrice[i];

  /* Free all */
  srt_free(plNfp);
  free_svector_size(pszFrequency, 0, iNumInstruments - 1, 64);
  free_svector_size(pszBasis, 0, iNumInstruments - 1, 64);
  srt_free(pdStrike);
  srt_free(pdBondStrike);
  free_svector_size(pszRecPay, 0, iNumInstruments - 1, 64);
  free_svector_size(pszRefRate, 0, iNumInstruments - 1, 64);
  free_svector_size(pszProductType, 0, iNumInstruments - 1, 64);
  srt_free(pdOptionPrice);
  srt_free(pdATMOptionPrice);

  /* Return the error message if any */
  return err;

} /* END Err SrtGrfnFlexiCapCalibrate(...) */
