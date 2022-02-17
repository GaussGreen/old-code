/* Standard C Libs */
#include "grf_h_all.h"
#include "grf_h_public.h"
#include "srt_h_all.h"
#include "srt_h_calibmidat.h"
#include "srtaccess.h"

/****************************************************************/
/* Calboot with swaptions for Midat  , and so on
 */
/****************************************************************/
Err SrtGrfnSwaptionCalboot(long *plStartDates, long *plEndDates,
                           String szFrequency, String szBasis,
                           String szRefRateCode, double *dStrike,
                           String szRecPay, long lNumDates, String szModelName,
                           /* LGM Parameters */
                           double dTau, SRT_Boolean bUseTwoFactor,
                           double dAlpha, double dGamma, double dRho,
                           /* GRFN Parameters for underlying */
                           String *pszGrfnParamNames,
                           String *pszGrfnParamValues, int iNumGrfnParams,
                           String szYieldCurveName,
                           Err (*pfGetVol)(double dStart, double dEnd,
                                           double dStrike, double dForward,
                                           double dSpread, double *pdBsVol),
                           String szVolType,

                           /* Outputs from this calibration */
                           String szUndName) {
  Err err = NULL;
  Date lToday;
  SrtCurvePtr sCurve;
  SrtGrfnParam sGrfnParams;
  SrtReceiverType *pRecPay;
  SwapDP *sdpArray = NULL;
  SrtUndPtr Underlying;
  char szRealModelName[64];
  char szPremium[32];

  String *pszFrequency, *pszBasis, *pszProductType, *pszRefRate;

  int iNumSigmaRows, iNumSigmaCols, iNumTauRows, iNumTauCols, iNumInstruments,
      *pDerivType, i;

  double **ppdSigmaCurve, **ppdTauCurve, *pdStrike, *pdBondStrike,
      *pdOptionPrice, *pdATMOptionPrice, forward, dForwardSpread, dCurrentVol;

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
  strcat(szUndName, szRealModelName);

  /* Creates an initial underlying with the input values for Tau  , Alpha  ,
   * Gamma and Rho */
  if (err = SrtInitIRUnd(szUndName, szYieldCurveName, szRealModelName,
                         iNumSigmaRows, iNumSigmaCols, ppdSigmaCurve,
                         iNumTauRows, iNumTauCols, ppdTauCurve, 0.0, dAlpha,
                         dGamma, dRho, 0.0, 0.0, 0.0, /* vasicek parms */
                         0, 0, NULL))
    return err;

  Underlying = lookup_und(szUndName);

  /* Free the memory allocated for the Sigma and Tau Cruves */
  free_dmatrix(ppdSigmaCurve, 0, iNumSigmaCols - 1, 0, iNumSigmaRows - 1);
  free_dmatrix(ppdTauCurve, 0, iNumTauCols - 1, 0, iNumTauRows - 1);

  /* Set and Overwrite defaults with user defined parameters */
  if (err = srt_f_set_GrfnParams(iNumGrfnParams, pszGrfnParamNames,
                                 pszGrfnParamValues, &sGrfnParams))
    return err;

  /* Gets the Yield Curve (as well as today) */
  sCurve = lookup_curve(szYieldCurveName);
  lToday = get_today_from_curve(sCurve);

  /* Select all the Fixing Dates on Today or after */
  iNumInstruments = lNumDates;

  /* Allocate memory for the calibration instruments */
  pszFrequency = svector_size(0, iNumInstruments - 1, 64);
  pszBasis = svector_size(0, iNumInstruments - 1, 64);
  pdStrike = (double *)srt_malloc(iNumInstruments * sizeof(double));
  pdBondStrike = (double *)srt_malloc(iNumInstruments * sizeof(double));
  pszRefRate = svector_size(0, iNumInstruments - 1, 64);
  pszProductType = svector_size(0, iNumInstruments - 1, 64);
  pdOptionPrice = (double *)srt_malloc(iNumInstruments * sizeof(double));
  pdATMOptionPrice = (double *)srt_malloc(iNumInstruments * sizeof(double));
  strcpy(szPremium, "PREMIUM");
  pRecPay = (int *)calloc(iNumInstruments, sizeof(int));
  pDerivType = (int *)calloc(iNumInstruments, sizeof(int));

  /* Fills in the Fields for the instruments */
  for (i = 0; i < iNumInstruments; i++) {
    forward = 0.0;
    sprintf(pszFrequency[i], szFrequency);
    sprintf(pszBasis[i], szBasis);
    pdStrike[i] = dStrike[i];
    pdBondStrike[i] = 1.0;
    sprintf(pszRefRate[i], szRefRateCode);

    if (err = interp_rec_pay(szRecPay, &pRecPay[i]))
      return err;

    if (err = interp_struct("SWAPTION", &pDerivType[i]))
      return err;

    if (err = swp_f_ForwardRate(plStartDates[i], plEndDates[i], pszFrequency[i],
                                pszBasis[i], szYieldCurveName, pszRefRate[i],
                                &forward))
      return err;

    dForwardSpread =
        swp_f_spread(plStartDates[i], plEndDates[i], pszRefRate[i]);

    err = pfGetVol(plStartDates[i], plEndDates[i], forward, forward,
                   dForwardSpread, &dCurrentVol);

    err = swp_f_Swaption(plStartDates[i], plEndDates[i], pszFrequency[i],
                         pszBasis[i], dCurrentVol, forward, szRecPay,
                         pszRefRate[i], szYieldCurveName, szPremium, szVolType,
                         &pdATMOptionPrice[i]);

    err = pfGetVol(plStartDates[i], plEndDates[i], pdStrike[i], forward,
                   dForwardSpread, &dCurrentVol);

    err = swp_f_Swaption(plStartDates[i], plEndDates[i], pszFrequency[i],
                         pszBasis[i], dCurrentVol, pdStrike[i], szRecPay,
                         pszRefRate[i], szYieldCurveName, szPremium, szVolType,
                         &pdOptionPrice[i]);
  }

  /* Calls the Calibration Routine */
  sdpArray = (SwapDP *)calloc(iNumInstruments, sizeof(SwapDP));

  for (i = 0; i < iNumInstruments; i++) {
    if (err = swp_f_initSwapDP(plStartDates[i], plEndDates[i], pszFrequency[i],
                               pszBasis[i], &sdpArray[i]))
      return err;

    sdpArray[i].spot_lag = get_spotlag_from_underlying(Underlying);
  }

  /* Call to the main function for calibration by bootstrap */

  err = srt_f_bootstrap_calibrate(Underlying, &sGrfnParams, sdpArray, pdStrike,
                                  pdBondStrike, pDerivType, pRecPay, pszRefRate,
                                  pdOptionPrice, pdATMOptionPrice, dTau, dAlpha,
                                  dGamma, dRho, &iNumInstruments);

  /* Free all */
  /* srt_free(plLongEndDates); */
  free_svector_size(pszFrequency, 0, iNumInstruments - 1, 64);
  free_svector_size(pszBasis, 0, iNumInstruments - 1, 64);
  srt_free(pdStrike);
  srt_free(pdBondStrike);
  /* free(pszRecPay); */
  /* free(pDerivType); */
  free_svector_size(pszRefRate, 0, iNumInstruments - 1, 64);
  free_svector_size(pszProductType, 0, iNumInstruments - 1, 64);
  srt_free(pdOptionPrice);
  srt_free(pdATMOptionPrice);

  /* Return the error message if any */
  return err;
}
