/* ------------------------------------------------------------------------
   FILENAME:  	srt_f_quanto_fct.c

   PURPOSE:     Function to compute FRA & Swap Adjustment.
   ------------------------------------------------------------------------ */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_ts.h"
#include "srt_h_ts_fx.h"

Err srt_f_quanto_fra(long FraStartDate, long FraEndDate, String FloatRefRate,
                     char *FxUndName, double *AdjFra) {

  Err err = NULL;
  SrtUndPtr sFxUnd, sDomUnd, sForUnd;
  SrtModelType eModelType;
  SrtUnderlyingType eUndType;
  TermStruct *psFxTs, *psDomTs, *psForTs;
  char *szDomUndName, *szForUndName;
  String szfor_crv;
  double G, H;
  double df_start, df_end, cov;
  Date *FixingDates, *StartDates, *EndDates, *PayDates;
  double *Coverage;
  int NumPayDates;
  int VectorSize;
  long today;
  double FraStartMat, FraEndMat;
  double LnFutureAdj, LnFxAdj, FutureAdj, FxAdj;
  SrtCompounding FraFreq;
  BasisCode FraBasis;
  SwapDP *FloatLegDP;
  int SpotLag;

  /*--------------------Get some FX informations
   * -------------------------------*/

  sFxUnd = lookup_und(FxUndName);
  if (sFxUnd == NULL)
    return serror("Undefined Underlying %s", FxUndName);

  eModelType = get_mdltype_from_fxund(sFxUnd);
  if (eModelType != FX_STOCH_RATES)
    return serror("The FX Underlying should be an FX_STOCH_RATES ");

  eUndType = get_underlying_type(sFxUnd);
  if (eUndType != FOREX_UND)
    return serror("Need an Fx Underlying in srt_f_quanto_fra !");

  err = get_underlying_ts(sFxUnd, &psFxTs);
  if (err)
    return serror("Can not get the FX term structure ");

  /*---------------Get some Domestic Ir informations -----------------------*/

  szDomUndName = get_domname_from_fxund(sFxUnd);
  sDomUnd = lookup_und(szDomUndName);
  if (!sDomUnd)
    return serror("Underlying %s not defined", szDomUndName);

  err = get_underlying_ts(sDomUnd, &psDomTs);
  if (err)
    return serror("Can not get the domestic term structure ");

  today = get_today_from_underlying(sDomUnd);

  /*---------------Get some Foreign Ir informations ------------------------*/

  szForUndName = get_forname_from_fxund(sFxUnd);
  sForUnd = lookup_und(szForUndName);
  if (!sForUnd)
    return serror("Underlying %s not defined", szForUndName);

  err = get_underlying_ts(sForUnd, &psForTs);
  if (err)
    return serror("Can not get the foreign term structure");

  szfor_crv = get_discname_from_underlying(sForUnd);

  /*---------------Compute the different elements for the
   * quanto-adjustment---------------*/

  FloatLegDP = (SwapDP *)malloc(sizeof(SwapDP));

  SpotLag = get_spotlag_from_underlying(sDomUnd);

  /* Get the details of the floating reference rate: compounding and basis */
  err = swp_f_get_ref_rate_details(FloatRefRate, &FraBasis, &FraFreq);
  if (err)
    return err;

  /* Initializations of the Floating Leg */

  err =
      swp_f_setSwapDP(FraStartDate, FraEndDate, FraFreq, FraBasis, FloatLegDP);

  if (err)
    return err;

  FloatLegDP->spot_lag = SpotLag;

  err = swp_f_make_FloatLegDatesAndCoverages(
      FloatLegDP, today, &PayDates, &NumPayDates, &FixingDates, &StartDates,
      &EndDates, &Coverage, &VectorSize);
  if (err)
    return err;

  cov = Coverage[VectorSize - 1];
  FraStartMat = (FixingDates[0] - today) * YEARS_IN_DAY;
  FraEndMat = (PayDates[NumPayDates - 1] - today) * YEARS_IN_DAY;

  /* Future Forward Adjustment */

  G_H_func(FraStartMat, psForTs, &G, &H);

  /* This is a Test */

  LnFutureAdj = 0;

  LnFutureAdj += -Lambda_func(FraStartMat, FraEndMat, psForTs) *
                 H_fd_func(FraStartMat, psDomTs, psForTs, psFxTs);

  LnFutureAdj += Lambda_func(FraStartMat, FraEndMat, psForTs) *
                 Lambda_func(FraStartMat, FraEndMat, psDomTs) *
                 Phi_fd_func(FraStartMat, psFxTs);

  LnFutureAdj += -pow(Lambda_func(FraStartMat, FraEndMat, psForTs) *
                          F_func(FraStartMat, psForTs),
                      2) *
                 G;

  LnFutureAdj += -Lambda_func(FraStartMat, FraEndMat, psForTs) *
                 F_func(FraStartMat, psForTs) * Psi_func(FraStartMat, psForTs) *
                 G;

  LnFutureAdj += Lambda_func(FraStartMat, FraEndMat, psForTs) *
                 F_func(FraStartMat, psForTs) * S_func(FraStartMat, psForTs);

  FutureAdj = exp(LnFutureAdj);

  /* Pure FX Adjustment */

  LnFxAdj = 0;

  LnFxAdj += Psi_func(FraStartMat, psForTs) * V_fx_func(FraStartMat, psFxTs);

  LnFxAdj += -W_fx_func(FraStartMat, psFxTs);

  LnFxAdj += -Psi_func(FraEndMat, psForTs) * V_fx_func(FraEndMat, psFxTs);

  LnFxAdj += W_fx_func(FraEndMat, psFxTs);

  FxAdj = exp(LnFxAdj);

  /*--------------Compute the Foreign Forward Rate and adjust it
   * ---------------------*/

  df_start = swp_f_df((Ddate)today, (Ddate)FraStartDate, szfor_crv);
  df_end = swp_f_df((Ddate)today, (Ddate)PayDates[NumPayDates - 1], szfor_crv);

  *AdjFra = ((df_start / df_end) * FutureAdj * FxAdj - 1) / cov;

  if (PayDates)
    srt_free(PayDates);
  if (FixingDates)
    srt_free(FixingDates);
  if (StartDates)
    srt_free(StartDates);
  if (EndDates)
    srt_free(EndDates);
  if (Coverage)
    srt_free(Coverage);
  if (FloatLegDP)
    free(FloatLegDP);

  return err;
}

Err srt_f_quanto_swap(long SwapStartDate, long SwapEndDate,
                      SrtCompounding FixedFreq, SrtBasisCode FixedBasis,
                      String FloatRefRate, char *FxUndName, double *AdjSwap) {

  Err err = NULL;
  SwapDP *FloatLegDP;
  SrtUndPtr sFxUnd, sDomUnd;
  SrtCrvPtr Crv;
  long Today;
  Date *FixingDates, *StartDates, *EndDates, *PayDates;
  double *Coverage;
  long NumPayDates;
  String szdom_crv;
  SrtCcyParam *ccy_param;
  String szDomUndName;
  long WeightsVectorSize, VectorSize;
  double *WeightsVector;
  double AdjFra;
  int i, SpotLag;
  SrtCompounding FloatFreq;
  BasisCode FloatBasis;

  /* double			SwapLevel  ,Weights; */

  /* Get the foreign yield crv */

  sFxUnd = lookup_und(FxUndName);
  if (!sFxUnd)
    return serror("FX Underlying %s not defined", FxUndName);

  szDomUndName = get_domname_from_fxund(sFxUnd);
  sDomUnd = lookup_und(szDomUndName);
  if (!sDomUnd)
    return serror("Underlying %s not defined", szDomUndName);
  szdom_crv = get_discname_from_underlying(sDomUnd);

  /*Get the Fra Start & End Dates */

  FloatLegDP = (SwapDP *)malloc(sizeof(SwapDP));
  Crv = lookup_curve(szdom_crv);
  Today = get_clcndate_from_yldcrv(Crv);

  ccy_param = get_ccyparam_from_yldcrv(Crv);

  SpotLag = ccy_param->spot_lag;

  /* Get the details of the floating reference rate: compounding and basis */
  err = swp_f_get_ref_rate_details(FloatRefRate, &FloatBasis, &FloatFreq);
  if (err)
    return err;

  /* Initializations of the Floating Leg */

  err = swp_f_setSwapDP(SwapStartDate, SwapEndDate, FloatFreq, FloatBasis,
                        FloatLegDP);

  if (err)
    return err;

  FloatLegDP->spot_lag = SpotLag;

  err = swp_f_make_FloatLegDatesAndCoverages(
      FloatLegDP, Today, &PayDates, &NumPayDates, &FixingDates, &StartDates,
      &EndDates, &Coverage, &VectorSize);
  if (err)
    return err;

  /* For the test procedure only at this stage */

  err = srt_SwapFraWeigths(SwapStartDate, SwapEndDate, FixedFreq, FixedBasis,
                           FloatFreq, FloatBasis, szdom_crv, &WeightsVector,
                           &WeightsVectorSize);

  if (err)
    return err;

  (*AdjSwap) = 0;
  for (i = 0; i < (VectorSize - 1); i++) {
    err = srt_f_quanto_fra(StartDates[i], EndDates[i], FloatRefRate, FxUndName,
                           &AdjFra);
    if (err)
      return err;

    (*AdjSwap) += WeightsVector[i] * AdjFra;
  }

  if (WeightsVector)
    srt_free(WeightsVector);
  if (PayDates)
    srt_free(PayDates);
  if (FixingDates)
    srt_free(FixingDates);
  if (StartDates)
    srt_free(StartDates);
  if (EndDates)
    srt_free(EndDates);
  if (Coverage)
    srt_free(Coverage);
  free(FloatLegDP);
  return err;
}

/*------------- COMPUTE THE WEIGHTS OF A SWAP (FOR GRFN
 * USE)---------------------------------*/

Err srt_SwapFraWeigths(Date SwapStartDate, Date SwapEndDate,
                       SrtCompounding FixedFreq, SrtBasisCode FixedBasis,
                       SrtCompounding FloatFreq, SrtBasisCode FloatBasis,
                       String YieldCurveName, double **WeightsVector,
                       long *WeightsVectorSize) {
  Err err = NULL;
  SrtCurvePtr Crv;
  SwapDP *FloatLegDP, *FixedLegDP;

  Date *FixingDates, *StartDates, *EndDates, *PayDates;
  double *Coverage, SwapLevel;
  int NumPayDates, VectorSize;

  Date Today, SpotDate;
  int SpotLag, ScheduleIndex;

  /*Allocations*/

  FixedLegDP = (SwapDP *)malloc(sizeof(SwapDP));
  FloatLegDP = (SwapDP *)malloc(sizeof(SwapDP));

  /* Get Classics informations from the Yield Curve Name */

  Crv = lookup_curve(YieldCurveName);
  Today = get_clcndate_from_yldcrv(Crv);
  SpotDate = get_spotdate_from_yldcrv(Crv);

  SpotLag = get_spotlag_from_curve(Crv);

  if (SwapStartDate < SpotDate)
    SwapStartDate = SpotDate;

  /* Initializations of the Two legs*/
  if (err = swp_f_setSwapDP(SwapStartDate, SwapEndDate, FixedFreq, FixedBasis,
                            FixedLegDP))
    return err;

  FixedLegDP->spot_lag = SpotLag;

  if (err = swp_f_setSwapDP(SwapStartDate, SwapEndDate, FloatFreq, FloatBasis,
                            FloatLegDP))
    return err;

  FloatLegDP->spot_lag = SpotLag;

  /* get the level from the fixed leg of the swap */
  if (err = swp_f_Level_SwapDP(FixedLegDP, YieldCurveName, &SwapLevel))
    return err;

  if (err = swp_f_make_FloatLegDatesAndCoverages(
          FloatLegDP, Today, &PayDates, &NumPayDates, &FixingDates, &StartDates,
          &EndDates, &Coverage, &VectorSize))
    return err;

  (*WeightsVector) = dvector(0, VectorSize - 1);

  for (ScheduleIndex = 0; ScheduleIndex < VectorSize; ScheduleIndex++)

    (*WeightsVector)[ScheduleIndex] =
        swp_f_df(Today, EndDates[ScheduleIndex], YieldCurveName) *
        Coverage[ScheduleIndex] / SwapLevel;

  *WeightsVectorSize = VectorSize;

  /* Free */
  free(FixingDates);
  free(StartDates);
  free(EndDates);
  free(PayDates);
  free(Coverage);

  free(FloatLegDP);
  free(FixedLegDP);

  return NULL;
}

/*------------- COMPUTE THE WEIGHTS OF A SWAP (FOR WESTMINSTER USE
 * ONLY)------------------------------------------*/

Err srt_Weigths(Date SwapStartDate, Date SwapEndDate, String FixedFreqStr,
                String FixedBasisStr, String FloatFreqStr, String FloatBasisStr,
                String YieldCurveName, double **WeightsVector,
                long *WeightsVectorSize) {
  Err err;
  long temp;

  SrtCompounding FloatComp, FixedComp;
  BasisCode FloatBasis, FixedBasis;

  strupper(FloatBasisStr);
  strupper(FixedBasisStr);

  if (err = interp_basis(FixedBasisStr, &FixedBasis))
    return err;

  if (err = interp_basis(FloatBasisStr, &FloatBasis))
    return err;

  if (err = interp_compounding(FloatFreqStr, &FloatComp))
    return err;
  if (err = interp_compounding(FixedFreqStr, &FixedComp))
    return err;

  srt_SwapFraWeigths(SwapStartDate, SwapEndDate, FixedComp, FixedBasis,
                     FloatComp, FloatBasis, YieldCurveName, WeightsVector,
                     &temp);
  *WeightsVectorSize = temp;

  return err;
}

/*COMPUTE THE WEIGHTS OF THE SWAP IN THE DECOMPOSITION OF THE FRA */

Err srt_FraSwapWeights(long SpotDate, long EndDate, SrtCompounding FixedFreq,
                       SrtBasisCode FixedBasis, SrtCompounding FloatFreq,
                       SrtBasisCode FloatBasis, String YieldCurveName,
                       /*OUTPUT*/
                       double ***FraSwapWeightsMatrix, long *FraSwapWeightsSize)

/*FraSwapWeightsSize IS THE NROW OR NLINE OF FraSwapWeightsMatrix */

{

  Err err = NULL;
  SwapDP *FloatLegDP;
  SrtCrvPtr CrvPtr;
  long today;
  long NumPayDates, NumEndDates, WeightsVectorSize;
  double *WeightsVector, *Coverage;
  Date *FixingDates, *StartDates, *EndDates, *PayDates;
  double **SwapFraWeightsMatrix = NULL;
  long i, j;
  long SwapTheoEndDate;

  FloatLegDP = (SwapDP *)malloc(sizeof(SwapDP));
  err = swp_f_setSwapDP(SpotDate, EndDate, FloatFreq, FloatBasis, FloatLegDP);
  if (err)
    return err;

  CrvPtr = lookup_curve(YieldCurveName);
  if (!CrvPtr)
    return serror("Yield Curve Not Found");
  today = get_clcndate_from_yldcrv(CrvPtr);

  err = swp_f_make_FloatLegDatesAndCoverages(
      FloatLegDP, today, &PayDates, &NumPayDates, &FixingDates, &StartDates,
      &EndDates, &Coverage, &NumEndDates);
  if (err)
    return err;

  SwapFraWeightsMatrix = dmatrix(0, NumEndDates - 1, 0, NumEndDates - 1);
  /*NumEndDates IS EQUAL TO THE NUMBER OF SWAP[i]*/
  if (!SwapFraWeightsMatrix)
    return serror("Allocation Memory Failure in srt_FraSwapWeights");

  for (i = 0; i < NumEndDates; i++) {
    SwapTheoEndDate = add_unit(DTOL(SpotDate), (i + 1) * (12 / (int)FloatFreq),
                               SRT_MONTH, NO_BUSDAY_CONVENTION);

    err = srt_SwapFraWeigths(SpotDate, SwapTheoEndDate, FixedFreq, FixedBasis,
                             FloatFreq, FloatBasis, YieldCurveName,
                             &WeightsVector, &WeightsVectorSize);
    /*THE FRA WEIGHTS IN THE DECOMPOSITION OF THE SWAP */
    if (err)
      return err;
    /*SwapFraWeightsMatrix STORES THE FRA WEIGHTS OF THE DECOMPOSITION OF
     * SWAP[i] */
    for (j = 0; j < WeightsVectorSize; j++) {
      SwapFraWeightsMatrix[i][j] = WeightsVector[j];
    }

    for (j = WeightsVectorSize; j < NumEndDates; j++) {
      SwapFraWeightsMatrix[i][j] = 0.0;
    }
  }

  (*FraSwapWeightsMatrix) =
      inverse_matrix(SwapFraWeightsMatrix, 0, NumEndDates - 1);
  if (!(*FraSwapWeightsMatrix))
    return serror("Can Not Inverse SwapFraWeightsMatrix in srt_FraSwapWeights");

  *FraSwapWeightsSize = WeightsVectorSize;

  return err;
}
