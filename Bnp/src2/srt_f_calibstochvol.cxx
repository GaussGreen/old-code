/* ===============================================================================

   FILENAME:  srt_f_calibstochvol

   PURPOSE	: Performs calibration for stochastic volatility
              LGM_STOCH_VOL        , CHEY_STOCH_VOL and CHEY_BETA_STOCH_VOL
              using fixed point algorithm

  Last Modified by Ezra Nahum July 1999
   ===============================================================================
*/

#include "grf_h_public.h"
#include "srt_h_all.h"
#include "srt_h_calib.h"
#include "srt_h_closedform.h"
#include "srt_h_grfclsdfrm.h"
#include "swp_h_all.h"
#include "utallhdr.h"
#ifdef PVMPI
#include "parallel.h"
#endif
#include "math.h"

static struct {
  GrfnCell **pshtSHEET;
  Date *pszEventDates;
  long lNumEventDates;
  long lNumRows;
  long lNumCols;
} STATIC_SPRDSHT;

#define FREE_STOCHVOLCALIB_MEMORY                                              \
  {                                                                            \
    if (pdMarketPrices)                                                        \
      free_dvector(pdMarketPrices, 1, (long)iNumInstruments);                  \
    if (beta)                                                                  \
      free_dvector(beta, 1, (long)iNumInstruments);                            \
    if (pdPrice_hat)                                                           \
      free_dvector(pdPrice_hat, 1, (long)iNumInstruments);                     \
    if (lambda)                                                                \
      free_dvector(lambda, 1, (long)iNumInstruments);                          \
    if (previousbeta)                                                          \
      free_dvector(previousbeta, 1, (long)iNumInstruments);                    \
    if (previouspdPrice_hat)                                                   \
      free_dvector(previouspdPrice_hat, 1, (long)iNumInstruments);             \
    if (pdCalibratedPrices)                                                    \
      free_dvector(pdCalibratedPrices, 1, (long)iNumInstruments);              \
    if (determin_model_resid_error)                                            \
      free_dvector(determin_model_resid_error, 1, (long)iNumInstruments);      \
    if (stoch_model_resid_error)                                               \
      free_dvector(stoch_model_resid_error, 1, (long)iNumInstruments);         \
    if (stoch_error)                                                           \
      free_dvector(stoch_error, 1, 20);                                        \
    if (impl_vol_mkt)                                                          \
      free_dvector(impl_vol_mkt, 1, (long)iNumInstruments);                    \
    if (impl_vol_det)                                                          \
      free_dvector(impl_vol_det, 1, (long)iNumInstruments);                    \
    if (impl_vol_stc)                                                          \
      free_dvector(impl_vol_stc, 1, (long)iNumInstruments);                    \
    pdMarketPrices = NULL;                                                     \
    beta = NULL;                                                               \
    lambda = NULL;                                                             \
    pdPrice_hat = NULL;                                                        \
    previousbeta = NULL;                                                       \
    previouspdPrice_hat = NULL;                                                \
    pdCalibratedPrices = NULL;                                                 \
    determin_model_resid_error = NULL;                                         \
    stoch_error = NULL;                                                        \
    impl_vol_mkt = NULL;                                                       \
    impl_vol_det = NULL;                                                       \
    impl_vol_stc = NULL;                                                       \
    stoch_model_resid_error = NULL;                                            \
    free_dvector(stoch_error, 1, 20);                                          \
    if (STATIC_SPRDSHT.pszEventDates)                                          \
      srt_free(STATIC_SPRDSHT.pszEventDates);                                  \
    if (STATIC_SPRDSHT.pshtSHEET)                                              \
      grfn_free_GrfnCellmatrix(STATIC_SPRDSHT.pshtSHEET,                       \
                               STATIC_SPRDSHT.lNumRows,                        \
                               STATIC_SPRDSHT.lNumCols);                       \
  }

/* global pointer to the array of computed calibrated model prices
   modified in levenberg_calib_funcs */
extern double *GlobalTheoPrices;

Err srt_f_calib_stoch_vol(
    SrtGrfnParam *psGrfnParams, SrtMdlType eModelType, SrtMdlDim eModelDim,
    SwapDP *psSwapDp, double *pdStrike, double *pdBondStrike,
    StructType *peOptionType, SrtReceiverType *peRecPay, double *pdPrice,
    double *pdVega, String *pszRefRateCode, int iNumInstruments,
    double *dFraMaturities, double **ppdCorrelationMatrix, long lNumTenors,
    SrtUndPtr sUndPtr, SrtCalibParam *psCalibParams, double **ppdSigmaCurve,
    long lNumSigmas, long lNumSigmaCols, double **ppdTauCurve, long lNumTaus,
    long lNumTauCols, double *pdChiSquare) {
  double *pdMarketPrices = NULL; /* market pdPrices */
  double *beta = NULL;           /* stochastic volatility model pdPrices */
  double *previousbeta =
      NULL; /* previous value of beta        , used to compute the lambdas */
  double *pdPrice_hat = NULL; /* adjusted "market" pdPrices */
  double *lambda =
      NULL; /* "optimal" coefficient in the fixed point algorithm */
  double *determin_model_resid_error =
      NULL; /* residual error deterministic volatility model */
  double *pdCalibratedPrices =
      NULL; /* pdPrices from deterministic volatility model calibrated
         on market pdPrices - step "m=0" */
  double *stoch_model_resid_error =
      NULL;               /* residual error stochastic volatility model */
  double deter_sum_error; /* "XI" for the deterministic volatility model - step
                             "m==0" */
  double stoch_sum_error; /* "XI" for the stochastic volatility model; */
  double *stoch_error = NULL;  /* XI stoch. vol across iterations*/
  double *impl_vol_mkt = NULL; /* BS implied volatility - market */
  double *impl_vol_det =
      NULL; /* BS implied volatility - deterministic volatility model*/
  double *impl_vol_stc =
      NULL; /* BS implied volatility - stochastic volatility model*/
  double *previouspdPrice_hat = NULL; /*BS implied volatility - target prices */
  SrtMdlType eDeterModelType;
  double rho;
  double vovol;
  double *pdTempPrice = NULL;
  SrtIOStruct *iolist;
  String und_name;
  Date today;
  /* Cheyette beta */
  double ch_beta;

  int m;
  int max_iter;
  int i;
  Date value_date;

  TermStruct *ts;

  char *obj_name = NULL;
  char *yc_name = NULL;

  Err err = NULL;

  SRT_Boolean initial_setting_force_mc;

  /* Allocate memory */
  if (!(pdMarketPrices = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }
  if (!(beta = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  if (!(pdPrice_hat = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  if (!(pdCalibratedPrices = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  if (!(determin_model_resid_error = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  if (!(stoch_model_resid_error = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  if (!(impl_vol_mkt = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  if (!(impl_vol_det = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  if (!(impl_vol_stc = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  if (!(previouspdPrice_hat = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  if (!(stoch_error = dvector(1, 20))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  if (!(previousbeta = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  if (!(lambda = dvector(1, (long)iNumInstruments))) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("memory allocation error");
  }

  /* initializes the rho and the vovol in the SigmaCurve*/

  m = 0;
  /* max_iter limits the time of the calibration stopping it after a while*/
  /*even if the criterion has not been satisfied*/
  max_iter = 8;
  i = 0;

  deter_sum_error = 0;
  stoch_sum_error = 0;

  value_date = get_today_from_underlying(sUndPtr);

  err = get_underlying_mdltype(sUndPtr, &eModelType);
  if (err) {
    FREE_STOCHVOLCALIB_MEMORY;
    return err;
  }

  while (i < iNumInstruments) {
    pdMarketPrices[i] = pdPrice[i];
    pdPrice_hat[i] = pdPrice[i];
    lambda[i] = 1;
    i++;
  }

  /* Get the initial vovol and rho structure */
  get_underlying_ts(sUndPtr, &ts);
  vovol = find_vovol(1, ts);
  rho = find_rho(1, ts);
  ch_beta = find_beta(1, ts);

  /* Get the yield curve from the underlying */
  yc_name = get_ycname_from_underlying(sUndPtr);

  /* Compute the BS implied volatilities for the market pdPrices */
  for (i = 0; i < iNumInstruments; i++) {
    err = swp_f_SwaptionImpliedVol_SwapDP(
        pdMarketPrices[i], &psSwapDp[i], pdStrike[i], peRecPay[i],
        pszRefRateCode[i], yc_name, SRT_LOGNORMAL, &impl_vol_mkt[i]);
    if (err) {
      FREE_STOCHVOLCALIB_MEMORY;
      return serror("error in implied vol");
    }
  }

  if (eModelType == LGM_STOCH_VOL)
    eDeterModelType = LGM;
  if (eModelType == CHEY_STOCH_VOL)
    eDeterModelType = CHEY;
  if (eModelType == CHEY_BETA_STOCH_VOL && vovol != 0)
    eDeterModelType = CHEY_BETA;

  /* The following condition recognizes that the model is a CHEY_BETA */
  /* if the vovol is equal to zero*/

  if (eModelType == CHEY_BETA_STOCH_VOL && vovol == 0) {
    eModelType = CHEY_BETA;
    eDeterModelType = CHEY_BETA;
  }

  /* Gets the initial (input) settings for the  discretisation */
  initial_setting_force_mc = psGrfnParams->force_mc;
  /* Creating the Grfn spreadsheet that will be used to price the deals*/

  und_name = get_underlying_name(sUndPtr);
  today = get_today_from_underlying(sUndPtr);

  err = new_grfn_SwapDParray_to_GrfnCells(
      &STATIC_SPRDSHT.lNumEventDates, &STATIC_SPRDSHT.pszEventDates,
      &STATIC_SPRDSHT.lNumRows, &STATIC_SPRDSHT.lNumCols,
      &STATIC_SPRDSHT.pshtSHEET, today, iNumInstruments, psSwapDp, pdStrike,
      pdBondStrike, peRecPay, peOptionType, und_name, pszRefRateCode);

  if (err) {
    FREE_STOCHVOLCALIB_MEMORY;
    return serror("error in static spreadsheet definition");
  }

  /*--------------------------------------------------------------------------------------------------*/
  /* Now starts the calibration and in particular        , the fixed point
   * algorithm*/
  /*--------------------------------------------------------------------------------------------------*/

  while (m <= max_iter) {
    smessage("Stochastic volatility model calibration. Start iteration: %d", m);

    /* Update model_type in the underlying */
    set_irund_mdltype(sUndPtr, eDeterModelType);

    /* Residual error from calibration of deterministic volatility model
    to be used as convergence criterion for stochastic volatility model.
    pdPrice the input caps and swaptions using the calibrated
    deterministic volatility model */

    /* For Calibration        , use the user defined Grfn Parameters for
     * Discretisation */
    psGrfnParams->force_mc = initial_setting_force_mc;

    /* We make sure that calibmain knows that we are doing a FIXED POINT algo*/

    psCalibParams->eAlgoType = FIXED_POINT;

    /* Perfom calibration of the deterministic volatility model */

    if (err = srt_f_calib_main(
            psGrfnParams, eDeterModelType, eModelDim, psSwapDp, pdStrike,
            pdBondStrike, peOptionType, peRecPay, pdPrice_hat, pdVega,
            pszRefRateCode, iNumInstruments, dFraMaturities,
            ppdCorrelationMatrix, lNumTenors, sUndPtr, psCalibParams,
            ppdSigmaCurve, lNumSigmas, lNumSigmaCols, ppdTauCurve, lNumTaus,
            lNumTauCols, pdChiSquare)) {
      FREE_STOCHVOLCALIB_MEMORY;
      return err;
    }

    /* Update model_type in the underlying */
    set_irund_mdltype(sUndPtr, eModelType);

    /* Insert RHO and VOVOL in the Calibrated TermStructure */
    err = srt_f_init_IRM_TermStruct(value_date, ppdSigmaCurve, 2, lNumSigmas,
                                    ppdTauCurve, 2, lNumTaus, eModelType,
                                    eModelDim, ch_beta, 0.0, 0.0, rho, vovol,
                                    /* BETAETA - Cheyette beta*/
                                    0.0, 0.0, /* vasicek parms */
                                    0, 0, NULL, &ts);

    /* Update ts in mktptr */
    set_irund_ts(sUndPtr, ts);

    /* pdPrice with stochastic volatility Cheyette/LGM using Monte Carlo */
    psGrfnParams->force_mc = SRT_YES;

    err = srt_f_IOstructcreate(&iolist, "aggregate");
    /* Call Grfn with the tableau from the Spreadsheet*/
    if (!err) {
      err = srt_f_grfn(sUndPtr, psGrfnParams, STATIC_SPRDSHT.lNumEventDates,
                       &STATIC_SPRDSHT.pszEventDates, &STATIC_SPRDSHT.lNumRows,
                       &STATIC_SPRDSHT.lNumCols, &STATIC_SPRDSHT.pshtSHEET, 0,
                       0, 0, 0, 0, iolist, 0, 0);
    } else {
      FREE_STOCHVOLCALIB_MEMORY;
      return serror("memory allocation error");
    }

    /* Computes the prices for the given set of parameters */
    pdTempPrice = dvector(0, iNumInstruments - 1);

    err = srt_f_IOstructgetcolpvs((*iolist), &pdTempPrice,
                                  &STATIC_SPRDSHT.lNumCols);
    if (err) {
      err = srt_f_IOstructfree(&iolist);
      /*free_underlying_ts(sUndPtr);*/
      FREE_STOCHVOLCALIB_MEMORY;
      return err;
    }

    for (i = 0; i < iNumInstruments; i++)
      pdPrice[i] = pdTempPrice[i];

    if (iolist)
      err = srt_f_IOstructfree(&iolist);

    pdTempPrice = NULL;
    /*free_underlying_ts(sUndPtr);*/
    if (err) {
      FREE_STOCHVOLCALIB_MEMORY;
      return serror("error in freeing iolist");
    }

    /* pdPrice the deals with the stochastic volatility model */
    for (i = 0; i < iNumInstruments; i++) {
      /* Update stochastic volatility model pdPrices */
      beta[i] = pdPrice[i];

      /* Computing lambda        , the "direction of the gradient" in our
       * algorithm*/
      if (m == 0) {
        lambda[i] = 1;
      } else if (m == 1) {
        lambda[i] = (pdPrice_hat[i] - previouspdPrice_hat[i]) /
                    (beta[i] - previousbeta[i]);
      } else
      /* we force the lambda to decrease */
      {
        if ((pdPrice_hat[i] - previouspdPrice_hat[i]) /
                (beta[i] - previousbeta[i]) <=
            lambda[i]) {
          lambda[i] = (pdPrice_hat[i] - previouspdPrice_hat[i]) /
                      (beta[i] - previousbeta[i]);
        }
      }

      /* in the case lambda <0        , the monte carlo error is too big to
      provide any info so we stop the algorithm by keeping the lambda at 0 */

      if (lambda[i] < 0) {
        lambda[i] = 0;
      }

      /* Update previousbeta and previouspdPrice_hat */
      previousbeta[i] = beta[i];
      previouspdPrice_hat[i] = pdPrice_hat[i];

      /* Update implied "market pdPrices" using the fixed point principle */
      pdPrice_hat[i] =
          pdPrice_hat[i] - lambda[i] * (beta[i] - pdMarketPrices[i]);

      if (pdPrice_hat[i] < 0)
        pdPrice_hat[i] =
            pdPrice_hat[i] + lambda[i] * (beta[i] - pdMarketPrices[i]);

      /* Compute the BS implied volatilities for the stochastic volatility model
       */
      GlobalTheoPrices[i + 1] = beta[i];
      err = swp_f_SwaptionImpliedVol_SwapDP(
          beta[i], &psSwapDp[i], pdStrike[i], peRecPay[i], pszRefRateCode[i],
          yc_name, SRT_LOGNORMAL, &impl_vol_stc[i]);
      if (err) {
        FREE_STOCHVOLCALIB_MEMORY;
        return serror("error in implied vol");
      }

      /* Compute stochastic vol residual error for mth iteration*/
      stoch_model_resid_error[i] = impl_vol_stc[i] - impl_vol_mkt[i];
      stoch_sum_error = stoch_sum_error +
                        stoch_model_resid_error[i] * stoch_model_resid_error[i];

    } /* end "for i" loop on instruments - stochastic*/

    stoch_sum_error = sqrt(stoch_sum_error / iNumInstruments);
    stoch_error[m] = stoch_sum_error;

    smessage("Stochastic volatility model calibration. Finish iteration: %d",
             m);

    /* Check whether convergence criterion is satisfied or maximum number
            of iterations has been reached */

    /*	if (stoch_sum_error <= deter_sum_error)*/
    if (stoch_sum_error <= 0.000001) {
      smessage("Convergence reached at m = %d iteration", m);
      smessage("Stochastic volatility model: CHI = %.8f", stoch_sum_error);
      m = max_iter + 1;
    } else if (m == max_iter) {
      smessage("Convergence criterion not satisfied at max_iter = %d iteration",
               m);
      smessage("Stochastic volatility model: XI =  %.8f", stoch_sum_error);
      m++;
    } else {
      smessage("at max_iter = %d iteration", m);
      smessage("Stochastic volatility model: XI =  %.8f", stoch_sum_error);
      stoch_sum_error = 0;
      m++;
    }

  } /* end of "m" iterations */

  /* Clear memory*/
  FREE_STOCHVOLCALIB_MEMORY;

  /* Returns the calibrated underlying */
  return err;

} /* END Err srt_f_calib_fixed_point(...)  */
