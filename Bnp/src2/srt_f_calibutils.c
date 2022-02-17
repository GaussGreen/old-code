/* ------------------------------------------------------------------------
        FILE NAME: srt_f_calibutils.c

        PURPOSE: utility routines for a calibration

        AUTHORS: O. Van Eyseren  , A. Savine
   ----------------------------------------------------------------------- */

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_calibutils.h"

#define SWAP(a, b)                                                             \
  {                                                                            \
    double tempr;                                                              \
    tempr = (a);                                                               \
    (a) = (b);                                                                 \
    (b) = tempr;                                                               \
  }

/* A numerically big number to prevent a search in a wrong direction */
#define SRT_BIG 1.0e+20

/* -----------------------------------------------------------------------
              FUNCTIONS TO FIND SIGMA OR TAU DATES THAT ARE
                    RELEVANT FOR A SET OF INSTRUMENTS
   ----------------------------------------------------------------------- */

Err find_dates_from_deals(SwapDP *sdparray, StructType *opt_type,
                          long numInstruments, SrtCcyParam *ccy_param,
                          double **sigDates, long *lNumSigmas,
                          double **tauDates, long *lNumTaus) {

  long i, j;
  DateList date_list;
  long TempNumTaus;
  long last_fixing;
  long last_period_start;
  double *temp_sig_dates = NULL;
  double *temp_tau_dates = NULL;

  /* Allocate a temporary array to stroe all dates */
  temp_sig_dates = dvector(0, numInstruments - 1);
  if (!temp_sig_dates)
    return serror("Allocation failure in find_sigdates_from_deals");
  temp_tau_dates = dvector(0, numInstruments - 1);
  if (!temp_tau_dates)
    return serror("Allocation failure in find_sigdates_from_deals");

  /* Extract the relevant fixing date for each deal */
  for (i = 0; i < numInstruments; i++) {
    last_fixing = add_unit(sdparray[i].end, -sdparray[i].spot_lag, SRT_BDAY,
                           MODIFIED_SUCCEEDING);
    temp_tau_dates[i] = last_fixing;

    /* A cap: gets the start date for the last period : one before last date */
    if ((opt_type[i] == CAPFLOOR) || (opt_type[i] == RESETCAPFLOOR)) {
      date_list =
          SwapDP_to_DateList(&(sdparray[i]), ccy_param->cash_bus_day_conv);
      last_period_start = date_list.date[date_list.len - 2];
      last_fixing = add_unit(last_period_start, -sdparray[i].spot_lag, SRT_BDAY,
                             MODIFIED_SUCCEEDING);
      temp_sig_dates[i] = last_fixing;

      free(date_list.date);
    } else
    /* A swaption : extract the fixing date of the underlying swap */
    {
      last_fixing = add_unit(sdparray[i].start, -sdparray[i].spot_lag, SRT_BDAY,
                             MODIFIED_SUCCEEDING);
      temp_sig_dates[i] = last_fixing;
    }
  } /* END of loop on calibration instruments */

  /* Sorts sig dates: from smallest (i=0) to biggest (I=numInstruments-1) */
  for (i = 0; i < numInstruments; i++) {
    for (j = i + 1; j < numInstruments; j++) {
      if (temp_sig_dates[j] < temp_sig_dates[i]) {
        SWAP(temp_sig_dates[j], temp_sig_dates[i]);
      }
    }
  }

  /* Sorts tau dates: from smallest (i=0) to biggest (I=numInstruments-1) */
  for (i = 0; i < numInstruments; i++) {
    for (j = i + 1; j < numInstruments; j++) {
      if (temp_tau_dates[j] < temp_tau_dates[i]) {
        SWAP(temp_tau_dates[j], temp_tau_dates[i]);
      }
    }
  }

  /* Sets initial value for lNumSigmas */
  *lNumSigmas = numInstruments;
  TempNumTaus = numInstruments;

  /* Removes sig dates that are the same */
  for (i = 0; i < *lNumSigmas - 1; i++) {
    while ((temp_sig_dates[i] == temp_sig_dates[i + 1]) &&
           (i < *lNumSigmas - 1)) {

      for (j = i + 1; j < *lNumSigmas; j++) {
        temp_sig_dates[j - 1] = temp_sig_dates[j];
      }
      (*lNumSigmas)--;
    }
  }

  /* Removes tau dates that are the same */
  for (i = 0; i < TempNumTaus - 1; i++) {
    while ((temp_tau_dates[i] == temp_tau_dates[i + 1]) &&
           (i < TempNumTaus - 1)) {

      for (j = i + 1; j < TempNumTaus; j++) {
        temp_tau_dates[j - 1] = temp_tau_dates[j];
      }
      TempNumTaus--;
    }
  }

  /* If there is not at least two instruments with the same sig dates we don't
   * need to find tau */
  if (TempNumTaus + (*lNumSigmas) > numInstruments)
    (*lNumTaus) = numInstruments - (*lNumSigmas);
  else
    *lNumTaus = TempNumTaus;

  /* Allocates sigDates accordingly and fills it */
  *sigDates = dvector(0, *lNumSigmas - 1);
  for (i = 0; i < *lNumSigmas; i++) {
    (*sigDates)[i] = temp_sig_dates[i];
  }

  /* Allocates tauDates accordingly and fills it */
  if (*lNumTaus) {
    *tauDates = dvector(0, *lNumTaus - 1);
    for (i = *lNumTaus - 1; i >= 0; i--)
      (*tauDates)[i] = temp_tau_dates[i];
  }

  /* Free whatever has to be freed */
  if (temp_sig_dates)
    free_dvector(temp_sig_dates, 0, numInstruments - 1);
  temp_sig_dates = NULL;
  if (temp_tau_dates)
    free_dvector(temp_tau_dates, 0, numInstruments - 1);
  temp_tau_dates = NULL;

  /* Return a success message */
  return NULL;

} /* END of find_sigdates_from_deals */

Err find_sigdates_from_deals(SwapDP *sdparray, StructType *opt_type,
                             long numInstruments, SrtCcyParam *ccy_param,
                             double **sigDates, long *lNumSigmas) {

  long i, j;
  DateList date_list;
  long last_fixing;
  long last_period_start;
  double *temp_dates;

  /* Allocate a temporary array to stroe all dates */
  temp_dates = dvector(0, numInstruments - 1);
  if (!temp_dates)
    return serror("Allocation failure in find_sigdates_from_deals");

  /* Extract the relevant fixing date for each deal */
  for (i = 0; i < numInstruments; i++) {
    /* A cap: gets the start date for the last period : one before last date */
    if ((opt_type[i] == CAPFLOOR) || (opt_type[i] == RESETCAPFLOOR)) {
      date_list =
          SwapDP_to_DateList(&(sdparray[i]), ccy_param->cash_bus_day_conv);
      last_period_start = date_list.date[date_list.len - 2];
      last_fixing = add_unit(last_period_start, -sdparray[i].spot_lag, SRT_BDAY,
                             MODIFIED_SUCCEEDING);
      temp_dates[i] = last_fixing;

      free(date_list.date);
    } else
    /* A swaption : extract the fixing date of the underlying swap */
    {
      last_fixing = add_unit(sdparray[i].start, -sdparray[i].spot_lag, SRT_BDAY,
                             MODIFIED_SUCCEEDING);
      temp_dates[i] = last_fixing;
    }
  } /* END of loop on calibration instruments */

  /* Sorts dates: from smallest (i=0) to biggest (I=numInstruments-1) */
  for (i = 0; i < numInstruments; i++) {
    for (j = i + 1; j < numInstruments; j++) {
      if (temp_dates[j] < temp_dates[i])
        SWAP(temp_dates[j], temp_dates[i]);
    }
  }

  /* Sets initial value for lNumSigmas */
  *lNumSigmas = numInstruments;

  /* Removes dates that are the same */
  for (i = 0; i < *lNumSigmas - 1; i++) {
    while ((temp_dates[i] == temp_dates[i + 1]) && (i < *lNumSigmas - 1)) {
      for (j = i + 1; j < *lNumSigmas; j++) {
        temp_dates[j - 1] = temp_dates[j];
      }
      (*lNumSigmas)--;
    }
  }
  /* Allocates sigDates accordingly and fills it */
  *sigDates = dvector(0, *lNumSigmas - 1);
  for (i = 0; i < *lNumSigmas; i++) {
    (*sigDates)[i] = temp_dates[i];
  }

  /* Free whatever has to be freed */
  free_dvector(temp_dates, 0, numInstruments - 1);

  /* Return a success message */
  return NULL;

} /* END of find_sigdates_from_deals */

/* -----------------------------------------------------------------------  */

Err find_taudates_from_deals(SwapDP *sdparray, long numInstruments,
                             SrtCcyParam *ccy_param, double **tauDates,
                             long *lNumTaus) {

  long i, j;
  long last_fixing;
  double *temp_dates = NULL;

  /* Allocate a temporary array to stroe all dates */
  temp_dates = dvector(0, numInstruments - 1);
  if (!temp_dates)
    return serror("Allocation failure in find_taudates_from_deals");

  /* Extract the relevant fixing date for each deal */
  for (i = 0; i < numInstruments; i++) {
    last_fixing = add_unit(sdparray[i].end, -sdparray[i].spot_lag, SRT_BDAY,
                           MODIFIED_SUCCEEDING);
    temp_dates[i] = last_fixing;
  }

  /* Sorts dates: from smallest (i=0) to biggest (I=numInstruments-1) */
  for (i = 0; i < numInstruments; i++) {
    for (j = i + 1; j < numInstruments; j++) {
      if (temp_dates[j] < temp_dates[i])
        SWAP(temp_dates[j], temp_dates[i]);
    }
  }

  /* Sets initial value for lNumTaus */
  *lNumTaus = numInstruments;

  /* Removes dates that are the same */
  for (i = 0; i < *lNumTaus - 1; i++) {
    while ((temp_dates[i] == temp_dates[i + 1]) && (i < *lNumTaus - 1)) {
      for (j = i + 1; j < *lNumTaus; j++) {
        temp_dates[j - 1] = temp_dates[j];
      }
      (*lNumTaus)--;
    }
  }
  /* Allocates tauDates accordingly and fills it */
  *tauDates = dvector(0, *lNumTaus - 1);
  for (i = 0; i < *lNumTaus; i++) {
    (*tauDates)[i] = temp_dates[i];
  }

  /* Free whatever has to be freed */
  free_dvector(temp_dates, 0, numInstruments - 1);

  /* Return a success message */
  return NULL;

} /* END of find_taudates_from_deals */

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
              FUNCTIONS TO FIND A GOOD STARTING POINT FOR THE
              CALIBRATION ALGORITHM USING A SOBOL PRESAMPLING
   ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
        Given a set of parameters[1..n] and a range in which they can evolve  ,
        minimises the function (*funcs)(double*) through a Sobol random
        draw with npaths.
        For calibration  , (*funcs) should be the global criteria
        -----------------------------------------------------------------------
 */

Err sobol_startingpoint(double *param, double **param_bound, long numParams,
                        long npath, double (*funcs)(double *),
                        double *criteria) {
  Err err = NULL;
  long i, j;
  double *sobol;
  double cur_criteria;
  double min_criteria;
  double *temp_param;

  sobol = dvector(1, numParams);
  if (!sobol)
    return serror("Memory allocation failure (1) in sobol_startingpoint");
  temp_param = dvector(1, numParams);
  if (!temp_param) {
    if (sobol)
      free_dvector(sobol, 1, numParams);
    return serror("Memory allocation failure (2) in sobol_startingpoint");
  }

  /* Copy the current parameters in a temporary array */
  for (i = 1; i <= numParams; i++) {
    temp_param[i] = param[i];
  }
  min_criteria = SRT_BIG;

  /* Initialise the Sobol sequence */
  err = sobol_init(1, npath, 0, 0, 1, numParams);
  if (err) {
    if (temp_param)
      free_dvector(temp_param, 1, numParams);
    if (sobol)
      free_dvector(sobol, 1, numParams);
    return serror("Error in SobInit");
  }

  smessage("Starting Sobol search for starting point ");
  smessage("");

  /* Runs the paths and keep the one that has the lowest criteria for (*funcs)
   */
  for (i = 1; i <= npath; i++) {
    smessage(" sobol path: %d", (int)i);
    err = sobol_vector(sobol, 1, numParams);
    if (err) {
      if (temp_param)
        free_dvector(temp_param, 1, numParams);
      if (sobol)
        free_dvector(sobol, 1, numParams);
      return err;
    }

    for (j = 1; j <= numParams; j++) {
      temp_param[j] = param_bound[1][j] +
                      (param_bound[2][j] - param_bound[1][j]) * sobol[j];
    }

    /* Computes the error with the current set of parameters */
    cur_criteria = funcs(temp_param);

    if (cur_criteria <= min_criteria) {
      min_criteria = cur_criteria;
      memcpy(&param[1], &temp_param[1], numParams * sizeof(double));
    }
    smessage(" best criteria so far: %.8f", min_criteria);
  }

  /* Free the variables used in Sobol */
  err = sobol_free();

  /* Free whatever has to be freed */
  free_dvector(sobol, 1, numParams);
  free_dvector(temp_param, 1, numParams);

  /* Sets the minimum criteria obtained so far */
  *criteria = min_criteria;

  /* Return a success message */
  smessage("");
  return NULL;

} /* END sobol_startingpoint */

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
              FUNCTIONS TO MOVE FROMSIG-TAU TO PARAM
                         (AND VICE-VERSA)
   ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
        From the input array of vols and taus (with dates)  , sets optimisation
        paremeters in a vector as required by Numerical Recipes in C
        All the arrays should already be initialised.
        The SigmaValues array should look like this (in input):
                          [0]	[1]		[2]		[3]
   [4]		[5] Date	Sig1	Beta1	Sig2	Beta2	Rho The
   TauValues array should look like this (in input): [0]	[1]
   [2] Date	Tau1	Tau2 The initial set of calibration parameters should
   look like this (in output): [ Sigmas ][ Tau ] [ Beta ]  [ Omega] [ Alpha] [
   Gamma ] [ Rho ]
    ----------------------------------------------------------------------- */

/* AS: for one factor  , eCalibType is supposed to be blind; forcing to blind in
   one factor case is done in srt_f_CalibrateAll
 */

Err from_sigtau_to_optparam(double **ppdSigmaValues, long lNumSigmas,
                            double **ppdTauValues, long lNumUsedTaus,
                            long lNumUsedBetas, long lNumUsedOmegas,
                            SrtCalibType eCalibType, SrtMdlDim eModelDim,
                            double *pdOptParams /* From [1] to [param_number] */
) {
  long i;

  /* The Sigma TS is the first set of calibration parameters */
  for (i = 0; i < lNumSigmas; i++) {
    pdOptParams[i + 1] = ppdSigmaValues[1][i];
  }

  /* The Tau TS determinses the second set of calibration parameters */
  for (i = 0; i < lNumUsedTaus; i++) {
    pdOptParams[lNumSigmas + i + 1] = 1.00 / ppdTauValues[1][i];
  }

  /* The third set of calibration parameters is given by BEta (power) */
  for (i = 0; i < lNumUsedBetas; i++) {
    pdOptParams[lNumSigmas + lNumUsedTaus + i + 1] =
        1.00 / ppdSigmaValues[2][i];
  }

  /* For the Two factor model  , more parameters might be calibrated */
  if (eModelDim == TWO_FAC) {
    /* The fourth set of calibration parameters is given by Omega (diff in
     * Betas) */
    for (i = 0; i < lNumUsedOmegas; i++) {
      pdOptParams[lNumSigmas + lNumUsedTaus + lNumUsedBetas + i + 1] =
          ppdSigmaValues[4][i] - ppdSigmaValues[2][i];
    }

    /* For a Global calibration  , the correlation has to be calibrated */
    if (eCalibType == GLOBAL_CALIB) {
      /* Alpha ( = Sig2/Sig1 ) */
      pdOptParams[lNumSigmas + lNumUsedTaus + lNumUsedBetas + lNumUsedOmegas +
                  1] = ppdSigmaValues[3][0] / ppdSigmaValues[1][0];

      /* Gamma  ( = Lam2 - Lam1) */
      pdOptParams[lNumSigmas + lNumUsedTaus + lNumUsedBetas + lNumUsedOmegas +
                  2] = 1.00 / ppdTauValues[2][0] - 1.00 / ppdTauValues[1][0];

      /* Rho */
      pdOptParams[lNumSigmas + lNumUsedTaus + lNumUsedBetas + lNumUsedOmegas +
                  3] = ppdSigmaValues[5][0];
    }
  }

  return NULL;

} /* END Err from_sigtau_to_optparam(...) */

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
        From the optimisation paremeters vector (as required by Numerical
        Recipes in C) sets the corresponding TermStrcuture of vols and taus
    in a single array ts_values[7][...]
        The initial set of calibration parameters should look like this (in
   input): [ Sigmas ][ Tau ] [ Beta ]  [ Omega] [ Alpha] [ Gamma ] [ Rho ] The
   SigmaValues array should look like this (in output): [0]	[1]
   [2]		[3]		[4]		[5] Date	Sig1	Beta1
   Sig2	Beta2	Rho The TauValues array should look like this (in output): [0]
   [1]    [2] Date	 Tau1	Tau2
        -----------------------------------------------------------------------
 */

Err from_optparam_to_sigtau(
    double *pdOptParams, /* From [1] to [param_number] */
    double **ppdSigmaValues, long lNumSigmas, double **ppdTauValues,
    long lNumTaus, SRT_Boolean bFreezeTau, SRT_Boolean bOneTau,
    double *dFixedTau, SRT_Boolean bFreezeBeta, SRT_Boolean bOneBeta,
    double dFixedBeta, SrtMdlDim eModelDim, SRT_Boolean bFreezeOmega,
    SRT_Boolean bOneOmega, double dFixedOmega, SrtCalibType eCalibType,
    double dFixedAlpha, double dFixedGamma, double dFixedRho) {
  long i;
  double dAlpha, dGamma, dRho;
  long lNumUsedTaus;
  long lNumUsedBetas;
  long lNumUsedOmegas;

  /* The Sigma 1 term structure is stored in the first set of calibration
   * parameters */
  for (i = 0; i < lNumSigmas; i++) {
    ppdSigmaValues[1][i] = pdOptParams[i + 1];
  }

  /* The Tau (=1/Lambda) Term structure is stored in the second set of
   * calibration parameters */
  if (bFreezeTau == SRT_YES) {
    /* Tau has not been optimised: it is the Fixed Tau */
    lNumUsedTaus = 0;
    for (i = 0; i < lNumTaus; i++) {
      ppdTauValues[1][i] = dFixedTau[i];
    }
  } else if (bOneTau == SRT_YES) {
    /* Only a Single Tau has been optimised */
    lNumUsedTaus = 1;
    for (i = 0; i < lNumTaus; i++) {
      ppdTauValues[1][i] = 1.0 / pdOptParams[lNumSigmas + 1];
    }
  } else if (bOneTau == SRT_NO) {
    /* The Entire tau term Strucutre has been optimised */
    lNumUsedTaus = lNumTaus;
    for (i = 0; i < lNumTaus; i++) {
      ppdTauValues[1][i] = 1.0 / pdOptParams[lNumSigmas + i + 1];
    }
  }

  /* The Beta 1 term structure is stored in the third set of calibration
   * parameters */
  if (bFreezeBeta == SRT_YES) {
    /* Beta has not been optimised: it is the Fixed Beta */
    lNumUsedBetas = 0;
    for (i = 0; i < lNumSigmas; i++) {
      ppdSigmaValues[2][i] = dFixedBeta;
    }
  } else if (bOneBeta == SRT_YES) {
    /* Only a Single Beta has been optimised */
    lNumUsedBetas = 1;
    for (i = 0; i < lNumSigmas; i++) {
      ppdSigmaValues[2][i] = pdOptParams[lNumSigmas + lNumUsedTaus + 1];
    }
  } else if (bOneBeta == SRT_NO) {
    /* The Entire Beta Term Structure has been optimised */
    lNumUsedBetas = lNumSigmas;
    for (i = 0; i < lNumSigmas; i++) {
      ppdSigmaValues[2][i] = pdOptParams[lNumSigmas + lNumUsedTaus + i + 1];
    }
  }

  /* For Two factor models: correlation parameters */
  if (eModelDim == TWO_FAC) {
    /* The Beta 2 ( = Beta 1 + Omega) term structure is stored in the fourth set
     * of calibration parameters */
    if (bFreezeOmega == SRT_YES) {
      /* Omega has not been optimised: it is the Fixed Omega */
      lNumUsedOmegas = 0;
      for (i = 0; i < lNumSigmas; i++) {
        ppdSigmaValues[4][i] = ppdSigmaValues[2][i] + dFixedBeta;
      }
    } else if (bOneOmega == SRT_YES) {
      /* Only a Single Omega has been optimised */
      lNumUsedOmegas = 1;
      for (i = 0; i < lNumSigmas; i++) {
        ppdSigmaValues[4][i] =
            ppdSigmaValues[2][i] +
            pdOptParams[lNumSigmas + lNumUsedTaus + lNumUsedBetas + 1];
        ;
      }
    } else if (bOneOmega == SRT_NO) {
      /* The Entire Omega Term Strucutre has been optimised */
      lNumUsedOmegas = lNumSigmas;
      for (i = 0; i < lNumSigmas; i++) {
        ppdSigmaValues[4][i] =
            ppdSigmaValues[2][i] +
            pdOptParams[lNumSigmas + lNumUsedTaus + lNumUsedBetas + i + 1];
        ;
      }
    }

    /* The Correlation Parameters (ALpha  , Gamma and Rho) */
    if (eCalibType == GLOBAL_CALIB) {
      /* The parameters have been calibrated */
      dAlpha = pdOptParams[lNumSigmas + lNumUsedTaus + lNumUsedBetas +
                           lNumUsedOmegas + 1];
      dGamma = pdOptParams[lNumSigmas + lNumUsedTaus + lNumUsedBetas +
                           lNumUsedOmegas + 2];
      dRho = pdOptParams[lNumSigmas + lNumUsedTaus + lNumUsedBetas +
                         lNumUsedOmegas + 3];

      /* Rebuilds the Sigma 2 and Rho Term Structures */
      for (i = 0; i < lNumSigmas; i++) {
        ppdSigmaValues[3][i] = dAlpha * pdOptParams[i + 1];
        ppdSigmaValues[5][i] = dRho;
      }

      /* Rebuilds the Tau 2 (Lam2 = Lam 1 + Gamma) Term Struct */
      for (i = 0; i < lNumTaus; i++) {
        ppdTauValues[2][i] = 1.0 / (dGamma + 1.0 / ppdTauValues[1][i]);
      }
    } else if (eCalibType == FIXED_CALIB) {
      /* Rebuilds the Sigma 2 and Rho Term Structures using the Fixed Alpha and
       * Rho*/
      for (i = 0; i < lNumSigmas; i++) {
        ppdSigmaValues[3][i] = dFixedAlpha * pdOptParams[i + 1];
        ppdSigmaValues[5][i] = dFixedRho;
      }

      /* Rebuilds the Tau 2 (Lam2 = Lam 1 + Gamma) Term Struct using the Fixed
       * Gamma */
      for (i = 0; i < lNumTaus; i++) {
        ppdTauValues[2][i] = 1.0 / (dFixedGamma + 1.0 / ppdTauValues[1][i]);
      }
    } /* END if (eCalibtype == FIXED_CALIB ) */
  }

  return NULL;

} /* END Err from_optparam_to_sigtau(...) */

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
              FUNCTIONS TO SET UP PARAMETERS NEEDED BY CALIBRATION
   ----------------------------------------------------------------------- */
/* -----------------------------------------------------------------------
        Check if parameters (pdOptimParams goes from [1] to [lNumParams]) are
        within the min/max ppdParamBounds ([0]: min; [1] : max)
   ----------------------------------------------------------------------- */

SRT_Boolean are_calib_parameters_within_band(double *pdOptimParams,
                                             long lNumParams,
                                             double **ppdParamBounds) {
  long i;

  for (i = 1; i <= lNumParams; i++) {
    if ((pdOptimParams[i] < ppdParamBounds[1][i]) ||
        (pdOptimParams[i] > ppdParamBounds[2][i]))
      return SRT_NO;
  }

  return SRT_YES;
}

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
        From the calibration parameters (default or input on spreadsheet)  ,
        gives the ppdParamBounds of the parameters in an array (for Sobol):
                [1][1..n] is the minimum value of the parameter
                [2][1..n] is the maximum value of the parameter
        -----------------------------------------------------------------------
 */

Err set_calib_parameters_bounds(SrtCalibParam *psCalibParams,
                                double **ppdParamBounds, long lNumSigmas,
                                long lNumUsedTaus, long lNumUsedBetas,
                                long lNumUsedOmegas, SrtMdlDim eModelDim) {

  int i;

  if (!ppdParamBounds)
    return serror("Array not defined in set_calib_parameters_bounds");

  /* Bounds for Sigma */
  for (i = 1; i <= lNumSigmas; i++) {
    ppdParamBounds[1][i] = psCalibParams->dSigmaMin;
    ppdParamBounds[2][i] = psCalibParams->dSigmaMax;
  }

  /* Bounds for Lambda */
  for (i = 1; i <= lNumUsedTaus; i++) {
    ppdParamBounds[1][lNumSigmas + i] = psCalibParams->dLambdaMin;
    ppdParamBounds[2][lNumSigmas + i] = psCalibParams->dLambdaMax;
  }

  /* Bounds for Beta */
  for (i = 1; i <= lNumUsedBetas; i++) {
    ppdParamBounds[1][lNumSigmas + lNumUsedTaus + i] = psCalibParams->dBetaMin;
    ppdParamBounds[2][lNumSigmas + lNumUsedTaus + i] = psCalibParams->dBetaMax;
  }

  /* For two factor model only */
  if (eModelDim == TWO_FAC) {
    /* Bounds for Omega */
    for (i = 1; i <= lNumUsedOmegas; i++) {
      ppdParamBounds[1][lNumSigmas + lNumUsedTaus + lNumUsedBetas + i] =
          psCalibParams->dBetaMin;
      ppdParamBounds[2][lNumSigmas + lNumUsedTaus + lNumUsedBetas + i] =
          psCalibParams->dBetaMax;
    }

    /* The correlation parameters have to be calibrated too when GLOBAL */
    if (psCalibParams->eCalibType == GLOBAL_CALIB) {
      /* Alpha */
      ppdParamBounds[1][lNumSigmas + lNumUsedTaus + lNumUsedBetas +
                        lNumUsedOmegas + 1] = psCalibParams->dAlphaMin;
      ppdParamBounds[2][lNumSigmas + lNumUsedTaus + lNumUsedBetas +
                        lNumUsedOmegas + 1] = psCalibParams->dAlphaMax;

      /* Beta */
      ppdParamBounds[1][lNumSigmas + lNumUsedTaus + lNumUsedBetas +
                        lNumUsedOmegas + 2] = psCalibParams->dBetaMin;
      ppdParamBounds[2][lNumSigmas + lNumUsedTaus + lNumUsedBetas +
                        lNumUsedOmegas + 2] = psCalibParams->dBetaMax;

      /* Rho */
      ppdParamBounds[1][lNumSigmas + lNumUsedTaus + lNumUsedBetas +
                        lNumUsedOmegas + 3] = psCalibParams->dRhoMin;
      ppdParamBounds[2][lNumSigmas + lNumUsedTaus + lNumUsedBetas +
                        lNumUsedOmegas + 3] = psCalibParams->dRhoMax;
    }

  } /* END if (eModelDim == TWO_FAC) */

  return NULL;

} /* END Err set_calib_parameters_bounds(..) */

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
        Sets pdMarketWeights for minimisation criteria as in LM for
   SIMPLEX/ANNEALING criteria = Sum(0<=i<n) [(y_i-y0_i)/pdMarketWeights_i]^2 For
   Option prices  , pdMarketWeights = pdMktVega * sqrt( n * dOptionsWeight ).
        For correlation  , pdMarketWeights = n / sqrt( 1- dOptionsWeight )
        The global criteria is therefore:
                                   dOptionsWeight^2 * 1/n * Sum(i)
   [(y_i-y0_i)/vega_i]^2
                                + (1-frd_weght)^2 * (1/N)^2 * Sum(j) [corr_i -
   mktcorr_i]^2
    ----------------------------------------------------------------------- */
Err define_criteria_weights(
    long lNumInstruments, long lNumTenors, SRT_Boolean bSmoothSigma,
    double dOptionsWeight, double dCorrelWeight, double bSmoothSigmaWeight,
    double *pdMktVega,       /* From [0] to [lNumInstruments-1] */
    double *pdMarketWeights) /* From [1] to [lNumData] */
{
  long i, lNumCorrel;

  /* If smooth sigma  , update options and correl weights */
  if (bSmoothSigma == SRT_YES) {
    dOptionsWeight *= 1.0 - bSmoothSigmaWeight;
    dCorrelWeight *= 1.0 - bSmoothSigmaWeight;
  }

  /* Market instruments weights */
  for (i = 1; i <= lNumInstruments; i++) {
    if (dOptionsWeight > 0.0) {
      pdMarketWeights[i] =
          pdMktVega[i - 1] * sqrt((double)lNumInstruments / dOptionsWeight);
    } else {
      pdMarketWeights[i] = SRT_BIG;
    }
  }

  /* Correlation weights */
  lNumCorrel = lNumTenors * lNumTenors;
  for (i = 1; i <= lNumCorrel; i++) {
    if (dCorrelWeight > 0.0) {
      pdMarketWeights[lNumInstruments + i] =
          sqrt((double)lNumCorrel / dCorrelWeight);
    } else {
      pdMarketWeights[lNumInstruments + i] = SRT_BIG;
    }
  }

  /* Smooth Sigma */
  if (bSmoothSigma == SRT_YES) {
    if (bSmoothSigmaWeight > 0.0) {
      pdMarketWeights[lNumInstruments + lNumCorrel + 1] =
          1.0 / sqrt(bSmoothSigmaWeight);
    } else {
      pdMarketWeights[lNumInstruments + lNumCorrel + 1] = SRT_BIG;
    }
  }

  return NULL;

} /* Err define_criteria_weights(...) */

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
        Sets a vector of market target for minimisation criteria:
                if (index)<=lNumInstruments:
                        set market pdMktPrice of caliration instrument
                if (index)>lNumInstruments:
                        set market correlation
    ----------------------------------------------------------------------- */

Err define_market_targets(
    long lNumInstruments, SRT_Boolean bSmoothSigma,
    double *mkt_price,             /* From [0] to [lNumInstruments-1] */
    double **ppdCorrelationMatrix, /* From [0] to [lNumData-1] */
    long lNumTenors, long lNumData,
    double *pdMarketTargets) /* From [1] to [lNumData] */
{
  long i;
  long j;

  /* mkt instruments */
  for (i = 1; i <= lNumInstruments; i++) {
    pdMarketTargets[i] = mkt_price[i - 1];
  }
  /* correl */
  for (i = 0; i < lNumTenors; i++) {
    for (j = 0; j < lNumTenors; j++) {
      pdMarketTargets[lNumInstruments + i * lNumTenors + j + 1] =
          ppdCorrelationMatrix[i][j];
    }
  }
  /* smooth sig */
  if (bSmoothSigma == SRT_YES) {
    pdMarketTargets[lNumInstruments + lNumTenors * lNumTenors + 1] = 0.0;
  }

  return NULL;

} /* END Err define_market_targets(...) */

/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------

          FOR QUICK TRANSFORMATIONS FORM A MODEL TS TO AN LGM ONE

   ----------------------------------------------------------------------------
 */

/* From a FULL model SigmaCurve to a FULL LGM one */
Err from_model_ts_to_LGM_ts(SrtMdlDim eModelDim, double **ppdModelSigmaCurve,
                            long lNumSigmas, String szYieldCurveName,
                            double **ppdLgmSigmaCurve) {
  Err err = NULL;
  double dFra;
  double dBeta;
  int i;

  /* Copies the Model Sigma Curve into the LGM One */
  for (i = 0; i < lNumSigmas; i++) {
    /* Date */
    ppdLgmSigmaCurve[0][i] = ppdModelSigmaCurve[0][i];

    /* Sigma */
    ppdLgmSigmaCurve[1][i] = ppdModelSigmaCurve[1][i];

    /* Beta */
    ppdLgmSigmaCurve[2][i] = 0.0;

    if (eModelDim == TWO_FAC) {
      /* Sigma 2 */
      ppdLgmSigmaCurve[3][i] = ppdModelSigmaCurve[3][i];

      /* Beta 2 */
      ppdLgmSigmaCurve[4][i] = 0.0;

      /* Rho */
      ppdLgmSigmaCurve[5][i] = ppdModelSigmaCurve[5][i];
    }
  }

  /* Rebuilds a fake LGM TS by a quick Normal approximation of the Volatility */
  for (i = 0; i < lNumSigmas; i++) {
    /* Computes the FRA at the volatility date */
    dFra = swp_f_zr((Ddate)(ppdModelSigmaCurve[0][i]),
                    (Ddate)(ppdModelSigmaCurve[0][i]) + 14.0, szYieldCurveName);

    /* Gets the Beta */
    dBeta = (Ddate)(ppdModelSigmaCurve[2][i]);

    /* The rough normal vol is the Beta vol times the FRA power Beta */
    ppdLgmSigmaCurve[1][i] *= pow(dFra, dBeta);

    if (eModelDim == TWO_FAC) {
      dBeta = (Ddate)(ppdModelSigmaCurve[4][i]);
      ppdLgmSigmaCurve[3][i] *= pow(dFra, dBeta);
    }
  }

  /* Return a success message */
  return err;

} /* END Err from_model_ts_to_LGM_ts(...) */

/* --------------------------------------------------------------------------------
 */

/* From a FULL (rough) LGM SigmaCurve to a FULL Model one (inverse of previous)
 */

Err from_LGM_ts_to_model_ts(SrtMdlDim eModelDim, double **ppdLgmSigmaCurve,
                            long lNumSigmas, String szYieldCurveName,
                            double **ppdModelSigmaCurve) {
  Err err = NULL;
  double dFra;
  double dBeta;
  int i;

  /* Copies the Lgm Sigma Curve into the Model One */
  for (i = 0; i < lNumSigmas; i++) {
    /* Sigma */
    ppdModelSigmaCurve[1][i] = ppdLgmSigmaCurve[1][i];

    if (eModelDim == TWO_FAC) {
      /* Sigma 2 */
      ppdModelSigmaCurve[3][i] = ppdLgmSigmaCurve[3][i];

      /* Rho */
      ppdModelSigmaCurve[5][i] = ppdLgmSigmaCurve[5][i];
    }
  }

  /* Rebuilds a Model TS by a the inverse of the quick Normal approximation of
   * the Volatility */
  for (i = 0; i < lNumSigmas; i++) {
    /* Computes the FRA at the volatility date */
    dFra = swp_f_zr((Ddate)(ppdModelSigmaCurve[0][i]),
                    (Ddate)(ppdModelSigmaCurve[0][i]) + 14.0, szYieldCurveName);

    /* Gets the Beta */
    dBeta = (Ddate)(ppdModelSigmaCurve[2][i]);

    /* The rough model vol is the normal vol divided the FRA power Beta */
    ppdModelSigmaCurve[1][i] /= pow(dFra, dBeta);

    if (eModelDim == TWO_FAC) {
      dBeta = (Ddate)(ppdModelSigmaCurve[4][i]);
      ppdModelSigmaCurve[3][i] /= pow(dFra, dBeta);
    }
  }

  /* Return a success message */
  return err;

} /* END Err from_model_ts_to_LGM_ts(...) */

/* --------------------------------------------------------------------------------
 */

Err srt_f_SquareRootSymMatrix(double **A, long dim, double **result) {
  Err err = NULL;
  long i, j;
  double *eigen_val;
  double **eigen_vec;

  eigen_val = dvector(0, dim);
  eigen_vec = dmatrix(0, dim, 0, dim);

  err = diagonalise_symmetric_matrix(A, dim, eigen_val, eigen_vec);
  if (err) {
    free_dvector(eigen_val, 0, dim);
    eigen_val = NULL;
    free_dmatrix(eigen_vec, 0, dim, 0, dim);
    eigen_vec = NULL;
    return err;
  }

  for (i = 0; i < dim; i++)
    eigen_val[i] = (eigen_val[i] > 0) ? sqrt(eigen_val[i]) : 0.0;

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      result[i][j] += eigen_vec[i][j] * eigen_val[j];

  if (eigen_val)
    free_dvector(eigen_val, 0, dim);
  eigen_val = NULL;
  if (eigen_vec)
    free_dmatrix(eigen_vec, 0, dim, 0, dim);
  eigen_vec = NULL;

  return err;
}

Err srt_f_HellingerDist(double **HistSquare, int flag1, double **ModelCorr,
                        int flag2, long dim, double *result) {
  Err err = NULL;
  double temp;
  long i, j;
  double **pdSquareRoot1;
  double **pdSquareRoot2;

  /* Flag == 0 means the input is already a Square root Matrix of a Sym Matrix
   */
  pdSquareRoot1 = dmatrix(0, dim, 0, dim);
  if (flag1 == 1) {
    err = srt_f_SquareRootSymMatrix(ModelCorr, dim, pdSquareRoot1);
    if (err) {
      free_dmatrix(pdSquareRoot1, 0, dim, 0, dim);
      pdSquareRoot1 = NULL;
      return err;
    }
  } else {
    for (i = 0; i < dim; i++)
      for (j = 0; j < dim; j++)
        pdSquareRoot1[i][j] = HistSquare[i][j];
  }

  pdSquareRoot2 = dmatrix(0, dim, 0, dim);
  if (flag2 == 1) {
    err = srt_f_SquareRootSymMatrix(ModelCorr, dim, pdSquareRoot2);
    if (err) {
      free_dmatrix(pdSquareRoot1, 0, dim, 0, dim);
      pdSquareRoot1 = NULL;
      free_dmatrix(pdSquareRoot2, 0, dim, 0, dim);
      pdSquareRoot2 = NULL;
      return err;
    }

  } else {
    for (i = 0; i < dim; i++)
      for (j = 0; j < dim; j++)
        pdSquareRoot2[i][j] = ModelCorr[i][j];
  }

  /* Compute the diff of the two matrix and put the result into pdSquareRoot1 */
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      pdSquareRoot1[i][j] -= pdSquareRoot2[i][j];

  /* Compute the Trace of the Square of pdSquareRoot1 */
  temp = 0.0;
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      temp += pdSquareRoot1[i][j] * pdSquareRoot1[j][i];

  (*result) = temp;

  if (pdSquareRoot1)
    free_dmatrix(pdSquareRoot1, 0, dim, 0, dim);
  pdSquareRoot1 = NULL;
  if (pdSquareRoot2)
    free_dmatrix(pdSquareRoot2, 0, dim, 0, dim);
  pdSquareRoot2 = NULL;

  return err;
}
