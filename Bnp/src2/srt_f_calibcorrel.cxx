#include "SrtAccess.h"
#include "srt_h_all.h"
#include "srt_h_calibparams.h"

/*------------------------------------------------------------------------------
                                                DEFINITIONS
-------------------------------------------------------------------------------*/

#define XL_MAX_UND_NAME_LEN 36
#define SRT_BIG 1.0e+20
#define FX_ZERO_SHIFT 1.0e-06
#define FX_PROP_SHIFT 1.0e-03

/* ---------------------------------------------------------------------------
                                        STRUCTURES ' DEFINITIONS
--------------------------------------------------------------------------
*/

static struct {

  long ltoday;
  /* Optimisation Parameters Structures */
  int iNumParams;

  int iNumLocVols;
  long *plFXVolCrvDates;

  int iNumCorrDates;
  double *pdCorrDates; /*Must be *double for srt_f_init_Corr_TermStruct */
  int iNumCorrel;

  /* Underlyings Details */

  SrtModelType eModelType;
  char *szFXUndName;
  char *szDomUndName;
  char *szForUndName;

  String **ppszUndNames;
  /* Needs it for srt_f_init_Corr_TermStruct */

  SRT_Boolean bCalibCorrel;
  /* =1 calibrate the correlation =0 otherwise */

  double **ppdInitCorrels;
  /*One Vector of Correls per CorrDate - Needs It For The Initialisation */

  double dFXSpot;

} STATIC_PARAMS;

static struct {
  int iNumInstrs;
  long *plFXOptMDates;
  double *pdFXOptMats;

  double *pdFXOptMktVols;
  double *pdWeights;

} STATIC_INSTRUMENT;

/*-----------------------------------------------------------------------------
                                                        PROTOTYPES
-------------------------------------------------------------------------------*/

/*----------------------Functions ' Cores
 * -----------------------------------------*/

/* This function set the static structures: STATIC_PARAMS & STATIC_INSTRUMENT */

static Err
set_static_structures_for_fx_calib(char *szFXUndName, long *plFXOptMDates,
                                   double *pdFXOptMktVols, double *pdWeights,
                                   int iNumParams, int iNumLocVols,
                                   double **ppdFXVolCrv, int iNumInstrs) {
  long ltoday;
  SrtUndPtr sFXUnd;

  SrtCorrLstVal *psCorrval;
  SrtCorrLstPtr sCorrlist;
  SrtListAtom *psCorr_lc;
  long i;
  double time;
  Err err = NULL;

  /*----	FILL IN STATIC_INSTRUMENT STRUCTURE -------------------------*/

  STATIC_INSTRUMENT.iNumInstrs = iNumInstrs;

  sFXUnd = lookup_und(szFXUndName);
  ltoday = get_today_from_underlying(sFXUnd);

  STATIC_INSTRUMENT.pdFXOptMats = dvector(0, iNumInstrs - 1);
  STATIC_INSTRUMENT.pdWeights = dvector(1, iNumInstrs);
  STATIC_INSTRUMENT.plFXOptMDates = lvector(0, iNumInstrs - 1);

  if (!(STATIC_INSTRUMENT.pdFXOptMats && STATIC_INSTRUMENT.pdWeights &&
        STATIC_INSTRUMENT.plFXOptMDates)) {
    return serror("STATIC structure list improperly initialised");
  } else {
    for (i = 0; i < iNumInstrs; i++) {
      STATIC_INSTRUMENT.pdWeights[i + 1] = pdWeights[i + 1];
      STATIC_INSTRUMENT.plFXOptMDates[i] = plFXOptMDates[i];
      STATIC_INSTRUMENT.pdFXOptMats[i] =
          (plFXOptMDates[i] - ltoday) * YEARS_IN_DAY;
    }
  }

  STATIC_INSTRUMENT.pdFXOptMktVols = dvector(1, iNumInstrs);
  if (!STATIC_INSTRUMENT.pdFXOptMktVols) {
    return serror("STATIC structure list improperly initialised");
  } else {
    for (i = 1; i <= iNumInstrs; i++) {
      STATIC_INSTRUMENT.pdFXOptMktVols[i] = pdFXOptMktVols[i];
    }
  }
  /*-----------END OF FILL IN STATIC_INSTRUMENT STRUCTURE
   * -------------------------*/

  /*-----------FILL IN STATIC_PARAMS STRUCTURE -------------------------*/

  STATIC_PARAMS.ltoday = ltoday;
  STATIC_PARAMS.iNumParams = iNumParams;
  STATIC_PARAMS.iNumLocVols = iNumLocVols;
  STATIC_PARAMS.plFXVolCrvDates = lngvector(0, iNumLocVols - 1);
  if (!STATIC_PARAMS.plFXVolCrvDates) {
    return serror("STATIC structure list improperly initialised");
  } else {
    for (i = 0; i < iNumLocVols; i++) {
      STATIC_PARAMS.plFXVolCrvDates[i] = (long)ppdFXVolCrv[0][i];
    }
  }

  sCorrlist = srt_f_GetTheCorrelationList();
  if (!sCorrlist->head->element)
    return serror("correlation list improperly initialised");

  /* Get the number of Correlation Dates */
  psCorr_lc = sCorrlist->head;
  STATIC_PARAMS.iNumCorrDates = 0;
  while (psCorr_lc != NULL) {
    STATIC_PARAMS.iNumCorrDates++;
    psCorr_lc = psCorr_lc->next;
  }

  /* Get the Correlation Dates */
  STATIC_PARAMS.pdCorrDates = dvector(0, STATIC_PARAMS.iNumCorrDates - 1);
  if (!STATIC_PARAMS.pdCorrDates)
    return serror("STATIC structure list improperly initialised");

  i = 0;
  psCorr_lc = sCorrlist->head;
  while (psCorr_lc != NULL) {
    psCorrval = (SrtCorrLstVal *)psCorr_lc->element->val.pval;
    STATIC_PARAMS.pdCorrDates[i] = (long)psCorrval->date;
    psCorr_lc = psCorr_lc->next;
    i++;
  }

  /*Get the Number of Correlations */
  psCorr_lc = sCorrlist->head;
  psCorrval = (SrtCorrLstVal *)psCorr_lc->element->val.pval;
  STATIC_PARAMS.iNumCorrel = psCorrval->ncorrel;
  if (STATIC_PARAMS.iNumCorrel != 3)
    return serror("Wrong number of correlation in input.");

  /*Get the Underlyings'Name */

  STATIC_PARAMS.eModelType = get_mdltype_from_fxund(sFXUnd);

  STATIC_PARAMS.szFXUndName = (char *)malloc(32 * sizeof(char));
  strcpy(STATIC_PARAMS.szFXUndName, szFXUndName);

  sFXUnd = lookup_und(szFXUndName);
  if (sFXUnd == NULL)
    return serror("Undefined Underlying %s", szFXUndName);

  STATIC_PARAMS.szDomUndName = (char *)malloc(32 * sizeof(char));
  strcpy(STATIC_PARAMS.szDomUndName, get_domname_from_fxund(sFXUnd));

  STATIC_PARAMS.szForUndName = (char *)malloc(32 * sizeof(char));
  strcpy(STATIC_PARAMS.szForUndName, get_forname_from_fxund(sFXUnd));

  STATIC_PARAMS.ppszUndNames = smatrix_size(0, 2, 0, 1, 128);

  strcpy(STATIC_PARAMS.ppszUndNames[0][0], STATIC_PARAMS.szFXUndName);

  strcpy(STATIC_PARAMS.ppszUndNames[0][1], STATIC_PARAMS.szDomUndName);

  strcpy(STATIC_PARAMS.ppszUndNames[1][0], STATIC_PARAMS.szFXUndName);

  strcpy(STATIC_PARAMS.ppszUndNames[1][1], STATIC_PARAMS.szForUndName);

  strcpy(STATIC_PARAMS.ppszUndNames[2][0], STATIC_PARAMS.szDomUndName);

  strcpy(STATIC_PARAMS.ppszUndNames[2][1], STATIC_PARAMS.szForUndName);

  /*Get the Initial Correlations */

  STATIC_PARAMS.ppdInitCorrels = dmatrix(0, STATIC_PARAMS.iNumCorrel - 1, 0,
                                         STATIC_PARAMS.iNumCorrDates - 1);

  if (!STATIC_PARAMS.ppdInitCorrels)
    return serror("STATIC structure list improperly initialised");

  for (i = 0; i < STATIC_PARAMS.iNumCorrDates; i++) {
    time = (STATIC_PARAMS.pdCorrDates[i] - ltoday) * YEARS_IN_DAY;

    err = srt_f_get_corr_from_CorrList(sCorrlist, STATIC_PARAMS.szDomUndName,
                                       STATIC_PARAMS.szFXUndName, time,
                                       &STATIC_PARAMS.ppdInitCorrels[0][i]);

    err = srt_f_get_corr_from_CorrList(sCorrlist, STATIC_PARAMS.szForUndName,
                                       STATIC_PARAMS.szFXUndName, time,
                                       &STATIC_PARAMS.ppdInitCorrels[1][i]);

    err = srt_f_get_corr_from_CorrList(sCorrlist, STATIC_PARAMS.szForUndName,
                                       STATIC_PARAMS.szDomUndName, time,
                                       &STATIC_PARAMS.ppdInitCorrels[2][i]);
  }

  /* Get the FXUndSpot */

  STATIC_PARAMS.dFXSpot = get_spot_from_fxund(sFXUnd);

  return NULL;

} /* END OF set_static_structures_for_fx_calib */

/* This function gets an fx vol curve from optim parameters */

static void from_optparam_to_fx_vol_curve(double *pdOptimParams,
                                          double **ppdFXVolCrv) {

  int i, j;

  for (i = 0; i < STATIC_PARAMS.iNumLocVols; i++) {
    /* Get the dates */
    ppdFXVolCrv[0][i] = STATIC_PARAMS.plFXVolCrvDates[i];

    /*Get the local vols */
    if (STATIC_PARAMS.bCalibCorrel == SRT_YES) {
      j = i + STATIC_PARAMS.iNumCorrDates * STATIC_PARAMS.iNumCorrel;
      ppdFXVolCrv[1][i] = pdOptimParams[j + 1];
    } else {
      ppdFXVolCrv[1][i] = pdOptimParams[i + 1];
    }
  }

} /* End of from_optparam_to_fx_vol_curve */

/* This function gets the correlation tensor from the optim parameters */

static void from_optparam_to_correl(double *pdOptimParams,
                                    double ***pppdCorrTensor) {

  int i, j, k;

  for (k = 0; k < STATIC_PARAMS.iNumCorrDates; k++) {
    for (i = 0; i < STATIC_PARAMS.iNumCorrel; i++)
      pppdCorrTensor[i][i][k] = 1.0;
  }

  for (k = 0; k < STATIC_PARAMS.iNumCorrDates; k++) {
    for (i = 0; i < STATIC_PARAMS.iNumCorrel; i++) {
      for (j = 0; j < i; j++) {
        pppdCorrTensor[i][j][k] =
            pdOptimParams[k * STATIC_PARAMS.iNumCorrel + i + j];
        pppdCorrTensor[j][i][k] = pppdCorrTensor[i][j][k];
      }
    }
  } /* i=0 DomUnd        , i=1 ForUnd         , i=2 FXUnd */

} /* End of from_optparam_to_correl */

/* This function test if the the FX Vol are positive */
SRT_Boolean are_fx_vol_positive(double **ppdFXVolCrv) {

  SRT_Boolean are_positive;
  long i;

  are_positive = SRT_YES;

  for (i = 0; i < STATIC_PARAMS.iNumLocVols; i++) {
    if (ppdFXVolCrv[1][i] <= 0.0) {
      are_positive = SRT_NO;
      return are_positive;
    }
  }

  return are_positive;

} /*End of are_fx_vol_positive */

/* This function tests if the different correlations matrix are definite
 * positive */

static SRT_Boolean are_correl_matrix_definite_positive(double ***pppdCorrel) {

  long i, j, k;
  double **ppdEigenvector = NULL;
  double *pdEigenval = NULL;
  double **ppdCorrMatrixToTest = NULL;
  SRT_Boolean are_definite_positive = SRT_YES;
  Err err = NULL;

  ppdCorrMatrixToTest = dmatrix(0, (STATIC_PARAMS.iNumCorrel - 1), 0,
                                (STATIC_PARAMS.iNumCorrel - 1));

  ppdEigenvector = dmatrix(0, (STATIC_PARAMS.iNumCorrel - 1), 0,
                           (STATIC_PARAMS.iNumCorrel - 1));
  pdEigenval = dvector(0, (STATIC_PARAMS.iNumCorrel - 1));

  for (k = 0; k < STATIC_PARAMS.iNumCorrDates; k++) {
    for (i = 0; i < STATIC_PARAMS.iNumCorrel; i++) {
      for (j = 0; j < (i + 1); j++) {
        ppdCorrMatrixToTest[i][j] = pppdCorrel[i][j][k];
        ppdCorrMatrixToTest[j][i] = ppdCorrMatrixToTest[i][j];
      } /*End of loop on j */
    }   /*End of loop on i */

    err = diagonalise_symmetric_matrix(ppdCorrMatrixToTest,
                                       STATIC_PARAMS.iNumCorrel, pdEigenval,
                                       ppdEigenvector);
    for (j = 0; j < STATIC_PARAMS.iNumCorrel; j++) {
      if (pdEigenval[j] < 0) {
        if (pdEigenval)
          free_dvector(pdEigenval, 0, STATIC_PARAMS.iNumCorrel - 1);
        pdEigenval = NULL;

        if (ppdEigenvector) {
          free_dmatrix(ppdEigenvector, 0, STATIC_PARAMS.iNumCorrel - 1, 0,
                       STATIC_PARAMS.iNumCorrel - 1);
          ppdEigenvector = NULL;
        }

        if (ppdCorrMatrixToTest) {
          free_dmatrix(ppdCorrMatrixToTest, 0, STATIC_PARAMS.iNumCorrel - 1, 0,
                       STATIC_PARAMS.iNumCorrel - 1);
          ppdCorrMatrixToTest = NULL;
        }

        are_definite_positive = SRT_NO;
        return are_definite_positive;
      }
    }
  }

  if (pdEigenval)
    free_dvector(pdEigenval, 0, STATIC_PARAMS.iNumCorrel - 1);
  pdEigenval = NULL;

  if (ppdEigenvector) {
    free_dmatrix(ppdEigenvector, 0, STATIC_PARAMS.iNumCorrel - 1, 0,
                 STATIC_PARAMS.iNumCorrel - 1);
    ppdEigenvector = NULL;
  }

  if (ppdCorrMatrixToTest) {
    free_dmatrix(ppdCorrMatrixToTest, 0, STATIC_PARAMS.iNumCorrel - 1, 0,
                 STATIC_PARAMS.iNumCorrel - 1);
    ppdCorrMatrixToTest = NULL;
  }

  return are_definite_positive;

} /* End of are_correl_matrix_definite_positive */

/* This function tests if the optim params are FX Params  */

static SRT_Boolean are_optparams_fxparams(double *pdOptimParams) {

  double **ppdFXVolCrv = NULL;
  double ***pppdCorrTensor = NULL;
  SRT_Boolean are_fxparams = SRT_YES;
  Err err = NULL;

  /* Get the VolCurve from OptParams */

  ppdFXVolCrv = dmatrix(0, 1, 0, STATIC_PARAMS.iNumLocVols - 1);
  from_optparam_to_fx_vol_curve(pdOptimParams, ppdFXVolCrv);

  if (are_fx_vol_positive(ppdFXVolCrv) == SRT_NO)
    are_fxparams = SRT_NO;

  /* Get the Correl from OptParams */

  if (STATIC_PARAMS.bCalibCorrel == SRT_YES) {
    pppdCorrTensor = f3tensor(0, STATIC_PARAMS.iNumCorrel - 1, 0,
                              STATIC_PARAMS.iNumCorrel - 1, 0,
                              STATIC_PARAMS.iNumCorrDates - 1);

    from_optparam_to_correl(pdOptimParams, pppdCorrTensor);
    if (are_correl_matrix_definite_positive(pppdCorrTensor) == SRT_NO)
      are_fxparams = SRT_NO;
  }

  if (ppdFXVolCrv)
    free_dmatrix(ppdFXVolCrv, 0, 1, 0, STATIC_PARAMS.iNumLocVols - 1);
  if (pppdCorrTensor)
    free_f3tensor(pppdCorrTensor, 0, STATIC_PARAMS.iNumCorrel - 1, 0,
                  STATIC_PARAMS.iNumCorrel - 1, 0,
                  STATIC_PARAMS.iNumCorrDates - 1);

  return are_fxparams;

} /* End of are_optparams_fxparams */

/* This function get the correlation for the initialistaion of the IR from the
        correlation tensor */

static void get_correlation_from_correl_matrix(double ***pppdCorrTensor,
                                               double **ppdCorrelation) {
  int i, j, k;

  for (k = 0; k < STATIC_PARAMS.iNumCorrDates; k++) {
    for (i = 0; i < STATIC_PARAMS.iNumCorrel; i++) {
      for (j = 0; j < i; j++)
        ppdCorrelation[i + j - 1][k] = pppdCorrTensor[i][j][k];
    }
  } /* i=0 domund        ,i=1 forund        , i=2 fxund */

} /* End of get_correlation_from_correl_matrix */

static Err from_volcorrel_to_optparam(double **ppdFXVolCurve,
                                      double *pdOptimParams) {

  int i, j;
  int iIndex;

  iIndex = 0;
  if (STATIC_PARAMS.bCalibCorrel == SRT_YES) {
    for (j = 0; j < STATIC_PARAMS.iNumCorrDates; j++) {
      for (i = 0; i < STATIC_PARAMS.iNumCorrel; i++) {
        iIndex = j * STATIC_PARAMS.iNumCorrel + (i + 1);
        pdOptimParams[iIndex] = STATIC_PARAMS.ppdInitCorrels[i][j];
      }
    }

    for (i = (iIndex + 1); i < (STATIC_PARAMS.iNumParams + 1); i++) {
      j = i - STATIC_PARAMS.iNumCorrDates * STATIC_PARAMS.iNumCorrel - 1;
      pdOptimParams[i] = ppdFXVolCurve[1][j];
    }
  } else {
    for (i = 1; i < (STATIC_PARAMS.iNumParams + 1); i++) {
      pdOptimParams[i] = ppdFXVolCurve[1][i - 1];
    }
  }

  return NULL;
}

static Err from_optparam_to_fx_und(double *pdOptimParams) {

  Err err = NULL;
  SrtUndPtr undptr = NULL;
  SrtUndPtr dom_und, for_und;
  String dom_ccy, for_ccy;
  SrtMdlType dom_mdl_type;
  Date today;
  TermStruct *ts;
  SrtUndListPtr und_list;

  double ***pppdCorrTensor = NULL;
  double **ppdCorrelation = NULL;
  double **ppdFXVolCrv = NULL;

  if (are_optparams_fxparams(pdOptimParams) == SRT_NO)
    return err;

  if (STATIC_PARAMS.bCalibCorrel == SRT_YES) {

    pppdCorrTensor = f3tensor(0, STATIC_PARAMS.iNumCorrel - 1, 0,
                              STATIC_PARAMS.iNumCorrel - 1, 0,
                              STATIC_PARAMS.iNumCorrDates - 1);

    if (!pppdCorrTensor)
      return serror("Allocation Failure ");

    from_optparam_to_correl(pdOptimParams, pppdCorrTensor);

    ppdCorrelation = dmatrix(0, STATIC_PARAMS.iNumCorrel - 1, 0,
                             STATIC_PARAMS.iNumCorrDates - 1);
    if (!ppdCorrelation)
      return serror("Allocation Failure ");

    get_correlation_from_correl_matrix(pppdCorrTensor, ppdCorrelation);

    /* Replaces the current correlation matrix by the new one (with coeffcients
     * done inside) in UndInfo */
    err = SrtInitCorrelationMatrix(
        STATIC_PARAMS.iNumCorrDates, STATIC_PARAMS.iNumCorrel, ppdCorrelation,
        STATIC_PARAMS.pdCorrDates, STATIC_PARAMS.ppszUndNames);
    if (err)
      return err;
  }

  ppdFXVolCrv = dmatrix(0, 1, 0, STATIC_PARAMS.iNumLocVols - 1);
  from_optparam_to_fx_vol_curve(pdOptimParams, ppdFXVolCrv);

  dom_und = lookup_und(STATIC_PARAMS.szDomUndName);
  dom_ccy = get_underlying_ccy(dom_und);

  err = get_underlying_mdltype(dom_und, &dom_mdl_type);

  today = get_today_from_underlying(dom_und);

  for_und = lookup_und(STATIC_PARAMS.szForUndName);
  for_ccy = get_underlying_ccy(for_und);

  err = get_underlying_ccy(for_und);

  srt_f_destroy_und(STATIC_PARAMS.szFXUndName);

  err = srt_f_init_FX_TermStruct(
      today, ppdFXVolCrv, 2, STATIC_PARAMS.iNumLocVols,
      STATIC_PARAMS.eModelType, 0.0, 0.0, STATIC_PARAMS.szFXUndName,
      STATIC_PARAMS.szDomUndName, STATIC_PARAMS.szForUndName, &ts);
  if (err)
    return err;

  /*Get the underlying list and check it has been initialised */
  und_list = get_underlying_list();
  if (und_list == NULL)
    return serror("NO Underlying list");

  /* Puts the underlying in the Market List */

  err = srt_f_addundtolist(
      und_list, STATIC_PARAMS.szFXUndName, "FX_UND", dom_ccy, "FX_STOCH_RATES",
      STATIC_PARAMS.szDomUndName, STATIC_PARAMS.szForUndName, NULL, ts,
      STATIC_PARAMS.dFXSpot);
  if (err)
    return err;

  if (ppdFXVolCrv) {
    free_dmatrix(ppdFXVolCrv, 0, 1, 0, STATIC_PARAMS.iNumLocVols - 1);
    ppdFXVolCrv = NULL;
  }

  if (ppdCorrelation) {
    free_dmatrix(ppdCorrelation, 0, STATIC_PARAMS.iNumCorrel - 1, 0,
                 STATIC_PARAMS.iNumCorrDates - 1);
    ppdCorrelation = NULL;
  }

  if (pppdCorrTensor) {
    free_f3tensor(pppdCorrTensor, 0, STATIC_PARAMS.iNumCorrel - 1, 0,
                  STATIC_PARAMS.iNumCorrel - 1, 0,
                  STATIC_PARAMS.iNumCorrDates - 1);
    pppdCorrTensor = NULL;
  }

  return err;
}

/* This Function computes ONE BS Implied Volatility of an FX Option
         for a given set of optimisation parameters.
        This function is called by levenberg_calib_funcs
*/

static Err levenberg_fxvol_funcs(int iInstrIndex, double *pdOptimParams,
                                 double *pdFXBSVol)

{

  Err err = NULL;
  double dOptMat;

  if (are_optparams_fxparams(pdOptimParams) == SRT_NO) {
    *pdFXBSVol = SRT_BIG;
    return NULL;
  }

  err = from_optparam_to_fx_und(pdOptimParams);
  if (err)
    return err;

  dOptMat = STATIC_INSTRUMENT.pdFXOptMats[iInstrIndex - 1];
  err = srt_f_get_fx_implied_vol(dOptMat, STATIC_PARAMS.szFXUndName, pdFXBSVol);

  /*	(*pdFXBSVol)*=100; */

  if (*pdFXBSVol > DBL_MAX)
    serror("Calibration failed: get a FXVOL overflow: %f", *pdFXBSVol);

  return err;

} /*End of levenberg_fxvol_funcs */

static Err levenberg_calib_funcs(double iInstrIndex, double *pdOptimParams,
                                 double *value, double deriv[],
                                 int iNumParams) {
  long i;
  double shift;
  Err err = NULL;

  err = levenberg_fxvol_funcs((int)iInstrIndex, pdOptimParams, value);
  if (err)
    return err;

  /* Computes the derivatives of this value with respect to each parameters */
  for (i = 1; i <= iNumParams; i++) {
    if (pdOptimParams[i] == 0)
      shift = FX_ZERO_SHIFT;
    else
      shift = FX_PROP_SHIFT * pdOptimParams[i];
    pdOptimParams[i] += shift;

    /* Computes the shifted pdMktVol: stores it in deriv */

    err = levenberg_fxvol_funcs((int)iInstrIndex, pdOptimParams, &(deriv[i]));
    if (err)
      return err;

    deriv[i] -= (*value);
    deriv[i] /= shift;

    /*resets params to its initial value */

    pdOptimParams[i] -= shift;
  }

  /* Return a success string */
  return NULL;

} /* End of levenberg_calib_funcs */

static void free_calib_fx_core_memory(int iNumInstrs, int iNumLocVols,
                                      int iNumParams) {

  if (STATIC_INSTRUMENT.pdFXOptMats) {
    free_dvector(STATIC_INSTRUMENT.pdFXOptMats, 0, iNumInstrs - 1);
    STATIC_INSTRUMENT.pdFXOptMats = NULL;
  }

  if (STATIC_INSTRUMENT.pdFXOptMktVols) {
    free_dvector(STATIC_INSTRUMENT.pdFXOptMktVols, 1, iNumInstrs);
    STATIC_INSTRUMENT.pdFXOptMktVols = NULL;
  }

  if (STATIC_INSTRUMENT.pdWeights) {
    free_dvector(STATIC_INSTRUMENT.pdWeights, 1, iNumInstrs);
    STATIC_INSTRUMENT.pdWeights = NULL;
  }

  if (STATIC_PARAMS.pdCorrDates) {
    free_dvector(STATIC_PARAMS.pdCorrDates, 0, STATIC_PARAMS.iNumCorrDates - 1);
    STATIC_PARAMS.pdCorrDates = NULL;
  }

  if (STATIC_INSTRUMENT.plFXOptMDates) {
    free_lvector(STATIC_INSTRUMENT.plFXOptMDates, 0, iNumInstrs - 1);
    STATIC_INSTRUMENT.plFXOptMDates = NULL;
  }

  if (STATIC_PARAMS.plFXVolCrvDates) {
    free_lngvector(STATIC_PARAMS.plFXVolCrvDates, 0, iNumLocVols - 1);
    STATIC_PARAMS.plFXVolCrvDates = NULL;
  }

  if (STATIC_PARAMS.ppdInitCorrels) {
    free_dmatrix(STATIC_PARAMS.ppdInitCorrels, 0, STATIC_PARAMS.iNumCorrel - 1,
                 0, STATIC_PARAMS.iNumCorrDates - 1);
    STATIC_PARAMS.ppdInitCorrels = NULL;
  }
}

Err srt_f_calib_fx_core(char **pszFXCalibStrings, char **pszFXCalibValueStrings,
                        int iNumCalibParamss,

                        int iNumParams,

                        int iNumLocVols, double **ppdFXVolCrv,

                        int iNumInstrs, long *plFXOptMDates,
                        double *pdFXOptMktVols, /*From 1 to iNumInstrs */
                        double *pdWeights,      /*From 1 to iNumInstrs */

                        char *szFXUndName, double *ChiSquare)

{
  Err err = NULL;
  SrtUndPtr sFXUnd;
  SrtFXCalibParam sFXCalibParams;
  double *pdOptimParams = NULL;
  double *pdData = NULL;
  int i;

  err = srt_f_set_FXCalibParams(pszFXCalibStrings, pszFXCalibValueStrings,
                                iNumCalibParamss, &sFXCalibParams);
  if (err)
    return err;

  STATIC_PARAMS.bCalibCorrel = sFXCalibParams.bCalibCorrel;

  if (!(STATIC_PARAMS.bCalibCorrel == SRT_YES))
    iNumParams = iNumLocVols;

  err = set_static_structures_for_fx_calib(
      szFXUndName, plFXOptMDates, pdFXOptMktVols, pdWeights, iNumParams,
      iNumLocVols, ppdFXVolCrv, iNumInstrs);

  if (err)
    return err;

  /*-------------------------INITIALISATION-----------------------------------*/

  pdOptimParams = dvector(1, iNumParams);

  err = from_volcorrel_to_optparam(ppdFXVolCrv, pdOptimParams);

  /*----------------------- Calibration Launch------------------------------ */

  /*Safety measure: free the TermStructure attached to the underlying */

  sFXUnd = lookup_und(szFXUndName);

  if (sFXUnd == NULL)
    return serror("Undefined Underlying %s", szFXUndName);

  err = free_underlying_ts(sFXUnd);

  if (sFXCalibParams.lNumIter > 0) {
    if (sFXCalibParams.eAlgoType == LEVENBERG_MARQUARDT) {
      smessage("Calibration using Levenberg-Marquardt algorithm");
      smessage("");

      /* Sets a vector of x data: here it is 1        ,2        ,3....: the
       * market pdData index */

      iNumParams = STATIC_PARAMS.iNumParams;
      pdData = dvector(1, iNumParams);

      for (i = 1; i <= iNumParams; i++)
        pdData[i] = (double)(i);
      /* Call the Levenberg-Marquardt main routine with the right function */

      err = levenberg_marquardt(pdData, pdFXOptMktVols, pdWeights, iNumInstrs,
                                pdOptimParams, iNumParams,
                                sFXCalibParams.lNumIter, levenberg_calib_funcs,
                                ChiSquare);

      if (err) {
        free_calib_fx_core_memory(iNumInstrs, iNumLocVols, iNumParams);

        if (pdOptimParams)
          free_dvector(pdOptimParams, 1, iNumParams);
        pdOptimParams = NULL;

        if (pdData)
          free_dvector(pdData, 1, iNumParams);
        pdData = NULL;
      }
    }
  }

  free_calib_fx_core_memory(iNumInstrs, iNumLocVols, iNumParams);

  if (pdOptimParams)
    free_dvector(pdOptimParams, 1, iNumParams);
  pdOptimParams = NULL;

  if (pdData)
    free_dvector(pdData, 1, iNumParams);
  pdData = NULL;

  return err;
}
