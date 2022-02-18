
#ifndef MCEBOPTIMISATIONH
#define MCEBOPTIMISATIONH

#include "grf_h_all.h"

typedef struct
{
    /* GRFN Read Params */
    int iColPay;
    int iColBound;

    /* MCEB Settings */
    long lNbDates;
    int  iCallCurrent;
    int  iIsKO;
    int  iMultiIndex;
    int  iNbIndex;

    int    iDoSmoothing;
    double dCallSpread;

    /* Optimisation Outputs */
    double*  dBarrier;
    double** dCoefLin;

    /* Fees for Call */
    int iHasFees;
    int iColFees;
    int iPayAllFeeCol;

    /* Fwd IV calculation / adjustment */
    int     iCalcIV;
    int     iAdjustIV;
    int     iHasNumeraire;
    int     iColNumeraire;
    int     iAddMultFee;
    double* dMarketFwdIV;
    double* dModelFwdIV;
    double* dFee;

    /* All the extra informations */
    int iCalcOneTime;
    int iCalcOneTimePartial;
    int iCalcExeProba;

    double* dOneTimeCall;
    double* dOneTimePartial;
    double* dExeProba;

    /* Extra Parameters */
    int iDoInfos;
    int iRemoveLastOnLast;
    int iRescaleLinCoefs;
    int iKnockInCol;
    int iFindBestOptim;
    int iAddNonOptimisedForKI;

    /* File Name where to save all infos */
    char sFilePath[50];

} MCEBParams, *MCEBPARAMS;

void mceb_set_default_params(MCEBPARAMS sParams);

void mceb_copy_nondynamic_params(MCEBPARAMS sParamsCopy, MCEBPARAMS sParamsSrc);

Err mceb_allocate_params(MCEBPARAMS sParams, long lNbDates);

void mceb_free_params(MCEBPARAMS sParams);

void mceb_shift_extrainfos(MCEBPARAMS sParams);

Err mceb_adjust_fwdiv(
    double***  save_values,
    long       start_date_idx,
    long       end_date_idx,
    long       nb_paths,
    int*       optimise,
    MCEBPARAMS params);

Err mceb_allocate_savevalues_for_GRFN(
    long iNbPaths, long iNbEvent, MCEBPARAMS sParams, double**** dSaveValues);

void mceb_free_savevalues_for_GRFN(
    double*** dSaveValues, long iNbPaths, long iNbEvent, MCEBPARAMS sParams);

void mceb_fill_savevalues_from_GRFN(
    double** dSaveValues, double* dResEvt, long lSimulIdx, double dDF, MCEBPARAMS sParams);

Err find_and_optimise_boundary(
    double***  save_values,
    long       nb_dates,
    long       nb_paths,
    int*       optimise,
    MCEBPARAMS params,
    double*    value,
    double*    error);

Err find_linear_coef(
    long     p, /* number of variables = number of row in Transp(X) */
    long     n, /* number of observations */
    double*  y,
    double** x,
    double*  b,
    double** mat /* matrix p row, n col */);

Err find_best_optim_dates(
    double***  save_values,
    long       nb_dates,
    long       nb_target,
    long       nb_paths,
    int*       optimise,
    MCEBPARAMS params,
    double*    value,
    double*    error);

Err optimise_boundary_info(
    double*** save_values,
    long      nb_dates,
    long      nb_paths,
    int*      optimise,
    int       call_current,
    int       is_ko,
    double*   boundary,
    double*   value,
    double*   error,
    double**  infos); /* Col 1: One time Call, Col2: Exerc Proba, Col3: New PV */

#endif
