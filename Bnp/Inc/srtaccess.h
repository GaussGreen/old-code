/**********************************************************************
 *      Name: SrtUndAccess.h                                          *
 *  Function: Header for functions adding objects to SrtMkt structure *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 12/10/95                                                *
 *--------------------------------------------------------------------*
 *    Inputs:                                                         *
 *   Returns:                                                         *
 *   Globals:                                                         *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 12/10/95 FOS     Created for SORT5-GRFN3 port to NT for SPG        *
 **********************************************************************/

#ifndef SRTACCESS_H
#define SRTACCESS_H

#include "SrtAccessLon.h"
#include "SrtAccessTyo.h"
#include "grf_h_all.h"
#include "srt_h_all.h"

#ifdef __cplusplus
extern "C"
{
#endif

    char* SrtGrfnCHEYBETA(
        char*    underlying,
        int      numeventdates,
        long*    eventdates,
        long     tableauRows,
        long     tableauCols,
        char***  tableauStrings,
        int**    tableauMask,
        long     auxWidth,
        long*    auxLen,
        double** aux,
        /* param */
        int method, /* 0 MC, 1 MC adj, 2 PDE */
        int nstept,
        /* for MC */
        long         numpaths,
        SrtMCSamType gen_method,
        /* for PDE */
        int      nstepx,
        int      nstepphi,
        int*     nb_prod,
        double** prod_val);

    char* SrtGrfnCHEYBETApde(
        char*    underlying,
        int      numeventdates,
        long*    eventdates,
        long     tableauRows,
        long     tableauCols,
        char***  tableauStrings,
        int**    tableauMask,
        long     auxWidth,
        long*    auxLen,
        double** aux,
        int      nstept,
        int      nstepx,
        int      nstepphi,
        int*     nb_prod,
        double** prod_val);

    char* SrtGrfnMain(
        char*    domestic,
        int      numParams,
        char**   paramStrings,
        char**   valueStrings,
        int      numeventdates,
        long*    eventdates,
        long     tableauRows,
        long     tableauCols,
        char***  tableauStrings,
        int**    tableauMask,
        long     auxWidth,
        long*    auxLen,
        double** aux,
        double*  price,
        double*  stdev,
        double** grfn_cells,
        double** pay_report);

    char* SrtGrfnMainExFrontier(
        char*    domestic,
        int      numParams,
        char**   paramStrings,
        char**   valueStrings,
        int      numeventdates,
        long*    eventdates,
        long     tableauRows,
        long     tableauCols,
        char***  tableauStrings,
        int**    tableauMask,
        long     auxWidth,
        long*    auxLen,
        double** aux,
        double*  price,
        double*  stdev,
        double*  exfrontier,
        double** grfn_cells,
        double** pay_report);

    char* SrtGrfAmer(
        int     numParams,
        char**  paramStrings,
        char**  valueStrings,
        char*   undName,
        double  start,
        double  end,
        int     nfp,
        int     delay,
        char*   cpdStr,
        char*   basisStr,
        double* strikes,
        int     nStrikes,
        char*   recPayStr,
        double* price);

    char* SrtGrfMidat(
        int     numParams,
        char**  paramStrings,
        char**  valueStrings,
        char*   undName,
        int     num_exercise_dates,
        long*   exercise_dates,
        long*   exercise_start_dates,
        double* exercise_premiums,
        int     num_prod_dates,
        long*   prod_dates,
        double* prod_cfs,
        char*   recPayStr,
        double* price);

    char* SrtClosedForm(
        long     startLong,
        long     endLong,
        char*    compound,
        char*    basis,
        char*    recPay,
        char*    dealType,
        double   strike,
        double   bondStrike,
        char*    refRateCodeStr,
        int      mdlRows,
        char**   paramStrings,
        char**   valueStrings,
        char*    undName,
        double*  price,
        double** sigma_vega,
        long*    num_sig,
        double** tau_vega,
        long*    num_tau,
        char*    greekStr);

    Err SrtFutureClosedForm(
        long    lFutureDate,
        long    lStartDate,
        long    lNfpOrEnd,
        char*   szComp,
        char*   szBasis,
        char*   szRecPay,
        char*   szType,
        double  dStrike,
        double  dBondStrike,
        char*   szRefRateCode,
        int     iNumParams,
        char**  pszParamStrings,
        char**  pszValueStrings,
        char*   szUndName,
        double* pdFuturePrice);

    Err SrtNewCalibrateAll(
        String* pszGrfnParamNames,
        String* pszGrfnParamValues,
        long    lNumGrfnParams,
        long*   plStartDates,
        long*   plEndDatesOrNfp,
        String* pszFrequency,
        String* pszBasis,
        double* pdStrike,
        double* pdBondStrike,
        String* pszOptionType,
        String* pszRecPay,
        String* pszRefRateCode,
        double* pdPrice,
        double* pdVega,
        Err (*pfGetVol)(
            double  dStart,
            double  dEnd,
            double  dStrike,
            double  dForward,
            double  dSpread,
            double* pdBsVol),
        String*   LogNormStr,
        long      lNumInstruments,
        String*   pszFraTenors,
        double**  pdCorrelationMatrix,
        long      lNumTenors,
        String    szUndName,
        String    szNewName,
        String*   pszCalibParamNames,
        String*   pszCalibParamValues,
        long      lNumCalibParams,
        double*** pppdSigmaCurve,
        long*     plNumSigmas,
        long*     plNumSigmaCols,
        double*** pppdTauCurve,
        long*     plNumTaus,
        long*     plNumTauCols,
        double*   pdChiSquare,
        double**  ppdTheoPrices,
        int*      indexUsedInstr);

    char* WrapSrtBootstrap(
        int     numParams,
        char**  paramStrings,
        char**  valueStrings,
        int*    nInputs,
        long*   startDates,
        long*   endDates,
        char**  cpdStr,
        char**  basisStr,
        double* strike,
        double* bondStrike,
        char**  typeStr,
        char**  recPayStr,
        String* refRateStr,
        double* undPrice,
        Err (*pfGetVol)(
            double  dStart,
            double  dEnd,
            double  dStrike,
            double  dForward,
            double  dSpread,
            double* pdBsVol),
        String* LogNormStr,
        double  tau,
        double  alpha,
        double  beta,
        double  rho,
        char*   undName);

    char* SrtBootstrap(
        int     numParams,
        char**  paramStrings,
        char**  valueStrings,
        int*    nInputs,
        long*   startDates,
        long*   endDates,
        char**  cpdStr,
        char**  basisStr,
        double* strike,
        double* bondStrike,
        char**  typeStr,
        char**  recPayStr,
        String* refRateStr,
        double* undPrice,
        double* ATMprice,
        double  tau,
        double  alpha,
        double  beta,
        double  rho,
        char*   undName);

    Err srt_f_get_fx_implied_vol(double yr_to_exp, String fx_und_names, double* fx_bs_vol);

    char* SrtVersion();

    int SrtInit();

    int SrtClose();

    char* SrtInitIRUnd(
        char*    undName,
        char*    ycName,
        char*    model,
        int      volCrvRows,
        int      volCrvCols,
        double** volCrvVals,
        int      tauCrvRows,
        int      tauCrvCols,
        double** tauCrvVals,
        double   beta,
        double   alpha,
        double   gamma,
        double   rho,
        double   vovol,
        double   etaOrBeta2,
        double   vasicek_init_cond,
        int      vasicek_mean_rev_level_rows,
        int      vasicek_mean_rev_level_cols,
        double** vasicek_mean_rev_level_vals);

    char* SrtInitFXUnd(
        char*    undName,
        double   spot,
        char*    model,
        char*    domDiscName,
        char*    forDiscName,
        int      volCrvRows,
        int      volCrvCols,
        double** volCrvVals,
        /* OUTPUT */
        char* definite_name);

    char* SrtInitEQUnd(
        char*    undName,
        char*    model,
        double   spot,
        char*    ycName1,
        char*    ycName2,
        char*    repoName,
        int      volCrvRows,
        int      volCrvCols,
        double** volCrvVals,
        double   omega,
        double   beta,
        double   gamma,
        double   voldrift,
        double   vovol,
        double   rho);

    char* SrtInitDividends(
        long today, char* fwdName, char* ccy, int nRows, int nCols, double** fVals);

    char* SrtInitRepos(
        long today, char* repo_name, char* ccy, int nrows, int ncols, double** repo_vals);

    Err SrtInitCorrelationMatrix(
        int ndates, int ncorr, double** correl, double* dates, String** und_names);

#ifdef __cplusplus
}
#endif

#endif