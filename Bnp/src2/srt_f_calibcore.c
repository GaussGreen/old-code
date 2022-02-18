/* ========================================================================

  FILENAME:  srt_f_calibcore.c


      PURPOSE:   The core for any interest rate calibration.

                FUNCTION:  Err srt_f_calib_core(...)

========================================================================== */

/* Include files */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_closedform.h"
#include "srt_h_grfclsdfrm.h"
#include "swp_h_all.h"

/* global pointer to the array of computed calibrated model prices
   modified in levenberg_calib_funcs */
extern double* GlobalTheoPrices;

#ifdef PVMPI
#include "parallel.h"
#endif

#ifdef PVMPI
#define ID_SET_INSTRMTS 215
#define ID_SET_PARAMS 220
#define ID_FREE_MEMORY 230
#endif

/* A numerically big number to prevent a search in a wrong direction */
#define SRT_BIG 1.0e+20

/* Shifts to compute derivatives of pdMktPrice wrt parameters */
#define ZERO_SHIFT 1.0e-06
#define PROP_SHIFT 1.0e-03

/* ----------------------------------------------------------------------- */
#ifdef PVMPI
#ifdef _DEBUG
#pragma message("warning: remove trace stuff")
void MyTRACE(char*);
void MyTRACE1(char*, double d);
void MyTRACE1n(char*, double d);
void MyTRACE2(char*, double d1, double d2);
void MyTRACE3(char*, double d1, double d2, double d3);
#endif
#endif

/* -----------------------------------------------------------------------
PART 1
STATIC STRUCTURES DECLARATION
----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
A static structure to store parameters linked to the optimisation
routine used (number of parameters, ...)
----------------------------------------------------------------------- */
#ifndef PVMPI
static struct
{
    long     lNumData;
    double*  pdMarketTargets; /* From [1] to [lNumData] */
    double*  pdMarketWeights; /* From [1] to [lNumData] */
    long     lNumParams;
    double** ppdParamBounds; /* [1] is minimum ; [2] is maximum */

    double** ppdSigmaValues;
    long     lNumSigmas;
    long     lNumSigmaCols;
    double** ppdTauValues;
    long     lNumTaus;
    long     lNumTauCols;

    SRT_Boolean bFreezeTau;
    SRT_Boolean bOneTau;
    double*     dFixedTau;
    SRT_Boolean bFreezeBeta;
    SRT_Boolean bOneBeta;
    double      dFixedBeta;
    SRT_Boolean bFreezeOmega;
    SRT_Boolean bOneOmega;
    double      dFixedOmega;

    double dFixedAlpha;
    double dFixedGamma;
    double dFixedRho;

    SrtCalibType eCalibType;

    SrtMdlDim    eModelDim;
    SrtModelType eModelType;
    SRT_Boolean  bSmoothSigma; /* use smooth sigma criterion */

    char szUndName[64];
} STATIC_PARAMS;
#else
/* If I meet the "£$%^&*I( bloke who programmed this calibration...*/
PARAMS_STRUCT      STATIC_PARAMS;
#endif

/* -----------------------------------------------------------------------
A static structure to store market parameters neede in the prcing
functions used during the optimisation
----------------------------------------------------------------------- */

#ifndef PVMPI
static struct
{
    long             lClcnDate;
    long             lNumInstruments;
    SwapDP*          psSwapDp;
    double*          pdStrike;
    double*          pdBondStrike;
    SrtReceiverType* peRecPay;
    StructType*      peProductType;
    String*          pszRefRateCode;
    long             lNumTenors;
    double*          pdFraMaturities;
    SrtGrfnParam*    psGrfnParams;
} STATIC_INSTRUMENT;
#else
INSTRUMENTS_STRUCT STATIC_INSTRUMENT;
#endif
/* -----------------------------------------------------------------------
A static structure to store spreadsheet neede in the prcing
functions used during the optimisation
----------------------------------------------------------------------- */
static struct
{
    GrfnCell** pshtSHEET;
    Date*      pszEventDates;
    long       lNumEventDates;
    long       lNumRows;
    long       lNumCols;
} STATIC_SPRDSHT;

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
A static structure to store model price
functions used during the optimisation
----------------------------------------------------------------------- */
#ifndef PVMPI
static struct
{
#else
struct
{
#endif
    double*  pdModelValue;
    double** pdModelHessian;
} STATIC_MODELVALUES;

/* ----------------------------------------------------------------------- */
#ifdef PVMPI
struct prl_context PRL_CONTEXT;
#endif
/* -----------------------------------------------------------------------
PART 2
STATIC FUNCTIONS TO FILL IN STATIC STRUCTURES
----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
Stores all the relevant market information required for pricing during
the optimisation process
All the arrays here start at [0] and go up to [n-1]
----------------------------------------------------------------------- */

#ifdef PVMPI
Err set_static_instr_for_calib
#else
static Err set_static_instr_for_calib
#endif
    (long             lClcnDate,
     long             lNumInstruments,
     SwapDP*          psSwapDp,
     double*          pdStrike,
     double*          pdBondStrike,
     SrtReceiverType* peRecPay,
     StructType*      peProductType,
     String*          pszRefRateCode,
     long             lNumTenors,
     double*          pdFraMaturities,
     SrtGrfnParam*    psGrfnParams)
{
#ifdef PVMPI
    if (!PRL_CONTEXT.bLocalComp)
        SendInstruments(
            lClcnDate,
            lNumInstruments,
            psSwapDp,
            pdStrike,
            pdBondStrike,
            peRecPay,
            peProductType,
            pszRefRateCode,
            lNumTenors,
            pdFraMaturities,
            psGrfnParams);
#endif
    STATIC_INSTRUMENT.lClcnDate       = lClcnDate;
    STATIC_INSTRUMENT.lNumInstruments = lNumInstruments;
    STATIC_INSTRUMENT.psSwapDp        = psSwapDp;
    STATIC_INSTRUMENT.pdStrike        = pdStrike;
    STATIC_INSTRUMENT.pdBondStrike    = pdBondStrike;
    STATIC_INSTRUMENT.peRecPay        = peRecPay;
    STATIC_INSTRUMENT.peProductType   = peProductType;
    STATIC_INSTRUMENT.pszRefRateCode  = pszRefRateCode;
    STATIC_INSTRUMENT.lNumTenors      = lNumTenors;
    STATIC_INSTRUMENT.pdFraMaturities = pdFraMaturities;
    STATIC_INSTRUMENT.psGrfnParams    = psGrfnParams;

    return NULL;
}

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
Build the general spreadsheet neede in the prcing
functions used during the optimisation and the FREE function associated
----------------------------------------------------------------------- */

#ifdef PVMPI
Err set_static_sprdsht_for_calib
#else
static Err set_static_sprdsht_for_calib
#endif
    (String und_name, long today)
{
    Err err = NULL;
#ifdef PVMPI
    if (!PRL_CONTEXT.bLocalComp)
        SendSprdsht(und_name, today);
#endif
    err = new_grfn_SwapDParray_to_GrfnCells(
        &STATIC_SPRDSHT.lNumEventDates,
        &STATIC_SPRDSHT.pszEventDates,
        &STATIC_SPRDSHT.lNumRows,
        &STATIC_SPRDSHT.lNumCols,
        &STATIC_SPRDSHT.pshtSHEET,
        today,
        STATIC_INSTRUMENT.lNumInstruments,
        STATIC_INSTRUMENT.psSwapDp,
        STATIC_INSTRUMENT.pdStrike,
        STATIC_INSTRUMENT.pdBondStrike,
        STATIC_INSTRUMENT.peRecPay,
        STATIC_INSTRUMENT.peProductType,
        und_name,
        STATIC_INSTRUMENT.pszRefRateCode);
#ifdef PVMPI
    if (!PRL_CONTEXT.bLocalComp)
    {
#endif
        STATIC_MODELVALUES.pdModelValue = dvector(0, STATIC_PARAMS.lNumData - 1);
        STATIC_MODELVALUES.pdModelHessian =
            dmatrix(0, STATIC_PARAMS.lNumData - 1, 0, STATIC_PARAMS.lNumParams - 1);
#ifdef PVMPI
    }
#endif
    return err;
}

#ifdef PVMPI
Err free_space_for_static_sprdsht_for_calib
#else
static Err free_space_for_static_sprdsht_for_calib
#endif
    (long lNumRows, long lNumCols)
{
    if (STATIC_SPRDSHT.pszEventDates)
        srt_free(STATIC_SPRDSHT.pszEventDates);

    if (STATIC_SPRDSHT.pshtSHEET)
        grfn_free_GrfnCellmatrix(
            STATIC_SPRDSHT.pshtSHEET, STATIC_SPRDSHT.lNumRows, STATIC_SPRDSHT.lNumCols);

#ifdef PVMPI
    if (!PRL_CONTEXT.bLocalComp)
    {
#endif
        if (STATIC_MODELVALUES.pdModelValue)
            free_dvector(STATIC_MODELVALUES.pdModelValue, 0, STATIC_PARAMS.lNumData - 1);
        if (STATIC_MODELVALUES.pdModelHessian)
            free_dmatrix(
                STATIC_MODELVALUES.pdModelHessian,
                0,
                STATIC_PARAMS.lNumData - 1,
                0,
                STATIC_PARAMS.lNumParams - 1);
        STATIC_MODELVALUES.pdModelValue   = NULL;
        STATIC_MODELVALUES.pdModelHessian = NULL;
#ifdef PVMPI
    }
#endif
    return NULL;
}

/* -------------------------------------------------------------------------*/

/* -----------------------------------------------------------------------
Stores a few useful indexes and pdValues (mkt_price, pdMarketWeightss...) in a
static structure for calibration purposes.
----------------------------------------------------------------------- */

#ifdef PVMPI
Err set_static_params_for_calib
#else
static Err set_static_params_for_calib
#endif
    (SrtCalibType eCalibType,
     long         lNumData,
     double*      pdMarketWeights,
     double*      pdMarketTargets,
     long         lNumParams,
     long         lNumSigmas,
     long         lNumTaus,
     double**     ppdParamBounds,
     SRT_Boolean  bFreezeTau,
     SRT_Boolean  bOneTau,
     double**     dFixedTau,
     SRT_Boolean  bFreezeBeta,
     SRT_Boolean  bOneBeta,
     double       dFixedBeta,
     SRT_Boolean  bFreezeOmega,
     SRT_Boolean  bOneOmega,
     double       dFixedOmega,
     double       dFixedAlpha,
     double       dFixedGamma,
     double       dFixedRho,
     SrtModelType eModelType,
     SrtMdlDim    eModelDim,
     SRT_Boolean  bSmoothSigma,
     String       szUndName)
{
#ifdef PVMPI
    if (!PRL_CONTEXT.bLocalComp)
        SendParams(
            eCalibType,
            lNumData,
            pdMarketWeights,
            pdMarketTargets,
            lNumParams,
            lNumSigmas,
            lNumTaus,
            ppdParamBounds,
            bFreezeTau,
            bOneTau,
            dFixedTau,
            bFreezeBeta,
            bOneBeta,
            dFixedBeta,
            bFreezeOmega,
            bOneOmega,
            dFixedOmega,
            dFixedAlpha,
            dFixedGamma,
            dFixedRho,
            eModelType,
            eModelDim,
            bSmoothSigma,
            szUndName);
#endif /*PVMPI*/

    /* Calib peProductType */
    STATIC_PARAMS.eCalibType = eCalibType;

    /* How many parameters do we use for calibration */
    STATIC_PARAMS.lNumParams = lNumParams;
    /* How many market pdData are there  */
    STATIC_PARAMS.lNumData = lNumData;

    /* If calibration is done with a fixed tau, set its value */
    STATIC_PARAMS.dFixedTau  = (*dFixedTau);
    STATIC_PARAMS.bFreezeTau = bFreezeTau;
    STATIC_PARAMS.bOneTau    = bOneTau;

    /* If calibration is done with a fixed beta, set its value */
    STATIC_PARAMS.dFixedBeta  = dFixedBeta;
    STATIC_PARAMS.bFreezeBeta = bFreezeBeta;
    STATIC_PARAMS.bOneBeta    = bOneBeta;

    /* If calibration is done with a fixed omega, set its value */
    STATIC_PARAMS.dFixedOmega  = dFixedOmega;
    STATIC_PARAMS.bFreezeOmega = bFreezeOmega;
    STATIC_PARAMS.bOneOmega    = bOneOmega;

    /* If calibration is done with a fixed corr, set alpha and beta */
    STATIC_PARAMS.dFixedAlpha = dFixedAlpha;
    STATIC_PARAMS.dFixedGamma = dFixedGamma;
    STATIC_PARAMS.dFixedRho   = dFixedRho;

    /* For Annealing and Simplex, set the wight and target for each pdData */
    STATIC_PARAMS.pdMarketWeights = pdMarketWeights;
    STATIC_PARAMS.pdMarketTargets = pdMarketTargets;

    /* Sets the parameters ppdParamBounds */
    STATIC_PARAMS.ppdParamBounds = ppdParamBounds;

    /* Sets the model type and dimension */
    STATIC_PARAMS.eModelType = eModelType;
    STATIC_PARAMS.eModelDim  = eModelDim;

    /* Need to smooth sigma or not */
    STATIC_PARAMS.bSmoothSigma = bSmoothSigma;

    /* Stores the name of the underlying */
    strcpy(STATIC_PARAMS.szUndName, szUndName);

    /* Return a success message */
    return NULL;
}

/* ----------------------------------------------------------------------- */
#ifdef PVMPI
/* This function is needed on the slaves only */
/* Note we do not check for the bounds but if
param is in bound param+epsilon should also be,,,*/
void shift_param(double* param, int sign)
{
    *param += (*param) ? sign * PROP_SHIFT * (*param) : sign * ZERO_SHIFT;
}
#endif
/* ------------------------------------------------------------------------
This function allocate space for the Sigma and Tau pdValues stored in the
STATIC_PARAMS structure.
These pdValues will be populated by the call to  the
from_opparams_to_sigmatau function
------------------------------------------------------------------------ */

#ifdef PVMPI
Err allocate_space_for_static_sigtau_for_calib
#else
static Err allocate_space_for_static_sigtau_for_calib
#endif
    (SrtMdlDim eModelDim, long lNumSigmas, double* pdSigmaDates, long lNumTaus, double* pdTauDates)
{
    long i;
#ifdef PVMPI
    if (!PRL_CONTEXT.bLocalComp)
        SendParamForAllocation(eModelDim, lNumSigmas, pdSigmaDates, lNumTaus, pdTauDates);
#endif

    /* Sets static ts_lengths according to lNumSigmas and lNumTaus */
    STATIC_PARAMS.lNumSigmas = lNumSigmas;
    STATIC_PARAMS.lNumTaus   = lNumTaus;

    if (eModelDim == TWO_FAC)
    {
        STATIC_PARAMS.lNumSigmaCols = 6;
        STATIC_PARAMS.lNumTauCols   = 3;

        /* Allocate space for the static parameters Sigmas (6 cols) and Taus (3 cols) */
        STATIC_PARAMS.ppdSigmaValues = dmatrix(0, 5, 0, lNumSigmas - 1);
        STATIC_PARAMS.ppdTauValues   = dmatrix(0, 2, 0, lNumTaus - 1);
    }

    else if (eModelDim == ONE_FAC)
    {
        STATIC_PARAMS.lNumSigmaCols = 3;
        STATIC_PARAMS.lNumTauCols   = 2;

        /* Allocate space for the static parameters Sigmas (3 cols) and Taus (2 cols) */
        STATIC_PARAMS.ppdSigmaValues = dmatrix(0, 2, 0, lNumSigmas - 1);
        STATIC_PARAMS.ppdTauValues   = dmatrix(0, 1, 0, lNumTaus - 1);
    }

    /* Transfers the Dates  */
    for (i = 0; i < lNumSigmas; i++)
        (STATIC_PARAMS.ppdSigmaValues)[0][i] = pdSigmaDates[i];
    for (i = 0; i < lNumTaus; i++)
        (STATIC_PARAMS.ppdTauValues)[0][i] = pdTauDates[i];

    return NULL;
}

/* ----------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
This function freea the space allocated for the Sigma and Tau pdValues
stored in the STATIC_PARAMS structure.
----------------------------------------------------------------------- */

#ifdef PVMPI
Err free_space_for_static_sigtau_for_calib
#else
static Err free_space_for_static_sigtau_for_calib
#endif
    (SrtMdlDim eModelDim, long lNumSigmas, long lNumTaus)
{
#ifdef PVMPI
    if (!PRL_CONTEXT.bLocalComp)
        SendParamToDeallocate(eModelDim, lNumSigmas, lNumTaus);
#endif

    if (eModelDim == TWO_FAC)
    {
        /* Free space for the static parameters Sigmas (6 cols) and Taus (3 cols) */
        if (STATIC_PARAMS.ppdSigmaValues)
            free_dmatrix(STATIC_PARAMS.ppdSigmaValues, 0, 5, 0, lNumSigmas - 1);
        STATIC_PARAMS.ppdSigmaValues = NULL;
        if (STATIC_PARAMS.ppdTauValues)
            free_dmatrix(STATIC_PARAMS.ppdTauValues, 0, 2, 0, lNumTaus - 1);
        STATIC_PARAMS.ppdTauValues = NULL;
    }

    else if (eModelDim == ONE_FAC)
    {
        /* Free space for the static parameters Sigmas (3 cols) and Taus (2 cols) */
        if (STATIC_PARAMS.ppdSigmaValues)
            free_dmatrix(STATIC_PARAMS.ppdSigmaValues, 0, 2, 0, lNumSigmas - 1);
        STATIC_PARAMS.ppdSigmaValues = NULL;
        if (STATIC_PARAMS.ppdTauValues)
            free_dmatrix(STATIC_PARAMS.ppdTauValues, 0, 1, 0, lNumTaus - 1);
        STATIC_PARAMS.ppdTauValues = NULL;
    }

    if ((STATIC_PARAMS.bFreezeTau == SRT_YES) && (STATIC_PARAMS.bFreezeTau))
    {
        free_dvector(STATIC_PARAMS.dFixedTau, 0, lNumTaus - 1);
        STATIC_PARAMS.dFixedTau = NULL;
    }

    return NULL;
}

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
PART 4
PRICING FUNCTIONS AS REQUIRED BY THE NRC ROUTINES
(THEY ALL USE THE STATIC STRUCTURES)
----------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
Gices the correlation for the corresponding corr_index, given the
vector of optimisation parameters
(remember that pdFraMaturities starts at [0] )
This function is only called for a GLOBAL calibration:
alpha, beta and rho are constants (not time dependent)
If the correlation returned is not numerical, returns an error
------------------------------------------------------------------------ */

static Err corr_value(
    long    corr_index,
    double* pdOptParams, /* From [1] to [STATIC_PARAMS.lNumParams] */
    double* corr)
{
    long first_index;
    long second_index;

    corr_index -= 1;
    first_index  = DTOL(floor((double)corr_index / STATIC_INSTRUMENT.lNumTenors));
    second_index = corr_index - (first_index * STATIC_INSTRUMENT.lNumTenors);

    *corr = model_correl(
        STATIC_INSTRUMENT.pdFraMaturities[first_index],
        STATIC_INSTRUMENT.pdFraMaturities[second_index],
        pdOptParams[STATIC_PARAMS.lNumParams - 2], /* Alpha */
        pdOptParams[STATIC_PARAMS.lNumParams - 1], /* Beta */
        pdOptParams[STATIC_PARAMS.lNumParams]);    /* Rho */

    if (*corr > DBL_MAX)
        return serror("Calibration failed: get a corr overflow: %f", *corr);

    return NULL;
}

/* ------------------------------------------------------------------------ */
static Err model_value(
    long    data_index,
    double* pdOptParams, /* From [1] to [STATIC_PARAMS.lNumParams] */
    double* value);
/* --------------------------------------------------------------------
        Price the deals input as a calibration
        instrument from the spreadsheet, calling grfn.
        If the pdMktPrice returned is not numerical, returns an error
------------------------------------------------------------------------ */
#ifndef PVMPI
static Err aggregate_price_value
#else
Err        aggregate_price_value
#endif
    (double* pdOptParams, double* pdMktPrice)
{
    Err          err = NULL;
    SrtUndPtr    sUndPtr;
    SrtIOStruct* iolist;
    SrtMdlType   mdl_type;
    SrtMdlDim    mdl_dim;
    long         i;
    double*      pdTempPrice = NULL;
    TermStruct*  psTermStruct;

    /* Transforms the optimisation parameters in full TermStruct of sig,tau... */
    err = from_optparam_to_sigtau(
        pdOptParams, /* From [1] to [lNumParams] */
        STATIC_PARAMS.ppdSigmaValues,
        STATIC_PARAMS.lNumSigmas,
        STATIC_PARAMS.ppdTauValues,
        STATIC_PARAMS.lNumTaus,
        STATIC_PARAMS.bFreezeTau,
        STATIC_PARAMS.bOneTau,
        STATIC_PARAMS.dFixedTau,
        STATIC_PARAMS.bFreezeBeta,
        STATIC_PARAMS.bOneBeta,
        STATIC_PARAMS.dFixedBeta,
        STATIC_PARAMS.eModelDim,
        STATIC_PARAMS.bFreezeOmega,
        STATIC_PARAMS.bOneOmega,
        STATIC_PARAMS.dFixedOmega,
        STATIC_PARAMS.eCalibType,
        STATIC_PARAMS.dFixedAlpha,
        STATIC_PARAMS.dFixedGamma,
        STATIC_PARAMS.dFixedRho);
    if (err)
        return err;

    /* Initialises the TermStruct */
    err = srt_f_init_IRM_TermStruct(
        STATIC_INSTRUMENT.lClcnDate,
        STATIC_PARAMS.ppdSigmaValues,
        STATIC_PARAMS.lNumSigmaCols,
        STATIC_PARAMS.lNumSigmas,
        STATIC_PARAMS.ppdTauValues,
        STATIC_PARAMS.lNumTauCols,
        STATIC_PARAMS.lNumTaus,
        STATIC_PARAMS.eModelType,
        STATIC_PARAMS.eModelDim,
        0.0, /* Beta */
        0.0, /* Alpha */
        0.0, /* Gamma */
        0.0, /* Rho */
        0.0, /* Vovol */
        0.0, /* Eta */
        0.0, /* vasicek parms */
        0,
        0,
        NULL,
        &psTermStruct);
    if (err)
    {
        return err;
    }

    /* Attaches the termstruct to the underlying */
    sUndPtr = lookup_und(STATIC_PARAMS.szUndName);
    if (sUndPtr == NULL)
        return serror("Underlying %s not defined", STATIC_PARAMS.szUndName);
    set_irund_ts(sUndPtr, psTermStruct);

    /* Gets the Model type (LGM, Cheyette,...) */
    err = get_underlying_mdltype(sUndPtr, &mdl_type);

    /* Gets the Number of Factors in the Model */
    err = get_underlying_mdldim(sUndPtr, &mdl_dim);

    /* IF LGM or EtaBeta, there is a REAL closed form: go to it */
    if ((mdl_type == ETABETA || mdl_type == LGM))
    {
        for (i = 0; i < STATIC_INSTRUMENT.lNumInstruments; i++)
        {
            err = srt_f_closed_form(
                sUndPtr,
                mdl_type,
                mdl_dim,
                &(STATIC_INSTRUMENT.psSwapDp[i]),
                STATIC_INSTRUMENT.pdStrike[i],
                STATIC_INSTRUMENT.pdBondStrike[i],
                STATIC_INSTRUMENT.peRecPay[i],
                STATIC_INSTRUMENT.peProductType[i],
                STATIC_INSTRUMENT.pszRefRateCode[i],
                &(pdMktPrice[i]));

            if (err)
            {
                free_underlying_ts(sUndPtr);
                return err;
            }
        }
    }
    else
    {
        /* There is no closed form available: have to build a full grfn tableau */
        if (mdl_type == BDT)
        {
            free_underlying_ts(sUndPtr);
            return serror("Calibration failed: BDT model is not avaible");
        }

        err = srt_f_IOstructcreate(&iolist, "aggregate");

        /* Call Grfn with this tableau */
        if (!err)
        {
            err = srt_f_grfn(
                sUndPtr,
                STATIC_INSTRUMENT.psGrfnParams,
                STATIC_SPRDSHT.lNumEventDates,
                &STATIC_SPRDSHT.pszEventDates,
                &STATIC_SPRDSHT.lNumRows,
                &STATIC_SPRDSHT.lNumCols,
                &STATIC_SPRDSHT.pshtSHEET,
                0,
                0,
                0,
                0,
                0,
                iolist,
                0,
                0);
        }
        /* Computes the value (pdMktPrice or corr) for the given set of parameters */
        pdTempPrice = dvector(0, STATIC_INSTRUMENT.lNumInstruments - 1);

        err = srt_f_IOstructgetcolpvs((*iolist), &pdTempPrice, &STATIC_SPRDSHT.lNumCols);
        if (err)
        {
            err = srt_f_IOstructfree(&iolist);
            free_underlying_ts(sUndPtr);
            return err;
        }

        for (i = 0; i < STATIC_INSTRUMENT.lNumInstruments; i++)
            pdMktPrice[i] = pdTempPrice[i];

        if (iolist)
            err = srt_f_IOstructfree(&iolist);

        pdTempPrice = NULL;
        free_underlying_ts(sUndPtr);
        if (err)
            return err;
    }

    for (i = STATIC_INSTRUMENT.lNumInstruments; i < STATIC_PARAMS.lNumData; i++)
    {
        err = model_value(i + 1, pdOptParams, &pdMktPrice[i]);
        if (err)
            return err;
    }

    for (i = 0; i < STATIC_INSTRUMENT.lNumInstruments; i++)
    {
        if (pdMktPrice[i] > DBL_MAX)
            return serror("Calibration failed: get a pdMktPrice overflow: %f", *pdMktPrice);
    }

    return err;
}

/* ------------------------------------------------------------------------
Price the deal number (price_index - 1) input as a calibration
instrument from the spreadsheet, calling grfn.
If the pdMktPrice returned is not numerical, returns an error
------------------------------------------------------------------------ */

static Err price_value(long price_index, double* pdMktPrice)
{
    Err       err = NULL;
    SrtUndPtr sUndPtr;

    /* Gets the underlying */
    sUndPtr = lookup_und(STATIC_PARAMS.szUndName);

    err = srt_f_grfn_clsdfrm(
        sUndPtr,
        STATIC_INSTRUMENT.psGrfnParams,
        &(STATIC_INSTRUMENT.psSwapDp[price_index - 1]),
        STATIC_INSTRUMENT.pdStrike[price_index - 1],
        STATIC_INSTRUMENT.pdBondStrike[price_index - 1],
        STATIC_INSTRUMENT.peRecPay[price_index - 1],
        STATIC_INSTRUMENT.peProductType[price_index - 1],
        STATIC_INSTRUMENT.pszRefRateCode[price_index - 1],
        pdMktPrice);
    if (err)
        return err;

    if (*pdMktPrice > DBL_MAX)
        return serror("Calibration failed: get a pdMktPrice overflow: %f", *pdMktPrice);

    return err;
}

/* ----------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
Computes stdev of sigma function
------------------------------------------------------------------------ */

static double std_sig_value(void)
{
    long   i;
    double exp_sigma = 0.0, exp_sigma2 = 0.0;

    for (i = 1; i <= STATIC_PARAMS.lNumSigmas; i++)
    {
        exp_sigma += STATIC_PARAMS.ppdSigmaValues[1][i - 1];
        exp_sigma2 +=
            STATIC_PARAMS.ppdSigmaValues[1][i - 1] * STATIC_PARAMS.ppdSigmaValues[1][i - 1];
    }

    exp_sigma /= STATIC_PARAMS.lNumSigmas;
    exp_sigma2 /= STATIC_PARAMS.lNumSigmas;

    return sqrt(DMAX(exp_sigma2 - exp_sigma * exp_sigma, 0.0));
}

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
According to the data_index, gives the model value of:
if 1<=data_index <= num_instr: the instrument pdMktPrice
if data_index > num_instr: the correlation value
Parameters must have been checked before (within ppdParamBounds)
----------------------------------------------------------------------- */
static Err model_value(
    long    data_index,
    double* pdOptParams, /* From [1] to [STATIC_PARAMS.lNumParams] */
    double* value)
{
    Err    err = NULL;
    long   price_index;
    long   corr_index;
    double pdMktPrice;
    double corr;

    if (data_index <= STATIC_INSTRUMENT.lNumInstruments)
    {
        price_index = data_index;
        err         = price_value(price_index, &pdMktPrice);
        *value      = pdMktPrice;
    }
    else if (
        data_index <= STATIC_INSTRUMENT.lNumInstruments +
                          STATIC_INSTRUMENT.lNumTenors * STATIC_INSTRUMENT.lNumTenors)
    {
        corr_index = data_index - STATIC_INSTRUMENT.lNumInstruments;
        err        = corr_value(corr_index, pdOptParams, &corr);
        *value     = corr;
    }
    else
    {
        *value = std_sig_value();
    }

    return err;
}

/* ------------------------------------------------------------------------ */

/* -----------------------------------------------------------------------
For a given set of parameters pdOptParams[1..lNumParams], this function
computes the value of the global criteria for calibration, incorporating
both pricing and correlation (if required) as :
Sum(i) [(Price_i - mktprice_i)/pdMarketWeights_i]^2
+ Sum(j) [(Corr_i - mktcorr_i)/pdMarketWeights_i]^2
(i.e.dOptionsWeight * varerr to frd + (1-dOptionsWeight) * varerr to correl )

  Also, if the parameters go outside the ppdParamBounds, the function returns
  a SRT_BIG value that prevents the algorithm from using that point.

        THIS FUNCTION IS DECLARED IN THE NUMERICAL RECIPES IN C FORMAT FOR:
        SIMULATED ANNEALING
        SIMPLEX
        OUR SOBOL...
----------------------------------------------------------------------- */

static double global_criteria_funcs(
    double* pdOptParams) /* From [1] to [STATIC_PARAMS.lNumParams] */
{
    Err     err;
    long    i;
    double  diff;
    double  criteria;
    double* pdMktPrice;

    criteria = 0.0;

    /* Checks the current value of the parameters is alright; if not: SRT_BIG */
    if (are_calib_parameters_within_band(
            pdOptParams, STATIC_PARAMS.lNumParams, STATIC_PARAMS.ppdParamBounds) == SRT_NO)
    {
        criteria = SRT_BIG;
        return criteria;
    }

    pdMktPrice = dvector(0, STATIC_PARAMS.lNumData - 1);
    err        = aggregate_price_value(pdOptParams, pdMktPrice);
    if (err)
    {
        free_dvector(pdMktPrice, 0, STATIC_PARAMS.lNumData - 1);
        criteria = SRT_BIG;
        return criteria;
    }

    for (i = 1; i < STATIC_PARAMS.lNumData; i++)
    {
        diff = pdMktPrice[i - 1] - STATIC_PARAMS.pdMarketTargets[i];
        criteria +=
            (diff / STATIC_PARAMS.pdMarketWeights[i]) * (diff / STATIC_PARAMS.pdMarketWeights[i]);
    }

    free_dvector(pdMktPrice, 0, STATIC_PARAMS.lNumData - 1);
    pdMktPrice = NULL;

    return criteria;
}

/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
    This function computes one model value  for a given set of
    optimisation parameters.

      The TermStruct initialisation and freeing is done inside
      Also, if the parameters go outside the ppdParamBounds, the function returns
      a SRT_BIG value that prevents the algorithm from using that point.

            This function is called by levenberg_calib_funcs
    ---------------------------------------------------------------------- */

static Err levenberg_price_funcs(
    double  data_index,
    double  pdOptParams[], /* From [1] to [lNumParams] */
    double* value,
    int     lNumParams)
{
    Err         err;
    TermStruct* psTermStruct;
    SrtUndPtr   sUndPtr;

    /* Checks the current value of the parameters is alright; if not: SRT_BIG */
    if (are_calib_parameters_within_band(pdOptParams, lNumParams, STATIC_PARAMS.ppdParamBounds) ==
        SRT_NO)
    {
        *value = SRT_BIG;
        return NULL;
    }

    /* Transforms the optimisation parameters in full TermStruct of sig,tau... */
    err = from_optparam_to_sigtau(
        pdOptParams, /* From [1] to [lNumParams] */
        STATIC_PARAMS.ppdSigmaValues,
        STATIC_PARAMS.lNumSigmas,
        STATIC_PARAMS.ppdTauValues,
        STATIC_PARAMS.lNumTaus,
        STATIC_PARAMS.bFreezeTau,
        STATIC_PARAMS.bOneTau,
        STATIC_PARAMS.dFixedTau,
        STATIC_PARAMS.bFreezeBeta,
        STATIC_PARAMS.bOneBeta,
        STATIC_PARAMS.dFixedBeta,
        STATIC_PARAMS.eModelDim,
        STATIC_PARAMS.bFreezeOmega,
        STATIC_PARAMS.bOneOmega,
        STATIC_PARAMS.dFixedOmega,
        STATIC_PARAMS.eCalibType,
        STATIC_PARAMS.dFixedAlpha,
        STATIC_PARAMS.dFixedGamma,
        STATIC_PARAMS.dFixedRho);
    if (err)
        return err;

    /* Initialises the TermStruct */
    err = srt_f_init_IRM_TermStruct(
        STATIC_INSTRUMENT.lClcnDate,
        STATIC_PARAMS.ppdSigmaValues,
        STATIC_PARAMS.lNumSigmaCols,
        STATIC_PARAMS.lNumSigmas,
        STATIC_PARAMS.ppdTauValues,
        STATIC_PARAMS.lNumTauCols,
        STATIC_PARAMS.lNumTaus,
        STATIC_PARAMS.eModelType,
        STATIC_PARAMS.eModelDim,
        0.0, /* Beta */
        0.0, /* Alpha */
        0.0, /* Gamma */
        0.0, /* Rho */
        0.0, /* Vovol */
        0.0, /* Eta */
        0.0, /* vasicek parms */
        0,
        0,
        NULL,
        &psTermStruct);
    if (err)
    {
        return err;
    }

    /* Attaches the termstruct to the underlying */
    sUndPtr = lookup_und(STATIC_PARAMS.szUndName);
    if (sUndPtr == NULL)
        return serror("Underlying %s not defined", STATIC_PARAMS.szUndName);
    set_irund_ts(sUndPtr, psTermStruct);

    /* Computes the value (pdMktPrice or corr) for the given set of parameters */
    err = model_value((int)data_index, pdOptParams, value);
    if (err)
    {
        free_underlying_ts(sUndPtr);
        return err;
    }

    /* Free the TermStruct from the underlying */
    err = free_underlying_ts(sUndPtr);
    if (err)
        return err;

    /* Return a success string */
    return NULL;
}

/* ------------------------------------------------------------------------ */

/* -----------------------------------------------------------------------
    For a given set of parameters pdOptParams[1..lNumParams], this function
    evaluates in *value the mkt_value required (correlation or pdMktPrice) for
    the data_index input (i.e pdMktPrice or correlation number).
    The derivatives are also returned in deriv[1..lNumParams] for the
    same data_index.

      This function calls in a loop levenberg_price_funcs, that evaluates
      only the *value for a given set of parameters.

    THIS FUNCTION IS DECLARED IN THE NUMERICAL RECIPES IN C FORMAT FOR:
            LEVENBERG_MARQUARDT
            ----------------------------------------------------------------------- */

#ifdef PVMPI
Err levenberg_calib_funcs(
#else
static Err levenberg_calib_funcs(
#endif
    double  instr_index,   /* Equivalent of pdData[] */
    double  pdOptParams[], /* From [1] to [lNumParams] */
    double* value,
    double  deriv[],
    int     lNumParams)
{
    long   i;
    double shift;
    Err    err;

    /* Computes the value (pdMktPrice or corr) for the given set of parameters */

    err = levenberg_price_funcs(
        instr_index,
        pdOptParams, /* From [1] to [lNumParams] */
        value,
        lNumParams);
    if (err)
        return err;

    GlobalTheoPrices[(int)instr_index] = *value;

    /*	if (*value == SRT_BIG)
    {
    smessage("No valid pdMktPrice or correlation - check tau,alpha,gamma");
    }
*/

    /* Computes the derivatives of this value with respect to each parameter */
    for (i = 1; i <= lNumParams; i++)
    {
        if (pdOptParams[i] == 0.0)
            shift = ZERO_SHIFT;
        else
            shift = PROP_SHIFT * pdOptParams[i];
        pdOptParams[i] += shift;

        /* Computes the shifted pdMktPrice: stores it in deriv */
        err = levenberg_price_funcs(
            instr_index,
            pdOptParams, /* From [1] to [lNumParams] */
            &(deriv[i]),
            lNumParams);
        if (err)
            return err;
        /*
        if (*value == SRT_BIG)
        {
        smessage("No valid pdMktPrice or correlation - check tau,alpha,gamma");
        }

*/

        /* Derivative = (shift_price - pdMktPrice) / shift */
        deriv[i] -= *value;
        deriv[i] /= shift;

        /* resets params to its initial value */
        pdOptParams[i] -= shift;
    }

    /* Return a success string */
    return NULL;
}

/* -----------------------------------------------------------------------
For a given set of parameters pdOptParams[1..lNumParams], this function
evaluates in *value the mkt_value required (correlation or pdMktPrice) for
all data_index input (i.e pdMktPrice or correlation number).
The derivatives are also returned in deriv[1..lNumParams] for all data_index.

  THIS FUNCTION IS DECLARED IN THE NUMERICAL RECIPES IN C FORMAT FOR:
  LEVENBERG_MARQUARDT
----------------------------------------------------------------------- */

#ifdef PVMPI
Err levenberg_aggregate_calib_funcs(
#else
static Err levenberg_aggregate_calib_funcs(
#endif
    double  instr_index,   /* Equivalent of pdData[] */
    double  pdOptParams[], /* From [1] to [lNumParams] */
    double* value,
    double  deriv[],
    int     lNumParams)
{
    long               i, j;
    double             shift;
    Err                err;
    static SRT_Boolean bFirstTime = SRT_TRUE;
#ifndef PVMPI
    double* pdTempPrice;
#endif
    int nb_retries = 0;

    bFirstTime = (instr_index == 1) ? SRT_TRUE : SRT_FALSE;

    if (bFirstTime)
    {
        /* Checks the current value of the parameters is alright; if not: SRT_BIG */
        if (are_calib_parameters_within_band(
                pdOptParams, lNumParams, STATIC_PARAMS.ppdParamBounds) == SRT_NO)
        {
            for (i = 0; i < STATIC_PARAMS.lNumData; i++)
                STATIC_MODELVALUES.pdModelValue[i] = SRT_BIG;
            return NULL;
        }
#ifndef PVMPI
        pdTempPrice = dvector(0, STATIC_PARAMS.lNumData - 1);
        err         = aggregate_price_value(pdOptParams, pdTempPrice);
        if (err)
        {
            free_dvector(pdTempPrice, 0, STATIC_PARAMS.lNumData - 1);
            pdTempPrice = NULL;
            return err;
        }
        for (i = 0; i < STATIC_PARAMS.lNumData; i++)
        {
            STATIC_MODELVALUES.pdModelValue[i] = pdTempPrice[i];
            GlobalTheoPrices[i + 1]            = pdTempPrice[i];
        }
        for (i = 1; i <= STATIC_PARAMS.lNumParams; i++)
        {
            if (pdOptParams[i] == 0.0)
                shift = ZERO_SHIFT;
            else
                shift = PROP_SHIFT * pdOptParams[i];
            pdOptParams[i] += shift;
            /* Checks the current value of the parameters is alright; if not: SRT_BIG */
            if (are_calib_parameters_within_band(
                    pdOptParams, lNumParams, STATIC_PARAMS.ppdParamBounds) == SRT_NO)
            {
                free_dvector(pdTempPrice, 0, STATIC_PARAMS.lNumData - 1);
                pdTempPrice = NULL;
                return NULL;
            }
            err = aggregate_price_value(pdOptParams, pdTempPrice);
            if (err)
            {
                free_dvector(pdTempPrice, 0, STATIC_PARAMS.lNumData - 1);
                pdTempPrice = NULL;
                return err;
            }

            /* Calculate the Hessian once. */
            for (j = 0; j < STATIC_PARAMS.lNumData; j++)
            {
                STATIC_MODELVALUES.pdModelHessian[j][i - 1] =
                    (pdTempPrice[j] - STATIC_MODELVALUES.pdModelValue[j]) / shift;
            }

            err = aggregate_price_value(pdOptParams, pdTempPrice);
            if (err)
            {
                free_dvector(pdTempPrice, 0, STATIC_PARAMS.lNumData - 1);
                pdTempPrice = NULL;
                return err;
            }

            /* resets params to its initial value */
            pdOptParams[i] -= shift;
        }
        free_dvector(pdTempPrice, 0, STATIC_PARAMS.lNumData - 1);
        pdTempPrice = NULL;
#else
        /* Send the set of parameters and get the results back */
        nb_retries = 0;
        do
        {
            err =
                ParallelComputation(STATIC_PARAMS.lNumData, STATIC_PARAMS.lNumParams, pdOptParams);
            if (!strcmp(err, "ok"))
            {
                for (i = 0; i < STATIC_PARAMS.lNumParams; i++)
                {
                    if (pdOptParams[i + 1] == 0.0)
                        shift = ZERO_SHIFT;
                    else
                        shift = PROP_SHIFT * pdOptParams[i + 1];
                    for (j = 0; j < STATIC_PARAMS.lNumData; j++)
                    {
                        STATIC_MODELVALUES.pdModelHessian[j][i] -=
                            STATIC_MODELVALUES.pdModelValue[j];
                        STATIC_MODELVALUES.pdModelHessian[j][i] /= shift;
                    }
                }
                for (j = 0; j < STATIC_PARAMS.lNumData; j++)
                    GlobalTheoPrices[j + 1] = STATIC_MODELVALUES.pdModelValue[j];
                break;
            }
            else if (!strcmp(err, "Aborting !"))
            {
                smessage("All slaves failing. Aborting.");
                return err;
            }
            else if (!strcmp(err, "Same player shoot again !"))
            {
                smessage("Some errors. Trying again....");
                err = RemoteErrorRecovery();
                /* if the same data fails 3 times drop the whole thing and abort */
                if (++nb_retries > 3)
                    return "Too many successive iterations failed";
                if (err)
                    return err;
            }
        } while (1);
#endif
#ifdef PVMPI
#ifdef _DEBUG
        MyTRACE("");
        MyTRACE("Parameters: ");
        for (i = 0; i < STATIC_PARAMS.lNumParams; i++)
            MyTRACE1(" ", pdOptParams[i + 1]);
        MyTRACE("");
        MyTRACE("Prices: ");
        for (i = 0; i < STATIC_PARAMS.lNumData; i++)
            MyTRACE1(" ", STATIC_MODELVALUES.pdModelValue[i]);
        MyTRACE("");
        MyTRACE("Hessian: ");
        for (j = 0; j < STATIC_PARAMS.lNumParams; j++)
        {
            MyTRACE("");
            for (i = 0; i < STATIC_PARAMS.lNumData; i++)
                MyTRACE1(" ", STATIC_MODELVALUES.pdModelHessian[i][j]);
        }
        MyTRACE("");
#endif
#endif
    }
    /* End of the bFirstTime */

    *value = STATIC_MODELVALUES.pdModelValue[(long)instr_index - 1];
    for (i = 1; i <= STATIC_PARAMS.lNumParams; i++)
        deriv[i] = STATIC_MODELVALUES.pdModelHessian[(long)instr_index - 1][i - 1];
    /* Return a success string */
    return NULL;
}

/* ---------------------------------------------------------------------- */

/* ======================================================================== */

/* -----------------------------------------------------------------------
    PART 5
    THE CORE FUNCTION
    ------------------------------------------------------------------------ */

#define free_calib_core_memory                                                                     \
    {                                                                                              \
        if (pdMarketWeights)                                                                       \
            free_dvector(pdMarketWeights, 1, lNumData);                                            \
        if (pdMarketTargets)                                                                       \
            free_dvector(pdMarketTargets, 1, lNumData);                                            \
        if (ppdParamBounds)                                                                        \
            free_dmatrix(ppdParamBounds, 1, 2, 1, lNumParams);                                     \
        if (pdOptParams)                                                                           \
            free_dvector(pdOptParams, 1, lNumParams);                                              \
        free_space_for_static_sigtau_for_calib(eModelDim, lNumSigmas, lNumTaus);                   \
        free_space_for_static_sprdsht_for_calib(STATIC_SPRDSHT.lNumRows, STATIC_SPRDSHT.lNumCols); \
    }

Err srt_f_calib_core(
    SrtGrfnParam*    psGrfnParams,
    SrtModelType     eModelType,
    SrtMdlDim        eModelDim,
    SwapDP*          psSwapDp,
    double*          pdStrike,
    double*          pdBondStrike,
    StructType*      peProductType,
    SrtReceiverType* peRecPay,
    double*          pdMktPrice,
    double*          pdMktVega,
    String*          pszRefRateCode,
    int              lNumInstruments,
    double*          pdFraMaturities,
    double**         ppdCorrelationMatrix,
    long             lNumTenors,
    SrtUndPtr        sUndPtr,
    SrtCalibParam*   psCalibParams,
    double**         ppdSigmaValues,
    long             lNumSigmas,
    long             lNumSigmaCols,
    double**         ppdTauValues,
    long             lNumTaus,
    long             lNumTauCols,
    double*          pdChiSquare)
{
    long        i, j;
    long        lClcnDate;
    long        lNumParams;
    long        lNumData;               /* Number of pdData to be used in the criteria */
    double*     pdMarketWeights = NULL; /* For the minimisation criteria in LM */
    double*     pdMarketTargets = NULL; /* For the minimisation criteria in LM */
    double*     pdOptParams     = NULL;
    double*     pdData          = NULL;
    double**    ppdParamBounds  = NULL;
    double**    dFixedTau;
    double      dFixedBeta;
    double      dFixedOmega;
    double      dFixedAlpha;
    double      dFixedGamma;
    double      dFixedRho;
    double**    simplex_p = NULL;
    double*     pdValues  = NULL;
    double      best;
    long        n_func_eval;
    TermStruct* psTermStruct;
    Err         err;
    long        lNumUsedTaus;
    long        lNumUsedBetas;
    long        lNumUsedOmegas;
    String      szUndName;
    long        lStartTime;
    long        lEndTime;

    /* ------------------ CALIBRATION CHARACTERISATION ---------------------- */

    /* There is at least the Sigma Term Sturct to calibrate to the market prices of options */
    lNumParams   = lNumSigmas;
    lNumData     = lNumInstruments;
    dFixedTau    = (double**)malloc(sizeof(double*));
    dFixedTau[0] = NULL;

    /* Set the number of used taus and the value of fixed tau (if relevant) */
    if (psCalibParams->bFreezeTau == SRT_YES)
    {
        /* Calibration is done with a frozen Tau Struct */
        dFixedTau[0] = dvector(0, lNumTaus - 1);
        for (i = 0; i < lNumTaus; i++)
            dFixedTau[0][i] = ppdTauValues[1][i];
        /*dFixedTau = ppdTauValues[1][0];*/
        lNumUsedTaus = 0;
    }
    else
    {
        if (psCalibParams->bOneTau == SRT_YES)
            lNumUsedTaus = 1;
        else
            lNumUsedTaus = lNumTaus;
    }
    lNumParams += lNumUsedTaus;

    /* IF LGM or CHEYETTE (without Beta) : FreezeBeta makes sense */
    if ((eModelType == LGM) || (eModelType == CHEY))
    {
        psCalibParams->bFreezeBeta = SRT_YES;
    }

    /* Set the number of used betas and the value of fixed beta (if relevant) */
    if (psCalibParams->bFreezeBeta == SRT_YES)
    {
        /* Calibration is done with a frozen Beta */
        dFixedBeta    = ppdSigmaValues[2][0];
        lNumUsedBetas = 0;
    }
    else
    {
        if (psCalibParams->bOneBeta == SRT_YES)
            lNumUsedBetas = 1;
        else
            lNumUsedBetas = lNumSigmas;
    }
    lNumParams += lNumUsedBetas;

    /* For Two Factor models, some extra parameters might be required */
    if (eModelDim == TWO_FAC)
    {
        /* If not MIXED_BETA, FreezeOmega makes sense */
        if (eModelType != CHEY_BETA)
        {
            psCalibParams->bFreezeOmega = SRT_YES;
        }

        /* Set the number of used omegas and the value of fixed omega (if relevant) */
        if (psCalibParams->bFreezeOmega == SRT_YES)
        {
            /* Calibration is done with a frozen Omega ( = Beta 2 - Beta 1) */
            dFixedOmega    = ppdSigmaValues[4][0] - ppdSigmaValues[2][0];
            lNumUsedOmegas = 0;
        }
        else
        {
            if (psCalibParams->bOneOmega == SRT_YES)
                lNumUsedOmegas = 1;
            else
                lNumUsedOmegas = lNumSigmas;
        }

        lNumParams += lNumUsedOmegas;

        /* Sets the correlation parameters */
        if (psCalibParams->eCalibType == GLOBAL_CALIB)
        {
            lNumParams += 3;
            lNumData += lNumTenors * lNumTenors;
        }
        else if (psCalibParams->eCalibType == FIXED_CALIB)
        {
            dFixedAlpha = ppdSigmaValues[3][0] / ppdSigmaValues[1][0];
            dFixedGamma = 1.0 / ppdTauValues[2][0] - 1.0 / ppdTauValues[1][0];
            dFixedRho   = ppdSigmaValues[5][0];
        }

    } /* END if (eModelDim == TWO_FAC) */

    /* If the sigma smoothing is also a criteria, add one more pdData to "fit" (the sigam variance)
     */
    if (psCalibParams->bSmoothSigma == SRT_YES)
    {
        lNumData++;
    }

    /* ------------------ VARIABLES INITIALISATION ---------------------- */

    /* Sets weigths for minimisation criteria using vegas, options weights and sizes */
    pdMarketWeights = dvector(1, lNumData);
    err             = define_criteria_weights(
        lNumInstruments,
        lNumTenors,
        psCalibParams->bSmoothSigma,
        psCalibParams->dOptionsWeight,
        1.0 - psCalibParams->dOptionsWeight,
        psCalibParams->dSmoothSigmaWeight,
        pdMktVega,
        pdMarketWeights);
    if (err)
    {
        free_calib_core_memory;
        return err;
    }

    /* Sets a vector of market target for minimisation criteria  */
    pdMarketTargets = dvector(1, lNumData);
    err             = define_market_targets(
        lNumInstruments,
        psCalibParams->bSmoothSigma,
        pdMktPrice,
        ppdCorrelationMatrix,
        lNumTenors,
        lNumData,
        pdMarketTargets);
    if (err)
    {
        free_calib_core_memory;
        return err;
    }

    /* Sets an array of ppdParamBounds for minimisation criteria  */
    ppdParamBounds = dmatrix(1, 2, 1, lNumParams);
    err            = set_calib_parameters_bounds(
        psCalibParams,
        ppdParamBounds,
        lNumSigmas,
        lNumUsedTaus,
        lNumUsedBetas,
        lNumUsedOmegas,
        eModelDim);
    if (err)
    {
        free_calib_core_memory;
        return err;
    }

    /* Transform sig/tau in a vector of parameters as required by NR in C*/
    pdOptParams = dvector(1, lNumParams);
    err         = from_sigtau_to_optparam(
        ppdSigmaValues,
        lNumSigmas,
        ppdTauValues,
        lNumUsedTaus,
        lNumUsedBetas,
        lNumUsedOmegas,
        psCalibParams->eCalibType,
        eModelDim,
        pdOptParams);
    if (err)
    {
        free_calib_core_memory;
        return err;
    }

    /* If Sobol is required, gets the best initial value of the param vector */
    if ((psCalibParams->lNumSobolPaths > 0) && (psCalibParams->eAlgoType == LEVENBERG_MARQUARDT))
    {
        /* Gets the starting point in a vector of parameters */
        err = sobol_startingpoint(
            pdOptParams,
            ppdParamBounds,
            lNumParams,
            psCalibParams->lNumSobolPaths,
            global_criteria_funcs,
            pdChiSquare);
        if (err)
        {
            free_calib_core_memory;
            return err;
        }
    } /* END if (psCalibParams->lNumSobolPaths > 0 ) */

    /* ------------------ STATIC VARIABLES INITIALISATION ---------------------- */

    /* Sets all useful statics parameters (for calibration routines from Numerical Recipes) */
    szUndName = get_underlying_name(sUndPtr);
    err       = set_static_params_for_calib(
        psCalibParams->eCalibType,
        lNumData,
        pdMarketWeights,
        pdMarketTargets,
        lNumParams,
        lNumSigmas,
        lNumTaus,
        ppdParamBounds,
        psCalibParams->bFreezeTau,
        psCalibParams->bOneTau,
        dFixedTau,
        psCalibParams->bFreezeBeta,
        psCalibParams->bOneBeta,
        dFixedBeta,
        psCalibParams->bFreezeOmega,
        psCalibParams->bOneOmega,
        dFixedOmega,
        dFixedAlpha,
        dFixedGamma,
        dFixedRho,
        eModelType,
        eModelDim,
        psCalibParams->bSmoothSigma,
        szUndName);
    if (err)
    {
        free_calib_core_memory;
        return err;
    }

    /* Stores all the instruments used for calibration in a static structure */
    lClcnDate = get_today_from_underlying(sUndPtr);
    err       = set_static_instr_for_calib(
        lClcnDate,
        lNumInstruments,
        psSwapDp,
        pdStrike,
        pdBondStrike,
        peRecPay,
        peProductType,
        pszRefRateCode,
        lNumTenors,
        pdFraMaturities,
        psGrfnParams);
    if (err)
    {
        free_calib_core_memory;
        return err;
    }

    /* If Model is not LGM (no closed form avaible) then build the static SpreadSheet */
    if ((eModelType == LGM) || (psCalibParams->bAggregate == SRT_NO))
    {
        STATIC_SPRDSHT.pshtSHEET          = NULL;
        STATIC_SPRDSHT.pszEventDates      = NULL;
        STATIC_MODELVALUES.pdModelValue   = NULL;
        STATIC_MODELVALUES.pdModelHessian = NULL;
    }
    else
    {
        err = set_static_sprdsht_for_calib(
            get_underlying_name(sUndPtr), get_today_from_underlying(sUndPtr));
        if (err)
        {
            free_calib_core_memory;
            return err;
        }
    }
    /* Allocate space for the STATIC_PARAMS Sigmas and Taus pdValues */
    err = allocate_space_for_static_sigtau_for_calib(
        eModelDim, lNumSigmas, ppdSigmaValues[0], lNumTaus, ppdTauValues[0]);
    if (err)
    {
        free_calib_core_memory;
        return err;
    }

    /* ------------------------ CALIBRATION LAUNCH -------------------------------- */

    /* Safety measure: free the TermStruct attached to the underlying if any */
    err = free_underlying_ts(sUndPtr);
    if (err)
    {
        free_calib_core_memory;
        return err;
    }

    /* Gets the Starting Time */
    lStartTime = time(NULL);

    if (psCalibParams->lNumIter > 0)
    {
        /* Call the relevant algorithm for optimisation */

        /* --------------------- LEVENBERG - MARQUARDT -------------------------  */
        if (psCalibParams->eAlgoType == LEVENBERG_MARQUARDT)
        {
            smessage("Calibration using Levenberg-Marquardt algorithm");
            smessage("");

            /* Sets a vector of x pdData: here, it is 1,2,3..: the market pdData index */
            pdData = dvector(1, lNumData);
            for (i = 1; i <= lNumData; i++)
                pdData[i] = (double)(i);

            /* Call the Levenberg-Marquardt main routine with the right function */
            if ((eModelType == LGM) || (psCalibParams->bAggregate == SRT_NO))
            {
                err = levenberg_marquardt(
                    pdData,
                    pdMarketTargets,
                    pdMarketWeights,
                    lNumData,
                    pdOptParams,
                    lNumParams,
                    psCalibParams->lNumIter,
                    levenberg_calib_funcs,
                    pdChiSquare);
            }
            else
            {
                err = levenberg_marquardt(
                    pdData,
                    pdMarketTargets,
                    pdMarketWeights,
                    lNumData,
                    pdOptParams,
                    lNumParams,
                    psCalibParams->lNumIter,
                    levenberg_aggregate_calib_funcs,
                    pdChiSquare);
            }
            /* Frees whatever has to be freed */
            if (pdData)
                free_dvector(pdData, 1, lNumData);

        } /* END if (psCalibParams->algorithm == LEVENBERG_MARQUARDT) */

        /* --------------------- SIMULATED_ANNEALING --------------------------  */
        else if (psCalibParams->eAlgoType == SIMULATED_ANNEALING)
        {
            smessage("Calibration using Simulated-Annealing algorithm");
            smessage("");
            smessage(" initialising...");
            smessage(" please wait...");
            smessage("");

            /* Sets simplex starting points */
            simplex_p = dmatrix(1, lNumParams + 1, 1, lNumParams);
            pdValues  = dvector(1, lNumParams + 1);
            /* First point is starting point */
            for (i = 1; i <= lNumParams; i++)
                simplex_p[1][i] = pdOptParams[i];
            /* All the other points are scaling for each parameter */
            for (i = 2; i <= lNumParams + 1; i++)
            {
                for (j = 1; j <= lNumParams; j++)
                    simplex_p[i][j] = pdOptParams[j];
                simplex_p[i][i - 1] *= psCalibParams->dSimplexInitScale;
            }

            /* Sets pdValues of global criteria at starting points */
            for (i = 1; i <= lNumParams + 1; i++)
            {
                pdValues[i] = global_criteria_funcs(simplex_p[i]);
            }

            /* Call the Simulated Annealing  main routine with the right function */
            err = simulated_annealing(
                simplex_p,
                pdValues,
                lNumParams,
                psCalibParams->dSimplexTol,
                psCalibParams->lNumIter,
                global_criteria_funcs,
                psCalibParams->dMaxTemp,
                psCalibParams->dMinTemp,
                psCalibParams->dDecreaseFactor);

            /* Extract the optimised parameters from the simplex (in first vertex) */
            *pdChiSquare = pdValues[1];
            for (i = 1; i <= lNumParams; i++)
                pdOptParams[i] = simplex_p[1][i];

            /* Free whatever has to be freed */
            if (simplex_p)
                free_dmatrix(simplex_p, 1, lNumParams + 1, 1, lNumParams);
            if (pdValues)
                free_dvector(pdValues, 1, lNumParams + 1);

        } /* END if ( psCalibParams->algorithm == SIMULATED_ANNEALING ) */

        /* --------------------- SOBENBERG -------------------------  */
        else if (psCalibParams->eAlgoType == SOBENBERG)
        {
            smessage("Calibration using Sobenberg algorithm");
            smessage("");

            /* Sets a vector of x pdData: here, it is 1,2,3..: the market pdData index */
            pdData = dvector(1, lNumData);
            for (i = 1; i <= lNumData; i++)
                pdData[i] = (double)(i);

            /* Call the Levenberg-Marquardt main routine with the right function */
            err = sobenberg(
                pdData,
                pdMarketTargets,
                pdMarketWeights,
                lNumData,
                pdOptParams,
                ppdParamBounds[1],
                ppdParamBounds[2],
                lNumParams,
                psCalibParams->lNumSobenPoints,
                psCalibParams->lNumBestSobenPoints,
                psCalibParams->lNumClusters,
                psCalibParams->szSobenMethod,
                psCalibParams->lNumIter,
                levenberg_price_funcs,
                levenberg_calib_funcs,
                pdChiSquare);

            /* Frees whatever has to be freed */
            if (pdData)
                free_dvector(pdData, 1, lNumData);

        } /* END if (psCalibParams->algorithm == LEVENBERG_MARQUARDT) */

        /* --------------------------- SIMPLEX --------------------------------  */
        else if (psCalibParams->eAlgoType == SIMPLEX)
        {
            smessage("Calibration using Simplex algorithm");
            smessage("Please wait...");
            smessage("");

            /* Sets simplex starting points: a scaling from start for each parameter */
            simplex_p = dmatrix(1, lNumParams + 1, 1, lNumParams);
            pdValues  = dvector(1, lNumParams + 1);
            /* First point is starting point */
            for (i = 1; i <= lNumParams; i++)
                simplex_p[1][i] = pdOptParams[i];
            /* All the other points are scaling for each parameter */
            for (i = 2; i <= lNumParams + 1; i++)
            {
                for (j = 1; j <= lNumParams; j++)
                    simplex_p[i][j] = pdOptParams[j];
                simplex_p[i][i - 1] *= psCalibParams->dSimplexInitScale;
            }
            /* Sets pdValues of global criteria at starting points */
            for (i = 1; i <= lNumParams + 1; i++)
            {
                pdValues[i] = global_criteria_funcs(simplex_p[i]);
            }

            /* Call the Simplex main routine with the right function to minimise */
            err = amoeba(
                simplex_p,
                pdValues,
                (int)lNumParams,
                psCalibParams->dSimplexTol,
                psCalibParams->lNumIter,
                global_criteria_funcs,
                (int*)&n_func_eval);

            /* Extracts best optimised parameters from the vertexes in the simplex */
            best = pdValues[1];
            j    = 1;
            for (i = 2; i <= lNumParams; i++)
            {
                if (pdValues[i] <= best)
                    j = i;
            }
            *pdChiSquare = pdValues[j];
            for (i = 1; i <= lNumParams; i++)
                pdOptParams[i] = simplex_p[j][i];

            /* Free whatever has to be freed */
            if (simplex_p)
                free_dmatrix(simplex_p, 1, lNumParams + 1, 1, lNumParams);
            if (pdValues)
                free_dvector(pdValues, 1, lNumParams + 1);

        } /* END if ( psCalibParams->algorithm == SIMPLEX ) */

        else
        {
            return serror("Unknown calibration algorithm ...");
        }
    } /* END 	if (psCalibParams->lNumIter > 0) */
    if (err)
    {
        free_calib_core_memory;
        return err;
    }
    /* Frees the space allocated for the STATIC_PARAMS Sigmas and Taus */
    /*	err = free_space_for_static_sigtau_for_calib(
                    eModelDim,
                            lNumSigmas,
                            lNumTaus);
                            if (err)
                            {
                            free_calib_core_memory;
                            return err;
                            }
    */
    /* Sucess: rebuilds and replaces the term structure in the underlying */
    err = from_optparam_to_sigtau(
        pdOptParams, /* From [1] to [param_number] */
        ppdSigmaValues,
        lNumSigmas,
        ppdTauValues,
        lNumTaus,
        psCalibParams->bFreezeTau,
        psCalibParams->bOneTau,
        (*dFixedTau),
        psCalibParams->bFreezeBeta,
        psCalibParams->bOneBeta,
        dFixedBeta,
        eModelDim,
        psCalibParams->bFreezeOmega,
        psCalibParams->bOneOmega,
        dFixedOmega,
        psCalibParams->eCalibType,
        dFixedAlpha,
        dFixedGamma,
        dFixedRho);
    if (err)
    {
        free_calib_core_memory;
        return err;
    }

    err = srt_f_init_IRM_TermStruct(
        lClcnDate,
        ppdSigmaValues,
        lNumSigmaCols,
        lNumSigmas,
        ppdTauValues,
        lNumTauCols,
        lNumTaus,
        eModelType,
        eModelDim,
        0.0, /* Beta */
        0.0, /* Alpha */
        0.0, /* Gamma */
        0.0, /* Rho */
        0.0, /* Vovol */
        0.0, /* Eta */
        0.0,
        0,
        0,
        NULL,
        &psTermStruct);
    if (err)
    {
        free_calib_core_memory;
        return err;
    }

    /* Frees the memory allocated locally */
    free_calib_core_memory;

    /* Free the memory of dFixedTau */
    if (*dFixedTau)
        (*dFixedTau) = NULL;
    if (dFixedTau)
        dFixedTau = NULL;

    /* Attaches the Term Struct to the underlying */
    set_irund_ts(sUndPtr, psTermStruct);

    smessage("Calibration achieved");
    /* Gets the Ending Time */
    lEndTime = time(NULL);

    smessage("Calibration achieved in %.2f s", difftime(lEndTime, lStartTime));
    smessage("");

    /* Convert pdChiSquare into bp */
    *pdChiSquare = 100 * sqrt(*pdChiSquare);
    return NULL;
}
/* ======================================================================= */

#ifdef PVMPI
/* ======================================================================= */
/*
 *set_parallel_context
 *
 * Author: Cyril Godart
 * Date: 14/08/98
 *
 * Entry point to this file for a distant process.
 * Parameters are passed through any
 * method : message passing, memory sharing...
 * The set_static funcs are called so  the "static" storage allocation
 * specifier has not to be removed. Overhead is not critical since
 * this operation is done once for each calibration.
 */

Err set_parallel_context(int bLocalComp, char* szGroupName, int iNumCPUs, int* piCPUs, int Version)
{
    if (szGroupName)
    {
        PRL_CONTEXT.szGroupName = (char*)calloc(strlen(szGroupName) + 1, sizeof(char));
        strcpy(PRL_CONTEXT.szGroupName, szGroupName);
    }
    else
        PRL_CONTEXT.szGroupName = (char*)NULL;
    PRL_CONTEXT.bLocalComp = bLocalComp;
    PRL_CONTEXT.iNumCPUs   = iNumCPUs;
    PRL_CONTEXT.piCPUs     = piCPUs;
    PRL_CONTEXT.bAggregate = 0;
    return NULL;
}

void set_msg_tag(int iMsgTag)
{
    PRL_CONTEXT.iMsgTag = iMsgTag;
}

Err free_calib_memory(SrtMdlDim eModelDim, long lNumSigmas, long lNumTaus)
{
    free_space_for_static_sigtau_for_calib(eModelDim, lNumSigmas, lNumTaus);
    return NULL;
}

void free_weights_targets_and_bounds()
{
    if (STATIC_PARAMS.pdMarketWeights)
        free_dvector(STATIC_PARAMS.pdMarketWeights, 1, STATIC_PARAMS.lNumData);
    STATIC_PARAMS.pdMarketWeights = NULL;
    if (STATIC_PARAMS.pdMarketTargets)
        free_dvector(STATIC_PARAMS.pdMarketTargets, 1, STATIC_PARAMS.lNumData);
    STATIC_PARAMS.pdMarketTargets = NULL;
    if (STATIC_PARAMS.ppdParamBounds)
        free_dmatrix(STATIC_PARAMS.ppdParamBounds, 1, 2, 1, STATIC_PARAMS.lNumParams);
    STATIC_PARAMS.ppdParamBounds = NULL;
}

void free_instruments()
{
    if (STATIC_INSTRUMENT.lNumInstruments > 0)
    {
        if (STATIC_INSTRUMENT.psSwapDp)
            srt_free(STATIC_INSTRUMENT.psSwapDp);
        STATIC_INSTRUMENT.psSwapDp = NULL;
        if (STATIC_INSTRUMENT.pdStrike)
            srt_free(STATIC_INSTRUMENT.pdStrike);
        STATIC_INSTRUMENT.pdStrike = NULL;
        if (STATIC_INSTRUMENT.pdBondStrike)
            srt_free(STATIC_INSTRUMENT.pdBondStrike);
        STATIC_INSTRUMENT.pdBondStrike = NULL;
        if (STATIC_INSTRUMENT.peRecPay)
            srt_free(STATIC_INSTRUMENT.peRecPay);
        STATIC_INSTRUMENT.peRecPay = NULL;
        if (STATIC_INSTRUMENT.peProductType)
            srt_free(STATIC_INSTRUMENT.peProductType);
        STATIC_INSTRUMENT.peProductType = NULL;
        if (STATIC_INSTRUMENT.pszRefRateCode)
            free_svector_size(
                STATIC_INSTRUMENT.pszRefRateCode, 0, STATIC_INSTRUMENT.lNumInstruments - 1, 32);
        STATIC_INSTRUMENT.pszRefRateCode = NULL;
    }
    if (STATIC_INSTRUMENT.psGrfnParams)
        free(STATIC_INSTRUMENT.psGrfnParams);
    STATIC_INSTRUMENT.psGrfnParams = NULL;
    if (STATIC_INSTRUMENT.pdFraMaturities)
        free_vector(STATIC_INSTRUMENT.pdFraMaturities, 0, STATIC_INSTRUMENT.lNumTenors - 1);
    STATIC_INSTRUMENT.pdFraMaturities = NULL;
}
#endif
/* ------------------------------------------------------------------------ */
