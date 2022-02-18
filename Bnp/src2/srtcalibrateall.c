/* =======================================================================
        FILENAME:    srtnewcalibrateall.c

    AUTHORS:     Olivier Van Eyseren, Antoine Savine and JLX

        PURPOSE:     Top level routine for a all calibrations
                     This function should be called from any interface

        FUNCTION:    Err SrtNewCalibrateAll(...)
  ========================================================================== */
#include "math.h"
#include "srt_h_all.h"
#ifdef PVMPI
#define ID_CALIB_FUNC 225
#include "parallel.h"
#endif

/* global pointer to the array of computed calibrated model prices
   modified in levenberg_calib_funcs */
double* GlobalTheoPrices;

#define global_calib_free                                                               \
    {                                                                                   \
        if (peOptionType)                                                               \
            srt_free(peOptionType);                                                     \
        peOptionType = NULL;                                                            \
        if (peRecPay)                                                                   \
            srt_free(peRecPay);                                                         \
        peRecPay = NULL;                                                                \
        if (psSwapDp)                                                                   \
            srt_free(psSwapDp);                                                         \
        psSwapDp = NULL;                                                                \
        if (plDealSigmaDates)                                                           \
            free_dvector(plDealSigmaDates, 0, *plNumSigmas - 1);                        \
        plDealSigmaDates = NULL;                                                        \
        if (*pppdSigmaCurve)                                                            \
            free_dmatrix(*pppdSigmaCurve, 0, *plNumSigmaCols - 1, 0, *plNumSigmas - 1); \
        *pppdSigmaCurve = NULL;                                                         \
        if (plDealTauDates)                                                             \
            free_dvector(plDealTauDates, 0, *plNumTaus - 1);                            \
        plDealTauDates = NULL;                                                          \
        if (*pppdTauCurve)                                                              \
            free_dmatrix(*pppdTauCurve, 0, *plNumTauCols - 1, 0, *plNumTaus - 1);       \
        *pppdTauCurve = NULL;                                                           \
        if (pdFraMaturities)                                                            \
            free_vector(pdFraMaturities, 0, lNumTenors - 1);                            \
        pdFraMaturities = NULL;                                                         \
        err             = free_MCrandSrc();                                             \
    }

/**************************************************
 * change_underlying_name.
 * Change: Cyril Godart 28/09/98
 *         Made function so that parallel version
 *         of the code can make use of it on slaves.
 ****************************************************/
Err change_und_name(char* the_und)
{
    Err            err = NULL;
    char           tmp_str[256];
    SrtUndPtr      sUndPtr;
    SrtUndListPtr  und_list = NULL;
    SrtListObject* obj      = NULL;
#ifdef PVMPI
    if (!PRL_CONTEXT.bLocalComp)
        RequestChangeUndName(the_und);
#endif
    /* Change underlying name if needed */
    if (!strstr(the_und, "_CAL"))
    {
        strcpy(tmp_str, the_und);
        strcat(tmp_str, "_cal");
        strupper(tmp_str);
        strip_white_space(tmp_str);

        while (lookup_und(tmp_str))
        {
            err = srt_f_destroy_und(tmp_str);
            if (err)
            {
                return err;
            }
        }

        sUndPtr = lookup_und(the_und);

        strcpy(the_und, tmp_str);

        /* Get the underlying list and check if not empty */
        und_list = get_underlying_list();
        if (und_list == NULL)
        {
            return serror("No Underlying list defined: call SrtInit before");
        }

        /* Get the underlying == SrtListObject in the underlying_list */
        err = srt_f_lstgetobj(*und_list, get_underlying_name(sUndPtr), 0.0, &obj);
        if (err)
        {
            return serror("Could not find %s underlying in list", get_underlying_name(sUndPtr));
        }

        /* Change the name of the underlying as underlying and as object */
        strcpy(get_underlying_name(sUndPtr), tmp_str);
        strcpy(obj->name, tmp_str);
    }
    return NULL;
}

/* ------------------------------------------------------------------------
        FUNCTION: SrtCalibrateAll
        PURPOSE: calibrate model acoording to user's inputs:
                if GlobalCalibration:
                        a correlation matrix of FRA is required
                        correlation structure is supposed to be time independent
                        optimisation is made on:
                                sig1(t), lam1(t)
                                alpha, beta, rho
                        the criteria is a global one: prices + correlation
                if BlindCalibration:
                        no correlation is required
                        optimisation is made blindly on:
                                sig1(t), lam1(t)
                                alpha(t), beta(t), rho(t)
                        the criteria is just the based on market prices
                if FixedCalibration:
                        the correlation structure (alpha, beta, rho) is frozen
                        optimisation is made only on sig1(t), lam1(t)
                        the criteria is just based on market prices
        USAGE: to get to a :
                Global Calibration: DEFAULT
                        calib_param->type set to GLOBAL
                        the FRD matrix has to be defined
                Blind Calibration:
                        calib_param->type set to BLIND
                        the FRD matrix is set to NULL
                Fixed Calibration
                        calib_param->type set to FIXED
                        the FRD matrix is set to NULL
 ------------------------------------------------------------------------- */

#define CALIBRATE_MAX_STDEV 3.0
#define global_new_calib_free                                                           \
    {                                                                                   \
        if (peOptionTypeFree)                                                           \
            srt_free(peOptionTypeFree);                                                 \
        peOptionTypeFree = NULL;                                                        \
        if (peRecPayFree)                                                               \
            srt_free(peRecPayFree);                                                     \
        peRecPayFree = NULL;                                                            \
        if (psSwapDpFree)                                                               \
            srt_free(psSwapDpFree);                                                     \
        psSwapDpFree = NULL;                                                            \
        if (plDealSigmaDates)                                                           \
            free_dvector(plDealSigmaDates, 0, *plNumSigmas - 1);                        \
        plDealSigmaDates = NULL;                                                        \
        if (*pppdSigmaCurve)                                                            \
            free_dmatrix(*pppdSigmaCurve, 0, *plNumSigmaCols - 1, 0, *plNumSigmas - 1); \
        *pppdSigmaCurve = NULL;                                                         \
        if (plDealTauDates)                                                             \
            free_dvector(plDealTauDates, 0, *plNumTaus - 1);                            \
        plDealTauDates = NULL;                                                          \
        if (*pppdTauCurve)                                                              \
            free_dmatrix(*pppdTauCurve, 0, *plNumTauCols - 1, 0, *plNumTaus - 1);       \
        *pppdTauCurve = NULL;                                                           \
        if (pdFraMaturities)                                                            \
            free_vector(pdFraMaturities, 0, lNumTenors - 1);                            \
        pdFraMaturities = NULL;                                                         \
        free_MCrandSrc();                                                               \
    }

/* ------------------------------------------------------------------------
        FUNCTION: SrtNewCalibrateAll
        PURPOSE: calibrate model acoording to user's inputs:
                if GlobalCalibration:
                        a correlation matrix of FRA is required
                        correlation structure is supposed to be time independent
                        optimisation is made on:
                                sig1(t), lam1(t)
                                alpha, beta, rho
                        the criteria is a global one: prices + correlation
                if BlindCalibration:
                        no correlation is required
                        optimisation is made blindly on:
                                sig1(t), lam1(t)
                                alpha(t), beta(t), rho(t)
                        the criteria is just the based on market prices
                if FixedCalibration:
                        the correlation structure (alpha, beta, rho) is frozen
                        optimisation is made only on sig1(t), lam1(t)
                        the criteria is just based on market prices
        USAGE: to get to a :
                Global Calibration: DEFAULT
                        calib_param->type set to GLOBAL
                        the FRD matrix has to be defined
                Blind Calibration:
                        calib_param->type set to BLIND
                        the FRD matrix is set to NULL
                Fixed Calibration
                        calib_param->type set to FIXED
                        the FRD matrix is set to NULL
 ------------------------------------------------------------------------- */

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
    int*      indexUsedInstr,
    int*      InstrsIndexArrays)
{
    Err err = NULL;
    int i;

    SrtUndPtr         sUndPtr;
    SrtUnderlyingType eUnderlyingType;
    SrtMdlType        eModelType;
    SrtMdlDim         eModelDim;
    SrtCalibParam     calib_param;
    SrtGrfnParam      sGrfnParams;
    String            szYieldCurveName;
    SwapDP*           psSwapDp     = NULL;
    SwapDP*           psSwapDpFree = NULL;
    SrtCcyParam*      psCcyParams  = NULL;
    long              today;
    StructType*       peOptionType     = NULL;
    StructType*       peOptionTypeFree = NULL;
    SrtReceiverType*  peRecPay         = NULL;
    SrtReceiverType*  peRecPayFree     = NULL;
    double*           plDealSigmaDates = NULL;
    double*           plDealTauDates   = NULL;

    double*        pdFraMaturities = NULL;
    TermStruct*    ts              = NULL;
    char           the_und[256], tmp_str[256];
    SrtCurvePtr    yldcrv;
    SrtListObject* obj      = NULL;
    SrtUndListPtr  und_list = NULL;

    SrtDateList      date_list;
    long             last_period_start;
    long             first_period_start;
    long             lLastFixing;
    long             lFirstFixing;
    int              Spot_Lag;
    double           dForwardRate;
    double           dATMVol;
    double           dSpread;
    double           dForwardCash;
    double           dfStart, dfEnd;
    double           dLevel;
    double           dDelta;
    SrtDiffusionType eLogOrNorm;
    double           dStdDev;
    double           dDistance;
    double           dMaturity;
    int              iNumUsedInstruments;
    int              iNumTrueInstruments;
    int              j;

    /* global pointer to the array of computed calibrated model prices */
    GlobalTheoPrices = (*ppdTheoPrices);

    /* Gets the underlying */
    strcpy(the_und, szUndName);
    sUndPtr = lookup_und(the_und);
    if (!sUndPtr)
    {
        /* has not found the underlying: tries with _CAL */
        strcat(the_und, "_cal");
        strupper(the_und);
        strip_white_space(the_und);
        sUndPtr = lookup_und(the_und);
        if (!sUndPtr)
        {
            /* definitely, underlying does not exist */
            return serror("Unknown underlying");
        }
    }

    /* Checks the underlying is an interest rate */
    eUnderlyingType = get_underlying_type(sUndPtr);
    if (!ISUNDTYPE(sUndPtr, INTEREST_RATE_UND))
    {
        return serror("Cannot calibrate a non IR underlying");
    }

    /* Checks the model is compatible for calibration (LGM or CHEY) */
    if (err = get_underlying_mdltype(sUndPtr, &eModelType))
    {
        return err;
    }

    if (is_model_Cheyette_type(eModelType) == SRT_NO)
    {
        return serror("Can only calibrate LGM, CHEY, CHEYBETA (with STOCHVOL) so far");
    }

    /* Check the Term Structure depending on the model dimension */
    if (err = get_underlying_mdldim(sUndPtr, &eModelDim))
    {
        return err;
    }

    if (eModelDim == TWO_FAC)
    {
        if ((*plNumTauCols != 0) && (*plNumTauCols != 3))
        {
            return serror("Need a FULL Two Factor Tau TS for Calibration: Date - Tau1 - Tau2");
        }
        if ((*plNumSigmaCols != 0) && (*plNumSigmaCols != 6))
        {
            return serror(
                "Need a FULL Two Factor Sigma TS for Calibration: Date - Sig1 - Bet1 - Sig2 - Bet2 "
                "- Rho");
        }
    }
    else if (eModelDim == ONE_FAC)
    {
        if ((*plNumTauCols != 0) && (*plNumTauCols != 2))
        {
            return serror("Need a FULL One Factor Tau TS for Calibration : Date - Tau");
        }
        if ((*plNumSigmaCols != 0) && (*plNumSigmaCols != 3))
        {
            return serror("Need a FULL One Factor Sigma TS for Calibration : Date - Sig - Beta");
        }
    }
    else
    {
        return serror("Model should be 1D or 2D");
    }

    /* Sets the calibration parameters : default or user defined */
    err = srt_f_set_CalibParams(
        pszCalibParamNames, pszCalibParamValues, lNumCalibParams, &calib_param, eModelType);
    if (err)
    {
        return err;
    }

    /* Sets the Grfn Parameters  */
    err = srt_f_set_GrfnParams(lNumGrfnParams, pszGrfnParamNames, pszGrfnParamValues, &sGrfnParams);
    if (err)
    {
        return err;
    }

    /* New flag for the calib */
    if (calib_param.bAggregate == SRT_YES)
        sGrfnParams.calib = SRT_YES;

    /* Gets today from the underlying */
    szYieldCurveName = get_discname_from_underlying(sUndPtr);
    yldcrv           = lookup_curve(szYieldCurveName);
    psCcyParams      = get_ccyparam_from_yldcrv(yldcrv);
    today            = (Date)get_clcndate_from_curve(yldcrv);

    /* interpret strings */
    if (err = interp_diffusion_type(*LogNormStr, &eLogOrNorm))
    {
        return serror("Could not interpret the Type of the Vol in Market Vol Matrix");
    }

    /* Allocate an array of SwapDP and initialise them: start, end, comp, basis */
    psSwapDp     = (SwapDP*)srt_calloc(lNumInstruments, sizeof(SwapDP));
    psSwapDpFree = psSwapDp;
    if (!psSwapDp)
    {
        return serror("Memory allocation error (psSwapDp) in srt_f_CalibrateAll");
    }
    peOptionType     = (int*)srt_calloc(lNumInstruments, sizeof(int));
    peOptionTypeFree = peOptionType;
    if (!peOptionType)
    {
        srt_free(psSwapDpFree);
        return serror("Memory allocation error (peOptionType) in srt_f_CalibrateAll");
    }
    peRecPay     = (int*)srt_calloc(lNumInstruments, sizeof(int));
    peRecPayFree = peRecPay;
    if (!peRecPay)
    {
        srt_free(psSwapDpFree);
        srt_free(peOptionTypeFree);
        return serror("Memory allocation error (peRecPay) in srt_f_CalibrateAll");
    }

    for (i = 0; i < lNumInstruments; i++)
    {
        err = swp_f_initSwapDP(
            plStartDates[i], plEndDatesOrNfp[i], pszFrequency[i], pszBasis[i], &psSwapDp[i]);
        if (err)
        {
            global_new_calib_free;
            return err;
        }
        psSwapDp[i].spot_lag = psCcyParams->spot_lag;
    }

    /* Allocate a vector of StructType, SrtReceiverType and initialises them */
    for (i = 0; i < lNumInstruments; i++)
    {
        err = interp_rec_pay(pszRecPay[i], &peRecPay[i]);
        if (err)
        {
            global_new_calib_free;
            return (err);
        }

        err = interp_struct(pszOptionType[i], &peOptionType[i]);
        if (err)
        {
            global_new_calib_free;
            return (err);
        }
    } /* END loop on instruments */

    /****************************/
    /* Beginning of the CHANGES */
    /****************************/

    /* Sets an initial guess only for instruments with a fixing date after today */
    iNumUsedInstruments = 0;
    iNumTrueInstruments = lNumInstruments;

    for (i = 0; i < lNumInstruments; i++)
    {
        Spot_Lag = psSwapDp[iNumUsedInstruments].spot_lag;
        /* Builds the schedule of the Cap/Floor */
        date_list = SwapDP_to_DateList(&(psSwapDp[iNumUsedInstruments]), MODIFIED_SUCCEEDING);
        /* Extracts the relevant fixing date for the deal */

        if (peOptionType[iNumUsedInstruments] == CAPFLOOR)
        {
            /* Gets the last period start: not last date, the one before */
            /* Gets the first period start */
            first_period_start = psSwapDp[iNumUsedInstruments].start;
            lFirstFixing = add_unit(first_period_start, -Spot_Lag, SRT_BDAY, MODIFIED_SUCCEEDING);
            /* Calibration on the last caplet of the cap */
            last_period_start = date_list.date[date_list.len - 2];
            lLastFixing = add_unit(last_period_start, -Spot_Lag, SRT_BDAY, MODIFIED_SUCCEEDING);
        }
        else if (
            (peOptionType[iNumUsedInstruments] == SWAPTION) ||
            (peOptionType[iNumUsedInstruments] == BOND_OPTION)
            ///// added by Albert Wang 08/25/03 - begin
            || (peOptionType[iNumUsedInstruments] == SIMPLEMIDAT)
            ///// added by Albert Wang 08/25/03 - end
        )
        {
            /* The Sigma date is the fixng date of the swaption */
            lLastFixing = add_unit(
                psSwapDp[iNumUsedInstruments].start, -Spot_Lag, SRT_BDAY, MODIFIED_SUCCEEDING);
        }
        /* We only keep track of the instruments with a fixing date before or equal to today */
        if (lLastFixing > today)
        {
            /* In the case of the full Cap and if the first fixing is before today we modify the
             * start date and the relevant swapDP */
            if (peOptionType[iNumUsedInstruments] == CAPFLOOR)
            {
                int indexcap = 0;
                while ((lFirstFixing <= today) && (indexcap < date_list.len - 2))
                {
                    indexcap++;
                    first_period_start = date_list.date[indexcap];
                    lFirstFixing =
                        add_unit(first_period_start, -Spot_Lag, SRT_BDAY, MODIFIED_SUCCEEDING);
                }

                plStartDates[iNumUsedInstruments] = first_period_start;

                err = swp_f_initSwapDP(
                    first_period_start,
                    plEndDatesOrNfp[iNumUsedInstruments],
                    pszFrequency[iNumUsedInstruments],
                    pszBasis[iNumUsedInstruments],
                    &psSwapDp[iNumUsedInstruments]);
                if (err)
                {
                    srt_free(date_list.date);
                    global_new_calib_free;
                    return err;
                }
                psSwapDp[iNumUsedInstruments].spot_lag = psCcyParams->spot_lag;
                err                                    = swp_f_CapFloor_SwapDP(
                    &psSwapDp[iNumUsedInstruments],
                    pdStrike[iNumUsedInstruments],
                    pfGetVol,
                    peRecPay[iNumUsedInstruments],
                    pszRefRateCode[iNumUsedInstruments],
                    szYieldCurveName,
                    PREMIUM,
                    eLogOrNorm,
                    &pdPrice[iNumUsedInstruments]);
                if (err)
                {
                    srt_free(date_list.date);
                    global_new_calib_free;
                    return err;
                }

                err = swp_f_CapFloor_SwapDP(
                    &psSwapDp[iNumUsedInstruments],
                    pdStrike[iNumUsedInstruments],
                    pfGetVol,
                    peRecPay[iNumUsedInstruments],
                    pszRefRateCode[iNumUsedInstruments],
                    szYieldCurveName,
                    DELTA,
                    eLogOrNorm,
                    &dDelta);
                if (err)
                {
                    srt_free(date_list.date);
                    global_new_calib_free;
                    return err;
                }

                err =
                    swp_f_Level_SwapDP(&(psSwapDp[iNumUsedInstruments]), szYieldCurveName, &dLevel);

                dDistance = fabs(inv_cumnorm_fast(fabs(dDelta / dLevel)));
            }

            /* SWAPTION */
            if ((peOptionType[iNumUsedInstruments] == SWAPTION)
                ///// added by Albert Wang 08/25/03 - begin
                || (peOptionType[iNumUsedInstruments] == SIMPLEMIDAT)
                ///// added by Albert Wang 08/25/03 - end
            )
            {
                err = swp_f_ForwardRate_SwapDP(
                    &(psSwapDp[iNumUsedInstruments]),
                    szYieldCurveName,
                    pszRefRateCode[iNumUsedInstruments],
                    &dForwardRate);
                if (err)
                {
                    srt_free(date_list.date);
                    global_new_calib_free;
                    return err;
                }
                /* compute the spread */
                err =
                    swp_f_Level_SwapDP(&(psSwapDp[iNumUsedInstruments]), szYieldCurveName, &dLevel);
                dfStart      = swp_f_df(today, date_list.date[0], szYieldCurveName);
                dfEnd        = swp_f_df(today, date_list.date[date_list.len - 1], szYieldCurveName);
                dForwardCash = (dfStart - dfEnd) / dLevel;
                dSpread      = dForwardRate - dForwardCash;
                /*gets the ATM Vol */
                err = pfGetVol(
                    (double)DTOL(date_list.date[0]),
                    (double)DTOL(date_list.date[date_list.len - 1]),
                    dForwardRate,
                    dForwardRate,
                    dSpread,
                    &dATMVol);
                if (err)
                {
                    srt_free(date_list.date);
                    global_new_calib_free;
                    return err;
                }
                /* The maturity of the Option */
                dMaturity = (lLastFixing - today) * YEARS_IN_DAY;

                /* Computes the Distance betwenn Forward and Strike in Num Std Dev */
                dStdDev = sqrt(dMaturity) * dATMVol;

                if (eLogOrNorm == SRT_NORMAL)
                {
                    dDistance = fabs(dForwardRate - pdStrike[iNumUsedInstruments]) / dStdDev;
                }
                else
                {
                    dDistance = fabs(log(dForwardRate / pdStrike[iNumUsedInstruments])) / dStdDev;
                }
            }

            if (date_list.date)
                srt_free(date_list.date);

            /* If the Strike is more than MAX_STDEV away from the Forward, discard the instrument */
            if (dDistance > CALIBRATE_MAX_STDEV)
            {
                /* Removes this instrument from the list */
                for (j = iNumUsedInstruments; j < iNumTrueInstruments - 1; j++)
                {
                    psSwapDp[j]     = psSwapDp[j + 1];
                    pdStrike[j]     = pdStrike[j + 1];
                    pdBondStrike[j] = pdBondStrike[j + 1];
                    peOptionType[j] = peOptionType[j + 1];
                    peRecPay[j]     = peRecPay[j + 1];
                    pdPrice[j]      = pdPrice[j + 1];
                    pdVega[j]       = pdVega[j + 1];
                    strcpy(pszRefRateCode[j], pszRefRateCode[j + 1]);

                    /* For the output purpose we need to keep track of the details of the
                     * instruments in input */
                    plStartDates[j]    = plStartDates[j + 1];
                    plEndDatesOrNfp[j] = plEndDatesOrNfp[j + 1];
                    strcpy(pszFrequency[j], pszFrequency[j + 1]);
                    strcpy(pszBasis[j], pszBasis[j + 1]);
                    pdBondStrike[j] = pdBondStrike[j + 1];
                    strcpy(pszRecPay[j], pszRecPay[j + 1]);
                    strcpy(pszOptionType[j], pszOptionType[j + 1]);
                }
                continue;
            }
            else
            { /* instrument is accepted for calibration */
                indexUsedInstr[iNumUsedInstruments] = i;
            }

            /* Increment the Number of instuments by one */
            iNumUsedInstruments++;

        } /* END if (lLastFixing >clcndate ) */
        else
        {
            /* Removes this instrument from the list */
            psSwapDp++;
            pdStrike++;
            pdBondStrike++;
            peOptionType++;
            peRecPay++;
            pdPrice++;
            pdVega++;
            pszRefRateCode++;

            iNumTrueInstruments--;

            plStartDates++;
            plEndDatesOrNfp++;
            pszFrequency++;
            pszBasis++;
            pszRecPay++;
            pszOptionType++;
        }
    }

    /****************************/
    /* End of CHANGES           */
    /****************************/

    /* Gets the initialised underlying term struct */
    err = get_underlying_ts(sUndPtr, &ts);
    if (err)
    {
        global_new_calib_free;
        return err;
    }

    /* Checks if the Vol and Tau dates are user defined */
    if ((!(*pppdSigmaCurve) || (*plNumSigmas == 0)) && (!(*pppdTauCurve) || (*plNumTaus == 0)))
    {
        /* If only one TAU required, sets a default date */
        if (calib_param.bOneTau == SRT_YES)
        {
            *plNumTaus        = 1;
            plDealTauDates    = dvector(0, *plNumTaus - 1);
            plDealTauDates[0] = (double)today + 365.0;

            /* Allocation will be done inside with new plNumSigmas */
            err = find_sigdates_from_deals(
                psSwapDp,
                peOptionType,
                iNumUsedInstruments,
                psCcyParams,
                &plDealSigmaDates,
                plNumSigmas);
            if (err)
            {
                global_new_calib_free;
                return (err);
            }
        }
        else if (calib_param.bFreezeTau == SRT_YES)
        {
            if (eModelDim == ONE_FAC)
            {
                double* plTemSigmaDates = NULL;
                double* pdTempSigma     = NULL;
                double* pdTempBeta      = NULL;
                double* pdTempVovol     = NULL;
                double* pdTempRho       = NULL;
                double* pdTempMeanVol   = NULL;

                long    plTempNumSigmas;
                double* pdTempTau      = NULL;
                double* plTempTauDates = NULL;

                err = srt_f_display_IRM_OneFac_TermStruct(
                    ts,
                    &plTemSigmaDates,
                    &pdTempSigma,
                    &pdTempBeta,
                    &pdTempVovol,
                    &pdTempRho,
                    &pdTempMeanVol,
                    &plTempNumSigmas,
                    &plTempTauDates,
                    &pdTempTau,
                    plNumTaus);

                if (plTemSigmaDates)
                    free(plTemSigmaDates);
                if (pdTempSigma)
                    free(pdTempSigma);
                if (pdTempBeta)
                    free(pdTempBeta);
                if (pdTempVovol)
                    free(pdTempVovol);
                if (pdTempRho)
                    free(pdTempRho);
                if (pdTempMeanVol)
                    free(pdTempMeanVol);
                if (pdTempTau)
                    free(pdTempTau);

                plTemSigmaDates = NULL;
                pdTempSigma     = NULL;
                pdTempBeta      = NULL;
                pdTempTau       = NULL;
                pdTempVovol     = NULL;
                pdTempRho       = NULL;
                pdTempMeanVol   = NULL;

                plDealTauDates = dvector(0, *plNumTaus - 1);
                for (i = 0; i < *plNumTaus; i++)
                    plDealTauDates[i] = plTempTauDates[i];

                free(plTempTauDates);
                plTempTauDates = NULL;
            }
            else if (eModelDim == TWO_FAC)
            {
                double* plTemSigmaDates;
                double* pdTempSigma1;
                double* pdTempSigma2;
                double* pdTempBeta1;
                double* pdTempBeta2;
                double* pdTempRho;
                long    plTempNumSigmas;
                double* pdTempTau1;
                double* pdTempTau2;
                double* plTempTauDates;

                err = srt_f_display_IRM_TwoFac_TermStruct(
                    ts,
                    &plTemSigmaDates,
                    &pdTempSigma1,
                    &pdTempBeta1,
                    &pdTempSigma2,
                    &pdTempBeta2,
                    &pdTempRho,
                    &plTempNumSigmas,
                    &plTempTauDates,
                    &pdTempTau1,
                    &pdTempTau2,
                    plNumTaus);

                free(plTemSigmaDates);
                free(pdTempSigma1);
                free(pdTempSigma2);
                free(pdTempRho);
                free(pdTempBeta1);
                free(pdTempBeta2);
                free(pdTempTau1);
                free(pdTempTau2);
                plTemSigmaDates = NULL;
                pdTempSigma1    = NULL;
                pdTempSigma2    = NULL;
                pdTempRho       = NULL;
                pdTempBeta1     = NULL;
                pdTempBeta2     = NULL;
                pdTempTau1      = NULL;
                pdTempTau2      = NULL;

                plDealTauDates = dvector(0, *plNumTaus - 1);
                for (i = 0; i < *plNumTaus; i++)
                    plDealTauDates[i] = plTempTauDates[i];

                free(plTempTauDates);
                plTempTauDates = NULL;
            }
            else
            {
                global_new_calib_free;
                return (err);
            }
            /* Allocation will be done inside with new plNumSigmas */
            err = find_sigdates_from_deals(
                psSwapDp,
                peOptionType,
                iNumUsedInstruments,
                psCcyParams,
                &plDealSigmaDates,
                plNumSigmas);
            if (err)
            {
                global_new_calib_free;
                return (err);
            }
        }
        else
        {
            /* Allocation will be done inside with new *plNumTaus */
            err = find_dates_from_deals(
                psSwapDp,
                peOptionType,
                iNumUsedInstruments,
                psCcyParams,
                &plDealSigmaDates,
                plNumSigmas,
                &plDealTauDates,
                plNumTaus);
            if (err)
            {
                global_new_calib_free;
                return (err);
            }
        } /* END if (one_tau == SRT_NO) */
    }
    else if (!(*pppdSigmaCurve) || (*plNumSigmas == 0))
    {
        /* Allocation will be done inside with new plNumSigmas */
        err = find_sigdates_from_deals(
            psSwapDp,
            peOptionType,
            iNumUsedInstruments,
            psCcyParams,
            &plDealSigmaDates,
            plNumSigmas);
        if (err)
        {
            global_new_calib_free;
            return (err);
        }
    }
    else if (!(*pppdTauCurve) || (*plNumTaus == 0))
    {
        /* Allocation will be done inside with new plNumSigmas */
        err = find_taudates_from_deals(
            psSwapDp, iNumUsedInstruments, psCcyParams, &plDealTauDates, plNumTaus);
        if (err)
        {
            global_new_calib_free;
            return (err);
        }
    }

    /* Checks if the vol dates are user defined or compute them */
    if (!(*pppdSigmaCurve) || (*plNumSigmas == 0))
    {
        /* Allocate *pppdSigmaCurve and sets dates and values for sig from TS */
        if (eModelDim == TWO_FAC)
        {
            /* Allocate memory for a FULL Sigma Curve : Dates - Sig1 - Beta1 - Sig2 - Beta2- Rho */
            *plNumSigmaCols = 6;
            *pppdSigmaCurve = dmatrix(0, *plNumSigmaCols - 1, 0, *plNumSigmas - 1);
            if (!(*pppdSigmaCurve))
            {
                global_new_calib_free;
                return serror("Memory allocation error (*pppdSigmaCurve) in srt_f_CalibrateAll");
            }

            /* Populates the Full Sigma Curve according to the Initial TS */
            for (i = 0; i < *plNumSigmas; i++)
            {
                /*Sets the Dates according to the deal */
                (*pppdSigmaCurve)[0][i] = plDealSigmaDates[i];

                /* Set values for Sigma1 and Sigma 2 */
                err = find_2f_sig(
                    (plDealSigmaDates[i] - today) * YEARS_IN_DAY,
                    ts,
                    &(*pppdSigmaCurve)[1][i],
                    &(*pppdSigmaCurve)[3][i]);
                if (err)
                {
                    global_new_calib_free;
                    return (err);
                }

                /* Set values for Beta1 and Beta2 */
                err = find_tf_beta(
                    (plDealSigmaDates[i] - today) * YEARS_IN_DAY,
                    ts,
                    &(*pppdSigmaCurve)[2][i],
                    &(*pppdSigmaCurve)[4][i]);
                if (err)
                {
                    global_new_calib_free;
                    return (err);
                }

                /* Set value for Rho */
                err = find_2f_rho(
                    (plDealSigmaDates[i] - today) * YEARS_IN_DAY, ts, &(*pppdSigmaCurve)[5][i]);
                if (err)
                {
                    global_new_calib_free;
                    return (err);
                }
            } /* END for i loop on number of Sigma */

        } /* END if eModelDim == TWO_FAC */

        else if (eModelDim == ONE_FAC)
        {
            /* Allocate memory for a FULL Sigma Curve : Dates - Sig - Beta */
            *plNumSigmaCols = 3;
            *pppdSigmaCurve = dmatrix(0, *plNumSigmaCols - 1, 0, *plNumSigmas - 1);

            if (!(*pppdSigmaCurve))
            {
                global_new_calib_free;
                return serror("Memory allocation error (*pppdSigmaCurve) in srt_f_CalibrateAll");
            }

            /* Populates the Full Sigma Curve according to the Initial TS */
            for (i = 0; i < *plNumSigmas; i++)
            {
                /*Sets the Dates according to the deal */
                (*pppdSigmaCurve)[0][i] = plDealSigmaDates[i];

                /* Set values for Sigma */
                (*pppdSigmaCurve)[1][i] =
                    find_sig((plDealSigmaDates[i] - today) * YEARS_IN_DAY, ts);

                /* Set values for Beta  */
                (*pppdSigmaCurve)[2][i] =
                    find_beta((plDealSigmaDates[i] - today) * YEARS_IN_DAY, ts);

            } /* END for i loop on number of sigmas */

        } /* END if (eModelDim == ONE_FAC) */

    } /* END if ( !(*pppdSigmaCurve) ) */

    /* Checks if the tau dates are user defined or compute them */
    if (!(*pppdTauCurve) || (*plNumTaus == 0))
    {
        /* Allocate *pppdTauCurve and sets dates and values for tau */
        if (eModelDim == TWO_FAC)
        {
            /* Allocate memory for a FULL Tau Curve : Dates - Tau1 - Tau2 */
            *plNumTauCols = 3;
            *pppdTauCurve = dmatrix(0, *plNumTauCols - 1, 0, *plNumTaus - 1);
            if (!(*pppdTauCurve))
            {
                global_new_calib_free;
                return serror("Memory allocation error (*pppdTauCurve) in srt_f_CalibrateAll");
            }
            /* Populates the Full Tau Curve according to the Initial TS */
            for (i = 0; i < *plNumTaus; i++)
            {
                /*Sets the Dates according to the deals (or OneTau)  */
                (*pppdTauCurve)[0][i] = plDealTauDates[i];

                /* Set values for Tau1 and Tau2 */
                err = find_2f_tau(
                    (plDealTauDates[i] - today) * YEARS_IN_DAY,
                    ts,
                    &(*pppdTauCurve)[1][i],
                    &(*pppdTauCurve)[2][i]);
                if (err)
                {
                    global_new_calib_free;
                    return (err);
                }

            } /* END for i loop on Num Taus */

        } /* END if (eModelDim == TWO_FAC) */

        else if (eModelDim == ONE_FAC)
        {
            /* Allocate memory for a FULL Tau Curve : Dates - Tau */
            *plNumTauCols = 2;
            *pppdTauCurve = dmatrix(0, *plNumTauCols - 1, 0, *plNumTaus - 1);
            if (!(*pppdTauCurve))
            {
                global_new_calib_free;
                return serror("Memory allocation error (*pppdTauCurve) in srt_f_CalibrateAll");
            }
            /* Populates the Full Tau Curve according to the Initial TS */
            for (i = 0; i < *plNumTaus; i++)
            {
                /*Sets the Dates according to the deals (or OneTau)  */
                (*pppdTauCurve)[0][i] = plDealTauDates[i];

                /* Set values for Tau */
                (*pppdTauCurve)[1][i] = find_tau((plDealTauDates[i] - today) * YEARS_IN_DAY, ts);

            } /* END for i loop on Num Taus */

        } /* END if (eModelDim == ONE_FAC) */

    } /* END if ( !(*pppdTauCurve) ) */

    /* If options weight = 100%, set lNumTenors to 0 */
    if (eModelDim == ONE_FAC)
    {
        calib_param.dOptionsWeight = 1.0;
        lNumTenors                 = 0;
    }
    else
    {
        /* If calib is not global, set num_tenors to 0 and frd weight to 100% */
        if (calib_param.eCalibType != GLOBAL_CALIB)
        {
            calib_param.dOptionsWeight = 1.00;
            lNumTenors                 = 0;
        }
        else
        {
            /* End if (calib_param.type != GLOBAL_CALIB)) */
            if (calib_param.dOptionsWeight > 0.99999999)
            {
                lNumTenors                 = 0;
                calib_param.dOptionsWeight = 1.00;
            }
            else if (calib_param.dOptionsWeight < 0.00000001)
            {
                iNumUsedInstruments        = 0;
                calib_param.dOptionsWeight = 0.00;
            }
            else
            {
                calib_param.dOptionsWeight = sqrt(sin(0.50 * SRT_PI * calib_param.dOptionsWeight));
            } /* End if (calib_param.frd_weight > 0.99999999) */

        } /* END if (calib_param.eCalibType == GLOBAL_CALIB) */

    } /* END if eModelDim == TWO_FAC */

    /* Converts the FRA tenors into length of time from today to fixing */
    if ((lNumTenors > 0) && (pszFraTenors != NULL))
    {
        pdFraMaturities = dvector(0, lNumTenors - 1);
        for (i = 0; i < lNumTenors; i++)
        {
            err = interp_corr_tenor(pszFraTenors[i], &pdFraMaturities[i]);
            if (err)
            {
                global_new_calib_free;
                return err;
            }
        }
    }
    else
    {
        calib_param.dOptionsWeight = 1.00;
        lNumTenors                 = 0;
    } /* END if (lNumTenors >0) */

    if (lNumTenors + iNumUsedInstruments < 1)
    {
        global_new_calib_free;
        return serror("Nothing to calibrate !!!");
    }

    /* Stochastic Volatility Models */
    if (eModelType == LGM_STOCH_VOL || eModelType == CHEY_STOCH_VOL ||
        eModelType == CHEY_BETA_STOCH_VOL)
    {
        /* Main call to stochastic volatility calibration routine */
        err = srt_f_calib_stoch_vol(
            &sGrfnParams,
            eModelType,
            eModelDim,
            psSwapDp,
            pdStrike,
            pdBondStrike,
            peOptionType,
            peRecPay,
            pdPrice,
            pdVega,
            pszRefRateCode,
            iNumUsedInstruments,
            pdFraMaturities,
            pdCorrelationMatrix,
            lNumTenors,
            sUndPtr,
            &calib_param,
            *pppdSigmaCurve,
            *plNumSigmas,
            *plNumSigmaCols,
            *pppdTauCurve,
            *plNumTaus,
            *plNumTauCols,
            pdChiSquare);
        if (err)
        {
            global_new_calib_free;
            return err;
        }
        if ((err = change_und_name(the_und)))
        {
            global_new_calib_free;
            return err;
        }
    }
    else
    /* Call to the calibration routine for non stochastic volatility ones */
    {
        /* Calibration of the real model */
        err = srt_f_calib_main(
            &sGrfnParams,
            eModelType,
            eModelDim,
            psSwapDp,
            pdStrike,
            pdBondStrike,
            peOptionType,
            peRecPay,
            pdPrice,
            pdVega,
            pszRefRateCode,
            iNumUsedInstruments,
            pdFraMaturities,
            pdCorrelationMatrix,
            lNumTenors,
            sUndPtr,
            &calib_param,
            *pppdSigmaCurve,
            *plNumSigmas,
            *plNumSigmaCols,
            *pppdTauCurve,
            *plNumTaus,
            *plNumTauCols,
            pdChiSquare);
        if (err)
        {
            global_new_calib_free;
            return err;
        }

    } /* END 	if (eModelType != LGM_STOCH_VOL && eModelType != CHEY_STOCH_VOL
                    &&  eModelType != CHEY_BETA_STOCH_VOL ) */

    /* Change underlying name if needed */
    if (!strstr(the_und, "_CAL"))
    {
        strcpy(tmp_str, the_und);
        strcat(tmp_str, "_cal");
        strupper(tmp_str);
        strip_white_space(tmp_str);

        while (lookup_und(tmp_str))
        {
            err = srt_f_destroy_und(tmp_str);
            if (err)
            {
                global_new_calib_free;
                return err;
            }
        }

        sUndPtr = lookup_und(the_und);

        strcpy(the_und, tmp_str);

        /* Get the underlying list and check if not empty */
        und_list = get_underlying_list();
        if (und_list == NULL)
        {
            global_new_calib_free;
            return serror("No Underlying list defined: call SrtInit before");
        }

        /* Get the underlying == SrtListObject in the underlying_list */
        err = srt_f_lstgetobj(*und_list, get_underlying_name(sUndPtr), 0.0, &obj);
        if (err)
        {
            global_new_calib_free;
            return serror("Could not find %s underlying in list", get_underlying_name(sUndPtr));
        }

        /* Change the name of the underlying as underlying and as object */
        strcpy(get_underlying_name(sUndPtr), tmp_str);
        strcpy(obj->name, tmp_str);
    }

    /* Free memory */
    global_new_calib_free

        /* Return a success message */
        strcpy(szNewName, the_und);
    return NULL;

} /* END of srt_f_NewCalibrateAll */
