/**********************************************************************
 *      Name: SrtBootstrap.c                                          *
 *  Function: Performs bootstrap calibration                          *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 01/10/95                                                *
 *--------------------------------------------------------------------*
 *    Inputs:                                                         *
 *   Returns:                                                         *
 *   Globals:                                                         *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 01/10/95 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/
#include "SrtAccess.h"
#include "grf_h_public.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_calib.h"
#include "srt_h_grfclsdfrm.h"

#define FREE_WRAP_BOOTSTRAP_MEMORY                 \
    if (sdpArray)                                  \
        srt_free(sdpArray);                        \
    sdpArray = NULL;                               \
    if (ATMprice)                                  \
        free_dvector(ATMprice, 0, (*nInputs) - 1); \
    ATMprice = NULL;                               \
    if (rec_pay)                                   \
        srt_free(rec_pay);                         \
    rec_pay = NULL;                                \
    if (derivType)                                 \
        srt_free(derivType);                       \
    derivType = NULL;

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
    char*   undName)
{
    Err              err    = NULL;
    int              status = 0;
    SrtUndPtr        und;
    SrtGrfnParam     grfnparam;
    int              i;
    SwapDP*          sdpArray  = NULL;
    int*             derivType = NULL;
    SrtReceiverType* rec_pay   = NULL;
    Date             clcndate;
    String           szYieldCurveName;
    SrtCurvePtr      yldcrv;
    SrtDateList      date_list;
    long             last_period_start;
    long             lLastFixing;
    int              Spot_Lag;
    double*          ATMprice = NULL;
    double           dForwardRate;
    double           dATMVol;
    double           dSpread;
    double           dForwardCash;
    double           dfStart, dfEnd;
    double           dLevel;
    SrtDiffusionType eLogOrNorm;
    long             dStart, dEnd;
    /* Overwrite defaults with user defined parameters */

    if (err = srt_f_set_GrfnParams(numParams, paramStrings, valueStrings, &grfnparam))
    {
        return err;
    }

    if (undName)
    {
        und = lookup_und(undName);
    }
    else
    {
        return serror("No Underlying passed to SrtCalboot");
    }

    if (und)
    {
        if (!ISUNDTYPE(und, INTEREST_RATE_UND))
        {
            return serror("Underlying must be of type Interest Rate");
        }
    }
    else
    {
        return serror("Could not find underlying in market list");
    }
    /* Compute the Clean Date to compare with the Fixing Date of the instruments */
    clcndate         = get_today_from_underlying(und);
    szYieldCurveName = get_ycname_from_irund(und);
    yldcrv           = lookup_curve(szYieldCurveName);

    /* interpret strings */
    if (err = interp_diffusion_type(*LogNormStr, &eLogOrNorm))
    {
        return serror("Could not interpret the Type of the Vol in Market Vol Matrix");
    }
    rec_pay   = (int*)calloc(*nInputs, sizeof(int));
    derivType = (int*)calloc(*nInputs, sizeof(int));
    ATMprice  = dvector(0, (*nInputs) - 1);

    if (!rec_pay || !derivType || !ATMprice)
    {
        return serror("Memory allocation failure in Bootstrap");
    }

    for (i = 0; i < *nInputs; i++)
    {
        err = interp_rec_pay(recPayStr[i], &rec_pay[i]);
        if (err)
        {
            FREE_WRAP_BOOTSTRAP_MEMORY;
            return err;
        }
        err = interp_struct(typeStr[i], &derivType[i]);
        if (err)
        {
            FREE_WRAP_BOOTSTRAP_MEMORY;
            return err;
        }
    }

    sdpArray = (SwapDP*)calloc(*nInputs, sizeof(SwapDP));

    for (i = 0; i < *nInputs; i++)
    {
        err = swp_f_initSwapDP(startDates[i], endDates[i], cpdStr[i], basisStr[i], &sdpArray[i]);
        if (err)
        {
            FREE_WRAP_BOOTSTRAP_MEMORY;
            return err;
        }
        sdpArray[i].spot_lag = get_spotlag_from_underlying(und);
        Spot_Lag             = sdpArray[i].spot_lag;
        /* Builds the schedule of the Cap/Floor */
        date_list = SwapDP_to_DateList(&(sdpArray[i]), MODIFIED_SUCCEEDING);
        /* Extracts the relevant fixing date for the deal */

        if (derivType[i] == CAPFLOOR)
        {
            /* Gets the last period start: not last date, the one before */
            /* Calibration on the last caplet of the cap */
            last_period_start = date_list.date[date_list.len - 2];
            lLastFixing = add_unit(last_period_start, -Spot_Lag, SRT_BDAY, MODIFIED_SUCCEEDING);
        }
        else if ((derivType[i] == SWAPTION) || (derivType[i] == BOND_OPTION))
        {
            /* The Sigma date is the fixng date of the swaption */
            lLastFixing = add_unit(sdpArray[i].start, -Spot_Lag, SRT_BDAY, MODIFIED_SUCCEEDING);
        }
        /* We only keep track of the instruments with a fixing date before or equal to today */
        if (lLastFixing > clcndate)
        {
            /* Compute the forward rate and the spread for the forward period underlying the
             * instrument */
            err = swp_f_ForwardRate_SwapDP(
                &(sdpArray[i]), szYieldCurveName, refRateStr[i], &dForwardRate);
            /* Compute The ATM price */
            if (derivType[i] == CAPFLOOR)
            {
                err = swp_f_CapFloor_SwapDP(
                    &(sdpArray[i]),
                    dForwardRate,
                    pfGetVol,
                    rec_pay[i],
                    refRateStr[i],
                    szYieldCurveName,
                    PREMIUM,
                    eLogOrNorm,
                    &ATMprice[i]);
                if (err)
                {
                    srt_free(date_list.date);
                    FREE_WRAP_BOOTSTRAP_MEMORY;
                    return err;
                }
            }
            else if ((derivType[i] == SWAPTION) || (derivType[i] == BOND_OPTION))
            {
                err     = swp_f_Level_SwapDP(&(sdpArray[i]), szYieldCurveName, &dLevel);
                dfStart = swp_f_df(clcndate, date_list.date[0], szYieldCurveName);
                dfEnd   = swp_f_df(clcndate, date_list.date[date_list.len - 1], szYieldCurveName);
                dForwardCash = (dfStart - dfEnd) / dLevel;
                dSpread      = dForwardRate - dForwardCash;
                dStart       = DTOL(date_list.date[0]);
                dEnd         = DTOL(date_list.date[date_list.len - 1]);
                err          = pfGetVol(
                    (double)dStart, (double)dEnd, dForwardRate, dForwardRate, dSpread, &dATMVol);
                if (err)
                {
                    srt_free(date_list.date);
                    FREE_WRAP_BOOTSTRAP_MEMORY;
                    return err;
                }
                err = swp_f_Swaption_SwapDP(
                    &(sdpArray[i]),
                    dATMVol,
                    dForwardRate,
                    rec_pay[i],
                    refRateStr[i],
                    szYieldCurveName,
                    PREMIUM,
                    eLogOrNorm,
                    &ATMprice[i]);
                if (err)
                {
                    srt_free(date_list.date);
                    FREE_WRAP_BOOTSTRAP_MEMORY;
                    return err;
                }
            }
        }

        srt_free(date_list.date);
    }

    /* Call to the main function for calibration by bootstrap */
    err = srt_f_bootstrap_calibrate(
        und,
        &grfnparam,
        sdpArray,
        strike,
        bondStrike,
        derivType,
        rec_pay,
        refRateStr,
        undPrice,
        ATMprice,
        tau,
        alpha,
        beta,
        rho,
        nInputs);

    /* Need to return the Term Structure as a 2D range. Read it out of und and */
    /* place it into a return XLOPER. */
    FREE_WRAP_BOOTSTRAP_MEMORY;

    /* Return a success message */
    return err;

} /* END Err WrapSrtBootstrap(...) */

/* The new function for the Flexi Cap and Knock in Swap */

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
    char*   undName)
{
    Err              err    = NULL;
    int              status = 0;
    SrtUndPtr        und;
    SrtGrfnParam     grfnparam;
    int              i;
    SwapDP*          sdpArray = NULL;
    int*             derivType;
    SrtReceiverType* rec_pay;

    /* Overwrite defaults with user defined parameters */

    if (err = srt_f_set_GrfnParams(numParams, paramStrings, valueStrings, &grfnparam))
    {
        return err;
    }

    if (undName)
    {
        und = lookup_und(undName);
    }
    else
    {
        return serror("No Underlying passed to SrtCalboot");
    }

    if (und)
    {
        if (!ISUNDTYPE(und, INTEREST_RATE_UND))
        {
            return serror("Underlying must be of type Interest Rate");
        }
    }
    else
    {
        return serror("Could not find underlying in market list");
    }

    /* interpret strings */

    rec_pay   = (int*)calloc(*nInputs, sizeof(int));
    derivType = (int*)calloc(*nInputs, sizeof(int));

    if (!rec_pay || !derivType)
    {
        return serror("Memory allocation failure in Bootstrap");
    }

    for (i = 0; i < *nInputs; i++)
    {
        err = interp_rec_pay(recPayStr[i], &rec_pay[i]);
        if (err)
        {
            return err;
        }
        err = interp_struct(typeStr[i], &derivType[i]);
        if (err)
        {
            return err;
        }
    }

    sdpArray = (SwapDP*)calloc(*nInputs, sizeof(SwapDP));

    for (i = 0; i < *nInputs; i++)
    {
        err = swp_f_initSwapDP(startDates[i], endDates[i], cpdStr[i], basisStr[i], &sdpArray[i]);
        if (err)
            return err;

        sdpArray[i].spot_lag = get_spotlag_from_underlying(und);
    }

    /* Call to the main function for calibration by bootstrap */
    err = srt_f_bootstrap_calibrate(
        und,
        &grfnparam,
        sdpArray,
        strike,
        bondStrike,
        derivType,
        rec_pay,
        refRateStr,
        undPrice,
        ATMprice,
        tau,
        alpha,
        beta,
        rho,
        nInputs);

    /* Need to return the Term Structure as a 2D range. Read it out of und and */
    /* place it into a return XLOPER. */

    if (sdpArray)
    {
        srt_free(sdpArray);
    }

    srt_free(rec_pay);
    srt_free(derivType);

    /* Return a success message */
    return err;

} /* END Err SrtBootstrap(...) */