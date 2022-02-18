/* ============================================================================
   FILENAME:  srt_f_bootstrap.c

   PURPOSE:   generate volatility structure for an IR model using bootstrap
              method - Assume that the mean-reversion factor lambda (or tau)
              is fixed, as well as the correlation for a two factor
   ============================================================================ */

#include "grf_h_public.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_grfclsdfrm.h"
#include "stdio.h"
#include "stdlib.h"

#define DERIV_INTERVAL 0.0001
#define BOOTSTRAP_INT_CHE 0.05
#define ACCURACY 0.002
#define MAXITERNEWTON 20
#define TINY_SIGMA 1.0e-08
#define BOOTSTRAP_MAX_STDEV 4.0

static double **        s_tsv = NULL, **s_tsl = NULL;
static long             s_sig_cols, s_tau_cols, s_tau_rows;
static int*             s_type;
static SrtReceiverType* s_rec_pay;
static double *         s_strike, *s_bond_strike, *s_mkt_price;
static String*          s_ref_rate_code;
static SrtUndPtr        s_und;
static SrtGrfnParam*    s_grfnparams;
static SrtMdlType       s_mdltype;
static SrtMdlDim        s_mdldim;
static double           s_fixedalpha, s_fixedgamma, s_fixedrho;
static double           s_fixedbeta;
static SwapDP*          s_sdparray;
static Date             s_clcndate;

static int s_lCurrentSigmaIndex;
static int s_global_num_inp;

/* --------------------------------------------------------------------------- */

#define FREE_GRFN_BOOTSTRAP_MEMORY                                     \
    {                                                                  \
        if (pdFixingDatesFree)                                         \
            free_dvector(pdFixingDatesFree, 0, (*num_inp) - 1);        \
        pdFixingDatesFree = NULL;                                      \
        if (pdGuessVolatilityFree)                                     \
            free_dvector(pdGuessVolatilityFree, 0, (*num_inp) - 1);    \
        pdGuessVolatilityFree = NULL;                                  \
        if (s_mdldim == ONE_FAC)                                       \
        {                                                              \
            if (s_tsv)                                                 \
                free_dmatrix(s_tsv, 0, 1, 0, iNumUsedInstruments - 1); \
            s_tsv = NULL;                                              \
            if (s_tsl)                                                 \
                free_dmatrix(s_tsl, 0, 1, 0, 0);                       \
            s_tsl = NULL;                                              \
        }                                                              \
        else if (s_mdldim == TWO_FAC)                                  \
        {                                                              \
            if (s_tsv)                                                 \
                free_dmatrix(s_tsv, 0, 3, 0, iNumUsedInstruments - 1); \
            s_tsv = NULL;                                              \
            if (s_tsl)                                                 \
                free_dmatrix(s_tsl, 0, 2, 0, 0);                       \
            s_tsl = NULL;                                              \
        }                                                              \
    }

/* --------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
   Function to be used with the Newton-Raphson routines:
   returns the value of the function f in point x
   In this case, f is the price of the asset, and x is the bucket vol s_tsv[1][s_lCurrentSigmaIndex]
   --------------------------------------------------------------------------- */
static Err func(double x, double* f_x)
{
    Err         err = NULL;
    double      answer;
    TermStruct* ts;
    double      old_tsv;

    /* Keep track of the old sigma value */
    old_tsv = s_tsv[1][s_lCurrentSigmaIndex];

    /* Set the new sigam value */
    s_tsv[1][s_lCurrentSigmaIndex] = x;

    /* If two factor model, work at fixed alpha */
    if (s_mdldim == TWO_FAC)
    {
        s_tsv[2][s_lCurrentSigmaIndex] = x * s_fixedalpha;
    }

    /* Initialise the Term Struct */
    err = srt_f_init_IRM_TermStruct(
        s_clcndate,
        s_tsv,
        s_sig_cols,
        s_lCurrentSigmaIndex + 1,
        s_tsl,
        s_tau_cols,
        s_tau_rows,
        s_mdltype,
        s_mdldim,
        s_fixedbeta,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0, /* vasicek parms */
        0,
        0,
        NULL,
        &ts);
    if (err)
        return err;

    /* Attach the TS to the underlying */
    set_irund_ts(s_und, ts);

    /* Price the deal using Grfn closedform */
    if (err = srt_f_grfn_clsdfrm(
            s_und,
            s_grfnparams,
            &s_sdparray[s_lCurrentSigmaIndex],
            s_strike[s_lCurrentSigmaIndex],
            s_bond_strike[s_lCurrentSigmaIndex],
            s_rec_pay[s_lCurrentSigmaIndex],
            s_type[s_lCurrentSigmaIndex],
            s_ref_rate_code[s_lCurrentSigmaIndex],
            &answer))
    {
        return err;
    }

    /* Free the TS */
    free_underlying_ts(s_und);

    /* Reset the TS to its original value */
    s_tsv[1][s_lCurrentSigmaIndex] = old_tsv;

    *f_x = answer;

    return NULL;
}

/* ----------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------- */

/* Temp New Function For FlexCap and Knock In Swap          */

Err srt_f_bootstrap_calibrate(
    SrtUndPtr        undptr,
    SrtGrfnParam*    grfnparam,
    SwapDP*          sdp,
    double*          pdStrikeRate,
    double*          bnd_strike,
    int*             type,
    SrtReceiverType* recpay,
    String*          pszRefRateCode,
    double*          price,
    double*          ATMprice,
    double           tau,
    double           alpha,
    double           gamma,
    double           rho,
    int*             num_inp)
{
    Err              err;
    int              i, j;
    double           bracket = BOOTSTRAP_INT_CHE, acc = ACCURACY;
    TermStruct*      ts;
    TermStruct*      initial_ts;
    SrtCurvePtr      yldcrv;
    SrtMdlType       mdl_type;
    SrtMdlDim        mdl_dim;
    SrtDateList      date_list;
    String           und_name;
    String           szYieldCurveName;
    SrtCcyParam*     ccy_param;
    int              spot_lag;
    long             last_period_start, lLastFixing;
    long             dEnd, dStart;
    double           short_rate;
    double           dBsImpliedVolatility;
    double           dForwardRate;
    double           dMaturity;
    double           tau_correction;
    double           a[3], b[3];
    long             lNumNewtonIterations;
    double           nstop;
    double           covar_adjustment;
    double           zero_vol_price;
    SRT_Boolean      bPreviousPriceNotFitted;
    Date             clcndate;
    int              iNumUsedInstruments;
    int              iNumTrueInstruments;
    double*          pdFixingDates         = NULL;
    double*          pdFixingDatesFree     = NULL;
    double*          pdGuessVolatility     = NULL;
    double*          pdGuessVolatilityFree = NULL;
    SrtDiffusionType eLogOrNorm;
    double           dModelVol;
    double           dStdDev;
    double           dDistance;
    double           beta1, beta2;

    /* Gets CCY and UND parameters */
    und_name         = get_underlying_name(undptr);
    szYieldCurveName = get_ycname_from_irund(undptr);
    yldcrv           = lookup_curve(szYieldCurveName);
    ccy_param        = get_ccyparam_from_yldcrv(yldcrv);
    spot_lag         = ccy_param->spot_lag;
    undptr           = lookup_und(und_name);
    clcndate         = get_today_from_underlying(undptr);
    get_underlying_mdltype(undptr, &mdl_type);
    get_underlying_mdldim(undptr, &mdl_dim);
    yldcrv = lookup_curve(szYieldCurveName);

    /* Description of Tau and Vol structure */
    s_lCurrentSigmaIndex = 0;
    s_tau_rows           = 1;

    /* Aloocate space for Tau (the tau date is taken from the last input instrument) */
    if (mdl_dim == ONE_FAC)
    {
        s_tsl      = dmatrix(0, 1, 0, 0);
        s_tau_cols = 2;
    }
    else if (mdl_dim == TWO_FAC)
    {
        s_tsl      = dmatrix(0, 2, 0, 0);
        s_tau_cols = 3;
    }

    /* Set the Beta (assuming it is fixed) for the underlying */
    err = get_underlying_ts(undptr, &initial_ts);
    if (err)
    {
        return err;
    }

    if (mdl_dim == ONE_FAC)
    {
        s_fixedbeta = find_beta((Date)clcndate, initial_ts);
    }
    else if (mdl_dim == TWO_FAC)
    {
        err = find_tf_beta((Date)clcndate, initial_ts, &beta1, &beta2);
        if (err)
        {
            return err;
        }
        s_fixedbeta = (beta1 + beta2) / 2.;
    }

    /* Set the Number of Instruments equal to (*num_inp) by default */
    iNumTrueInstruments = (*num_inp);

    /* Sets the Tau value */
    s_tsl[0][0] = (double)(sdp[(*num_inp) - 1].end);
    s_tsl[1][0] = tau;
    if (mdl_dim == TWO_FAC)
    {
        s_tsl[2][0] = 1.0 / (1.0 / tau + gamma);
    }

    /* Allocate memory for the fixing dates */
    pdFixingDates         = dvector(0, (*num_inp) - 1);
    pdGuessVolatility     = dvector(0, (*num_inp) - 1);
    pdFixingDatesFree     = pdFixingDates;
    pdGuessVolatilityFree = pdGuessVolatility;

    /* Sets the covariance adjustment for the volatility guess */
    covar_adjustment = 1.0;
    if (mdl_dim == TWO_FAC)
    {
        covar_adjustment = sqrt(1.0 + alpha * (alpha + 2.0 * rho));
    }

    /* Sets the Volatility Type according to the model */
    if (mdl_type == LGM)
    {
        eLogOrNorm = SRT_NORMAL;
    }
    else
    {
        eLogOrNorm = SRT_LOGNORMAL;
    }

    /* Sets an initial guess only for instruments with a fixing date after today */
    iNumUsedInstruments = 0;
    for (i = 0; i < (*num_inp); i++)
    {
        /* Builds the schedule of the Cap/Floor */
        date_list = SwapDP_to_DateList(&(sdp[iNumUsedInstruments]), MODIFIED_SUCCEEDING);
        dStart    = date_list.date[0];
        dEnd      = date_list.date[date_list.len - 1];

        /* Extracts the relevant fixing date for the deal */
        if (type[iNumUsedInstruments] == CAPFLOOR)
        {
            /* Gets the last period start: not last date, the one before */
            /* Calibration on the last caplet of the cap */
            last_period_start = date_list.date[date_list.len - 2];
            lLastFixing = add_unit(last_period_start, -spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
        }
        else if (
            (type[iNumUsedInstruments] == SWAPTION) || (type[iNumUsedInstruments] == BOND_OPTION))
        {
            /* The Sigma date is the fixng date of the swaption */
            lLastFixing = add_unit(date_list.date[0], -spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
        }
        srt_free(date_list.date);
        /* Stores the Fixing Date */
        pdFixingDates[iNumUsedInstruments] = lLastFixing;

        /* We only keep track of the instruments with a fixing date before or equal to today */
        if (lLastFixing > clcndate)
        {
            /* Compute the forward rate for the forward period underlying the instrument */
            err = swp_f_ForwardRate_SwapDP(
                &(sdp[iNumUsedInstruments]),
                szYieldCurveName,
                pszRefRateCode[iNumUsedInstruments],
                &dForwardRate);
            if (err)
            {
                FREE_GRFN_BOOTSTRAP_MEMORY;
                return err;
            }

            /* Get the Black implied volatility from the market price */
            dBsImpliedVolatility = 1;
            if (type[iNumUsedInstruments] == SWAPTION)
            {
                err = swp_f_SwaptionImpliedVol_SwapDP(
                    ATMprice[iNumUsedInstruments],
                    &(sdp[iNumUsedInstruments]),
                    dForwardRate,
                    recpay[iNumUsedInstruments],
                    pszRefRateCode[iNumUsedInstruments],
                    szYieldCurveName,
                    eLogOrNorm,
                    &dBsImpliedVolatility);
                if (err)
                {
                    FREE_GRFN_BOOTSTRAP_MEMORY;
                    return err;
                }
            }
            else
            {
                err = swp_f_CapFloorImpVol_SwapDP(
                    ATMprice[iNumUsedInstruments],
                    &(sdp[iNumUsedInstruments]),
                    dForwardRate,
                    recpay[iNumUsedInstruments],
                    pszRefRateCode[iNumUsedInstruments],
                    szYieldCurveName,
                    eLogOrNorm,
                    &dBsImpliedVolatility);
                if (err)
                {
                    FREE_GRFN_BOOTSTRAP_MEMORY;
                    return err;
                }
            }

            /* The maturity of the Option */
            dMaturity = (lLastFixing - clcndate) * YEARS_IN_DAY;

            /* Computes the Distance betwenn Forward and Strike in Num Std Dev */
            dStdDev = sqrt(dMaturity) * dBsImpliedVolatility;
            if (eLogOrNorm == SRT_NORMAL)
            {
                dDistance = fabs(dForwardRate - pdStrikeRate[iNumUsedInstruments]) / dStdDev;
            }
            else
            {
                dDistance = fabs(log(dForwardRate / pdStrikeRate[iNumUsedInstruments])) / dStdDev;
            }

            /* If the Strike is more than MAX_STDEV away from the Forward, discard the instrument */
            if (dDistance > BOOTSTRAP_MAX_STDEV)
            {
                /* Removes this instrument from the list */
                for (j = iNumUsedInstruments; j < iNumTrueInstruments - 1; j++)
                {
                    sdp[j]          = sdp[j + 1];
                    pdStrikeRate[j] = pdStrikeRate[j + 1];
                    bnd_strike[j]   = bnd_strike[j + 1];
                    type[j]         = type[j + 1];
                    recpay[j]       = recpay[j + 1];
                    price[j]        = price[j + 1];
                    ATMprice[j]     = ATMprice[j + 1];
                    strcpy(pszRefRateCode[j], pszRefRateCode[j + 1]);
                }
                continue;
            }

            /* If The Strike is Okay then compute the vol     */

            /* Get the Black implied volatility from the market price */
            dBsImpliedVolatility = 1;
            if (type[iNumUsedInstruments] == SWAPTION)
            {
                err = swp_f_SwaptionImpliedVol_SwapDP(
                    price[iNumUsedInstruments],
                    &(sdp[iNumUsedInstruments]),
                    pdStrikeRate[iNumUsedInstruments],
                    recpay[iNumUsedInstruments],
                    pszRefRateCode[iNumUsedInstruments],
                    szYieldCurveName,
                    eLogOrNorm,
                    &dBsImpliedVolatility);
                if (err)
                {
                    FREE_GRFN_BOOTSTRAP_MEMORY;
                    return err;
                }
            }
            else
            {
                err = swp_f_CapFloorImpVol_SwapDP(
                    price[iNumUsedInstruments],
                    &(sdp[iNumUsedInstruments]),
                    pdStrikeRate[iNumUsedInstruments],
                    recpay[iNumUsedInstruments],
                    pszRefRateCode[iNumUsedInstruments],
                    szYieldCurveName,
                    eLogOrNorm,
                    &dBsImpliedVolatility);
                if (err)
                {
                    FREE_GRFN_BOOTSTRAP_MEMORY;
                    return err;
                }
            }
            /* Compute the IFR for the start date of the instrument */
            short_rate = swp_f_fwdcash(
                sdp[iNumUsedInstruments].start,
                sdp[iNumUsedInstruments].start + 14.0,
                sdp[iNumUsedInstruments].basis_code,
                szYieldCurveName);
            /*
                                    err =
               swp_f_Level_SwapDP(&(sdpArray[iNumUsedInstruments]),szYieldCurveName,&dLevel);
                                    dfStart = swp_f_df(clcndate,dStart,szYieldCurveName);
                                    dfEnd = swp_f_df(clcndate,dEnd,szYieldCurveName);
                                    dForwardCash = (dfStart - dfEnd) / dLevel;
            */
            /* Compute a tau correction for the volatility guess ( = exp(time/tau) ) */
            if (1 / tau < 1.0e-10)
                tau_correction = 1.0;
            else
                tau_correction = exp(1.0 / tau * dMaturity);

            /* Set the initial guess, using a rough tau and covariance correction */
            if ((mdl_type == CHEY) || (mdl_type == CHEY_BETA))
            {
                /* LogNormal Or Beta Case */
                dModelVol = dBsImpliedVolatility * dForwardRate / pow(short_rate, s_fixedbeta);
            }
            else if (mdl_type == LGM)
            {
                /* Normal case */
                dModelVol = dBsImpliedVolatility;
            }
            dModelVol /= covar_adjustment;

            pdGuessVolatility[iNumUsedInstruments] = tau_correction * dModelVol;

            /* Copy the value of parameters in the s_...pointers */

            /* Increment the Number of instuments by one */
            iNumUsedInstruments++;

        } /* END if (lLastFixing >clcndate ) */
        else
        {
            /* Removes this instrument from the list */
            sdp++;
            pdStrikeRate++;
            bnd_strike++;
            type++;
            recpay++;
            price++;
            ATMprice++;
            pszRefRateCode++;
            iNumTrueInstruments--;
        }

    } /* END 	for (i = 0; i < (*num_inp); i++) loop */

    if ((iNumTrueInstruments <= 0) || (iNumUsedInstruments <= 0))
    {
        FREE_GRFN_BOOTSTRAP_MEMORY;
        return serror("No instruments used for Calibration !");
    }

    /* Initialization: assign input parameters to global static variables */
    s_und           = undptr;
    s_clcndate      = clcndate;
    s_grfnparams    = grfnparam;
    s_sdparray      = sdp;
    s_strike        = pdStrikeRate;
    s_bond_strike   = bnd_strike;
    s_type          = type;
    s_rec_pay       = recpay;
    s_mkt_price     = price;
    s_ref_rate_code = pszRefRateCode;
    s_mdltype       = mdl_type;
    s_mdldim        = mdl_dim;

    /* Sets the Two Factor Parameters */
    s_fixedalpha = alpha;
    s_fixedgamma = gamma;
    s_fixedrho   = rho;

    /* The number of prices to fit */
    s_global_num_inp = iNumUsedInstruments;

    /* Allocate space for the sig dates and values and sets the number of columns */
    if (mdl_dim == ONE_FAC)
    {
        s_tsv      = dmatrix(0, 1, 0, iNumUsedInstruments - 1);
        s_sig_cols = 2;
    }
    else if (mdl_dim == TWO_FAC)
    {
        s_tsv      = dmatrix(0, 3, 0, iNumUsedInstruments - 1);
        s_sig_cols = 4;

        /* If a two factor model, set the fixed correlation in the TermStruct*/
        for (i = 0; i < iNumUsedInstruments; i++)
            s_tsv[3][i] = rho;
    }

    /* Initialisations before the fit */
    bPreviousPriceNotFitted = SRT_FALSE;
    s_lCurrentSigmaIndex    = -1;

    /* Free the TS */
    free_underlying_ts(undptr);

    /* Starts the loop on all market instruments */
    for (i = 0; i < iNumUsedInstruments; i++)
    {
        /* Increment the Sigma Row (in the loop) for the TermStructure IF REQUIRED */
        if (bPreviousPriceNotFitted == SRT_FALSE)
        {
            /* Move on to the next Sigma */
            s_lCurrentSigmaIndex++;
        }

        /* Set the relevant sigma date for the deal with a NULL value */
        s_tsv[0][s_lCurrentSigmaIndex] = pdFixingDates[i];
        s_tsv[1][s_lCurrentSigmaIndex] = TINY_SIGMA;

        /* Computes the model price with a zero vol  */
        err = func(TINY_SIGMA, &zero_vol_price);
        if (err)
            return err;

        /* Check this price is not higher than the market price to fit: if it is set vol to 0.0 and
         * continue */
        if (zero_vol_price > price[i])
        {
            s_tsv[1][s_lCurrentSigmaIndex] = TINY_SIGMA;
            bPreviousPriceNotFitted        = SRT_FALSE;
            continue;
        }

        /* Sets the initial guess as a value for sigma */
        s_tsv[1][s_lCurrentSigmaIndex] = pdGuessVolatility[i];

        /* A quick and usual Newton routine */
        a[0] = s_tsv[1][s_lCurrentSigmaIndex];
        err  = func(a[0], &b[0]);
        if (err)
        {
            FREE_GRFN_BOOTSTRAP_MEMORY;
            return err;
        }
        a[1] = a[0] * (1 + 0.05);
        err  = func(a[1], &b[1]);
        if (err)
        {
            FREE_GRFN_BOOTSTRAP_MEMORY;
            return err;
        }
        lNumNewtonIterations = 0;
        nstop                = 0.0;
        a[2]                 = a[1] * (1 + 0.05);
        while ((lNumNewtonIterations < MAXITERNEWTON) && (nstop < 1.0))
        {
            err = func(a[2], &b[2]);
            if (err)
            {
                FREE_GRFN_BOOTSTRAP_MEMORY;
                return err;
            }
            newton(price[i], 3, a, b, &nstop);
            if (a[2] < 0.0)
                a[2] = a[1] / 10;

            lNumNewtonIterations++;
        }

        /* If not over the number of iterations: set the vol , else skip this price */

        if (lNumNewtonIterations < MAXITERNEWTON)
        {
            s_tsv[1][s_lCurrentSigmaIndex] = a[2];

            /* Reset the second volatility if it is a Two Factor model */
            if (mdl_dim == TWO_FAC)
            {
                s_tsv[2][s_lCurrentSigmaIndex] = s_tsv[1][s_lCurrentSigmaIndex] * alpha;
            }
        }
        else
        {
            /* If the Newton Algorithm didn't converge put the Vol at TINY_SIGMA */
            s_tsv[1][s_lCurrentSigmaIndex] = TINY_SIGMA;

            /* Reset the second volatility if it is a Two Factor model */
            if (mdl_dim == TWO_FAC)
            {
                s_tsv[2][s_lCurrentSigmaIndex] = TINY_SIGMA * alpha;
            }
            /*	bPreviousPriceNotFitted = SRT_TRUE; */
        }

    } /* END loop on all prices */

    /* Rebuilds the full term structure for the underlying */
    if (err = srt_f_init_IRM_TermStruct(
            s_clcndate,
            s_tsv,
            s_sig_cols,
            s_lCurrentSigmaIndex + 1,
            s_tsl,
            s_tau_cols,
            s_tau_rows,
            mdl_type,
            s_mdldim,
            s_fixedbeta,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0, /* vasicek parms */
            0,
            0,
            NULL,
            &ts))

    {
        FREE_GRFN_BOOTSTRAP_MEMORY;
        return err;
    }

    /* Attaches the new term structure back to the underlying */
    set_irund_ts(s_und, ts);

    /* Increments the underlying ticker (it is a new one after all) */
    get_underlying_ticker(undptr) += 1;

    /* Free the memory */
    FREE_GRFN_BOOTSTRAP_MEMORY;

    /* Return a success message */
    return NULL;

} /* END Err srt_f_bootstrap_calibrate(...) */