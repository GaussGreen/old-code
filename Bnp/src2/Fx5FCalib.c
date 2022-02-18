/* ==========================================================================
   FILE_NAME:	Fx5FCalib.c

   PURPOSE:		Modelling of the spot FX vol by taking in consideration a LGM 2 factor
                                on the domestic and foreign market and a lognormal model on the Spot
   Fx

   DATE:		11/01/01

   AUTHOR:		L.C.
   ========================================================================== */

/*	All this functions are using a correlation matrix double **correl
        the meaning of this correlation is:
        index 0: domestic IR, factor 1
        index 1: domestic IR, factor 2
        index 2: foreign IR, factor 1
        index 3: foreign IR, factor 2
        index 4: spot Fx
 */

#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

/*	A set of static functions designed to calculate exponential integrals */

static double Phi_Func(double x, double T, double s, double t)
{
    double result;

    result = (exp(-x * (T - t)) - exp(-x * (T - s))) / x;

    return result;
}

static double Etha_Func(double x, double T, double s, double t)
{
    double result;

    result = (t - s - Phi_Func(x, T, s, t)) / x;

    return result;
}

static double B_Func(double x, double T, double s, double t)
{
    double result;

    result = -(t * t - s * s) / 2.0 +
             1.0 / x * ((t - 1.0 / x) * exp(-x * (T - t)) - (s - 1.0 / x) * exp(-x * (T - s)));

    return result;
}

static double Psi_Func(double x, double y, double T, double s, double t)
{
    double result;

    result = 1.0 / (x * y) *
             (t - s - Phi_Func(x, T, s, t) - Phi_Func(y, T, s, t) + Phi_Func(x + y, T, s, t));

    return result;
}

static double Psi2_Func(double x, double y, double Tx, double Ty, double s, double t)
{
    double result;

    result = 1.0 / (x * y) *
             (t - s - Phi_Func(x, Tx, s, t) - Phi_Func(y, Ty, s, t) +
              exp(-x * (Tx - Ty)) * Phi_Func(x + y, Ty, s, t));

    return result;
}

static double H_func(double lam, double t1, double t2)
{
    return exp(lam * t2) - exp(lam * t1);
}

/*	A set of useful functions to calculate variances */

static Err Coefs_Partial_Var(
    double   T,
    double   T1,
    double   T2,
    double   sig_dom1,
    double   sig_dom2,
    double   lda_dom1,
    double   lda_dom2,
    double   sig_for1,
    double   sig_for2,
    double   lda_for1,
    double   lda_for2,
    double** correl,
    double*  a,
    double*  b,
    double*  c)
{
    /* coefficient of sig_fx 2 */
    (*a) = T2 - T1;

    /* coefficient of sig_fx   */
    (*b) = -2.0 * (sig_for1 * correl[2][4] * Etha_Func(lda_for1, T, T1, T2) +
                   sig_for2 * correl[3][4] * Etha_Func(lda_for2, T, T1, T2)) +
           2.0 * (sig_dom1 * correl[0][4] * Etha_Func(lda_dom1, T, T1, T2) +
                  sig_dom2 * correl[1][4] * Etha_Func(lda_dom2, T, T1, T2));

    /* constant coefficient    */
    (*c) = sig_for1 * sig_for1 * Psi_Func(lda_for1, lda_for1, T, T1, T2) +
           sig_for2 * sig_for2 * Psi_Func(lda_for2, lda_for2, T, T1, T2) +
           sig_dom1 * sig_dom1 * Psi_Func(lda_dom1, lda_dom1, T, T1, T2) +
           sig_dom2 * sig_dom2 * Psi_Func(lda_dom2, lda_dom2, T, T1, T2) -
           2.0 * correl[0][2] * sig_dom1 * sig_for1 * Psi_Func(lda_for1, lda_dom1, T, T1, T2) -
           2.0 * correl[0][3] * sig_dom1 * sig_for2 * Psi_Func(lda_for2, lda_dom1, T, T1, T2) +
           2.0 * correl[0][1] * sig_dom1 * sig_dom2 * Psi_Func(lda_dom1, lda_dom2, T, T1, T2) -
           2.0 * correl[1][2] * sig_dom2 * sig_for1 * Psi_Func(lda_for1, lda_dom2, T, T1, T2) -
           2.0 * correl[1][3] * sig_dom2 * sig_for2 * Psi_Func(lda_for2, lda_dom2, T, T1, T2) +
           2.0 * correl[2][3] * sig_for2 * sig_for1 * Psi_Func(lda_for1, lda_for2, T, T1, T2);

    return NULL;
}

static Err Partial_Var(
    double   T,
    double   T1,
    double   T2,
    double   sig_dom1,
    double   sig_dom2,
    double   lda_dom1,
    double   lda_dom2,
    double   sig_for1,
    double   sig_for2,
    double   lda_for1,
    double   lda_for2,
    double   sig_fx,
    double** correl,
    double*  var)
{
    double a;
    double b;
    double c;
    Err    err = NULL;

    /* Get Coefs */
    err = Coefs_Partial_Var(
        T,
        T1,
        T2,
        sig_dom1,
        sig_dom2,
        lda_dom1,
        lda_dom2,
        sig_for1,
        sig_for2,
        lda_for1,
        lda_for2,
        correl,
        &a,
        &b,
        &c);

    if (err)
    {
        return err;
    }

    (*var) = sig_fx * (a * sig_fx + b) + c;

    return NULL;
}

/*	Fx Implied volatility */
/*  The rates maturity dates have to be merged ! */
Err Fx5DtsImpliedVol(
    double   opt_maturity,
    double   start_date,
    double   end_date,
    double*  maturity_rates,
    long     nbMat,
    double*  sig_curve_dom,
    double   lda_dom,
    double   alpha_dom,
    double   gamma_dom,
    double*  sig_curve_for,
    double   lda_for,
    double   alpha_for,
    double   gamma_for,
    double*  maturity_fx,
    double*  sig_curve_fx,
    long     nbrMat_fx,
    double** correl,
    double*  fx_vol)
{
    double sig_dom1, sig_dom2, lda_dom2, sig_for1, sig_for2, lda_for2, sig_fx;

    double T1, T2, t1, t2;
    double var;
    double var_partial;
    int    i, j;
    long   StartIndex, EndIndex, StartIndex2, EndIndex2;
    Err    err = NULL;

    lda_dom2 = lda_dom + gamma_dom;
    lda_for2 = lda_for + gamma_for;

    if (start_date > end_date)
    {
        err = "end_date before start_date in Fx5DtsImpliedVol";
        return err;
    }

    if (end_date == 0)
    {
        (*fx_vol) = 0;
        return err;
    }

    StartIndex = Get_Index(start_date, maturity_fx, nbrMat_fx);
    EndIndex   = Get_Index(end_date, maturity_fx, nbrMat_fx);

    var = 0;

    for (i = StartIndex; i < EndIndex + 1; i++)
    {
        if (i > StartIndex)
        {
            T1 = maturity_fx[i - 1];
        }
        else
        {
            /* First part */
            T1 = start_date;
        }

        if (i == EndIndex || StartIndex == EndIndex)
        {
            /* Last part */
            T2 = end_date;
        }
        else
        {
            T2 = maturity_fx[i];
        }

        StartIndex2 = Get_Index(T1, maturity_rates, nbMat);
        EndIndex2   = Get_Index(T2, maturity_rates, nbMat);

        sig_fx = sig_curve_fx[i];

        for (j = StartIndex2; j < EndIndex2 + 1; j++)
        {
            if (j > StartIndex2)
            {
                t1 = maturity_rates[j - 1];
            }
            else
            {
                /* First part */
                t1 = T1;
            }

            if (j == EndIndex2 || StartIndex2 == EndIndex2)
            {
                /* Last part */
                t2 = T2;
            }
            else
            {
                t2 = maturity_rates[j];
            }

            sig_dom1 = sig_curve_dom[j];
            sig_dom2 = sig_dom1 * alpha_dom;
            sig_for1 = sig_curve_for[j];
            sig_for2 = sig_for1 * alpha_for;

            err = Partial_Var(
                opt_maturity,
                t1,
                t2,
                sig_dom1,
                sig_dom2,
                lda_dom,
                lda_dom2,
                sig_for1,
                sig_for2,
                lda_for,
                lda_for2,
                sig_fx,
                correl,
                &var_partial);

            if (err)
            {
                return err;
            }

            var += var_partial;
        }
    }

    if (fabs(end_date - start_date) > 1.0e-08)
    {
        *fx_vol = sqrt(var / (end_date - start_date));
    }
    else
    {
        *fx_vol = 0.0;
    }

    return err;
}

/*	Calibration of a fx term structure to a set of fx options  */
Err Fx5DtsCalibration(
    double*  exercise_opt,
    double*  maturity_opt,
    double*  vol_opt,
    long     nbrOpt,
    double*  maturity_rates,
    long     nbrMat,
    double*  sig_curve_dom,
    double   lda_dom,
    double   alpha_dom,
    double   gamma_dom,
    double*  sig_curve_for,
    double   lda_for,
    double   alpha_for,
    double   gamma_for,
    double** correl,
    double** fx_vol_curve)

{
    double a, b, c, a_part, b_part, c_part, c2, delta;
    double sig_dom1, sig_dom2, lda_dom2, sig_for1, sig_for2, lda_for2;
    double cumvar, vol;
    double T1, T2, t1, t2;
    double VolMin;
    double val_opt;
    int    idxopt, i;
    long   StartIndex, EndIndex;

    Err err = NULL;

    /* loop on the number of options */
    (*fx_vol_curve) = NULL;
    (*fx_vol_curve) = (double*)calloc(nbrOpt, sizeof(double));
    if (!(*fx_vol_curve))
    {
        err = "Memory allocation error in Fx5DtsCalibration";
        goto FREE_RETURN;
    }

    lda_dom2 = lda_dom + gamma_dom;
    lda_for2 = lda_for + gamma_for;

    for (idxopt = 0; idxopt <= nbrOpt - 1; idxopt++)
    {
        T2      = exercise_opt[idxopt];
        val_opt = maturity_opt[idxopt];

        /* Get the cumulated variance till the previous maturity T1*/
        if (idxopt > 0)
        {
            T1 = exercise_opt[idxopt - 1];

            err = Fx5DtsImpliedVol(
                val_opt,
                0,
                T1,
                maturity_rates,
                nbrMat,
                sig_curve_dom,
                lda_dom,
                alpha_dom,
                gamma_dom,
                sig_curve_for,
                lda_for,
                alpha_for,
                gamma_for,
                exercise_opt,
                *fx_vol_curve,
                nbrOpt,
                correl,
                &vol);

            if (err)
            {
                goto FREE_RETURN;
            }

            cumvar = vol * vol * T1;
        }
        else
        {
            T1     = 0;
            cumvar = 0;
        }

        /* Calculate the coefficient of the last cumulative of the variance between T1 and T2 */

        a = b = c2 = 0;
        StartIndex = Get_Index(T1, maturity_rates, nbrMat);
        EndIndex   = Get_Index(T2, maturity_rates, nbrMat);

        for (i = StartIndex; i < EndIndex + 1; i++)
        {
            if (i > StartIndex)
            {
                t1 = maturity_rates[i - 1];
            }
            else
            {
                /* First part */
                t1 = T1;
            }

            if (i == EndIndex || StartIndex == EndIndex)
            {
                /* Last part */
                t2 = T2;
            }
            else
            {
                t2 = maturity_rates[i];
            }

            sig_dom1 = sig_curve_dom[i];
            sig_dom2 = sig_dom1 * alpha_dom;
            sig_for1 = sig_curve_for[i];
            sig_for2 = sig_for1 * alpha_for;

            err = Coefs_Partial_Var(
                val_opt,
                t1,
                t2,
                sig_dom1,
                sig_dom2,
                lda_dom,
                lda_dom2,
                sig_for1,
                sig_for2,
                lda_for,
                lda_for2,
                correl,
                &a_part,
                &b_part,
                &c_part);

            if (err)
            {
                goto FREE_RETURN;
            }

            a += a_part;
            b += b_part;
            c2 += c_part;
        }

        c = c2 + cumvar;
        /* substract value to match */
        c -= vol_opt[idxopt] * vol_opt[idxopt] * T2;

        /* just solve the second order equation */
        delta = b * b - 4 * a * c;
        /* delta < 0 or solutions are negatives */
        if ((delta < 0) || ((c > 0) && (b > 0)))
        {
            VolMin = sqrt((c2 + cumvar - b * b / (4 * a)) / val_opt) * 100;

            err = serror(
                "Cannot find solution to match the option %d (its vol should be > %.2f %%)",
                idxopt + 1,
                VolMin);

            goto FREE_RETURN;
        }

        /* if there is two positive solutions we are taking the smallest one */
        /* it would be easier to match the next one */
        if ((c > 0) && (b < 0))
        {
            (*fx_vol_curve)[idxopt] = (-b - sqrt(delta)) / (2. * a);
        }
        else
        {
            (*fx_vol_curve)[idxopt] = (-b + sqrt(delta)) / (2. * a);
        }
    }

FREE_RETURN:

    if (err)
    {
        if (*fx_vol_curve)
        {
            free(*fx_vol_curve);
            *fx_vol_curve = NULL;
        }
    }

    return err;
}

/*	Implied vol direct from underlying */
Err Fx5DImpliedVol(
    char*    fx_underlying,
    double** correl,
    double   val_time,
    double   start_time,
    double   end_time,
    double*  vol)
{
    long    sigma_n_dom, sigma_n_for, sigma_n_fx, nb_merge_dates;
    double *sigma_date_dom = NULL, *sigma_dom = NULL;
    double *sigma_date_for = NULL, *sigma_for = NULL;
    double *sigma_date_fx = NULL, *sigma_fx = NULL;

    double tau_dom, alpha_dom, gamma_dom, rho_dom;
    double tau_for, alpha_for, gamma_for, rho_for;

    double correl_dom_for, correl_dom_fx, correl_for_fx;
    double lda_dom, lda_for;

    double *sig_dom = NULL, *sig_for = NULL, *merge_dates = NULL;

    Err err = NULL;

    err = Get_FX_StochRate_TermStructures5F(
        fx_underlying,
        &sigma_date_dom,
        &sigma_dom,
        &sigma_n_dom,
        &tau_dom,
        &alpha_dom,
        &gamma_dom,
        &rho_dom,
        &sigma_date_for,
        &sigma_for,
        &sigma_n_for,
        &tau_for,
        &alpha_for,
        &gamma_for,
        &rho_for,
        &sigma_date_fx,
        &sigma_fx,
        &sigma_n_fx,
        &correl_dom_for,
        &correl_dom_fx,
        &correl_for_fx);
    if (err)
    {
        goto FREE_RETURN;
    }

    rho_dom = correl[0][1];
    rho_for = correl[2][3];

    lda_dom = 1.0 / tau_dom;
    lda_for = 1.0 / tau_for;

    err = merge_rates_ts(
        sigma_date_dom,
        sigma_dom,
        sigma_n_dom,
        sigma_date_for,
        sigma_for,
        sigma_n_for,
        &merge_dates,
        &sig_dom,
        &sig_for,
        &nb_merge_dates);

    if (err)
    {
        goto FREE_RETURN;
    }

    err = Fx5DtsImpliedVol(
        val_time,
        start_time,
        end_time,
        merge_dates,
        nb_merge_dates,
        sig_dom,
        lda_dom,
        alpha_dom,
        gamma_dom,
        sig_for,
        lda_for,
        alpha_for,
        gamma_for,
        sigma_date_fx,
        sigma_fx,
        sigma_n_fx,
        correl,
        vol);

    if (err)
    {
        goto FREE_RETURN;
    }

FREE_RETURN:

    if (sigma_date_dom)
    {
        free(sigma_date_dom);
    }

    if (sigma_dom)
    {
        free(sigma_dom);
    }

    if (sigma_date_for)
    {
        free(sigma_date_for);
    }

    if (sigma_for)
    {
        free(sigma_for);
    }

    if (sigma_date_fx)
    {
        free(sigma_date_fx);
    }

    if (sigma_fx)
    {
        free(sigma_fx);
    }

    if (sig_dom)
    {
        free(sig_dom);
    }

    if (sig_for)
    {
        free(sig_for);
    }

    if (merge_dates)
    {
        free(merge_dates);
    }

    return err;
}

Err Fill_5F_correl_matrix(
    double* maturity_rates,
    long    nbrMat,
    double* sig_curve_dom,
    double  lda_dom,
    double  alpha_dom,
    double  gamma_dom,
    double  rho_dom,
    double* sig_curve_for,
    double  lda_for,
    double  alpha_for,
    double  gamma_for,
    double  rho_for,
    double  correl_dom_for,
    double  correl_dom_fx,
    double  correl_for_fx,
    double* correl[])
{
    correl[0][1] = rho_dom;
    correl[0][2] = correl_dom_for;
    correl[0][3] = correl_dom_for;
    correl[0][4] = correl_dom_fx;
    correl[1][2] = correl_dom_for;
    correl[1][3] = correl_dom_for;
    correl[1][4] = correl_dom_fx;
    correl[2][3] = rho_for;
    correl[2][4] = correl_for_fx;

    return NULL;
}

/*	Fx calibration direct from underlying */
Err Fx5DCalibration(
    char*    dom_underlying,
    char*    for_underlying,
    double** correl,
    double*  exercise_opt,
    double*  maturity_opt,
    double*  vol_opt,
    long     nbropt,
    double** fx_vol_curve)
{
    long    sigma_n_dom, sigma_n_for;
    long    nb_merge_dates;
    double *sigma_date_dom = NULL, *sigma_dom = NULL, *sigma_date_for = NULL, *sigma_for = NULL,
           *merge_dates = NULL, *sig_dom = NULL, *sig_for = NULL;

    double tau_dom, lda_dom, alpha_dom, gamma_dom, rho_dom, tau_for, lda_for, alpha_for, gamma_for,
        rho_for;

    Err err = NULL;

    err = Get_LGM2F_TermStructure(
        dom_underlying,
        &sigma_date_dom,
        &sigma_dom,
        &sigma_n_dom,
        &tau_dom,
        &alpha_dom,
        &gamma_dom,
        &rho_dom);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = Get_LGM2F_TermStructure(
        for_underlying,
        &sigma_date_for,
        &sigma_for,
        &sigma_n_for,
        &tau_for,
        &alpha_for,
        &gamma_for,
        &rho_for);
    if (err)
    {
        goto FREE_RETURN;
    }

    rho_dom = correl[0][1];
    rho_for = correl[2][3];

    lda_dom = 1.0 / tau_dom;
    lda_for = 1.0 / tau_for;

    err = merge_rates_ts(
        sigma_date_dom,
        sigma_dom,
        sigma_n_dom,
        sigma_date_for,
        sigma_for,
        sigma_n_for,
        &merge_dates,
        &sig_dom,
        &sig_for,
        &nb_merge_dates);

    if (err)
    {
        goto FREE_RETURN;
    }

    err = Fx5DtsCalibration(
        exercise_opt,
        maturity_opt,
        vol_opt,
        nbropt,
        merge_dates,
        nb_merge_dates,
        sig_dom,
        lda_dom,
        alpha_dom,
        gamma_dom,
        sig_for,
        lda_for,
        alpha_for,
        gamma_for,
        correl,
        fx_vol_curve);

    if (err)
    {
        goto FREE_RETURN;
    }

FREE_RETURN:

    if (sigma_date_dom)
    {
        free(sigma_date_dom);
    }

    if (sigma_dom)
    {
        free(sigma_dom);
    }

    if (sigma_date_for)
    {
        free(sigma_date_for);
    }

    if (sigma_for)
    {
        free(sigma_for);
    }

    if (sig_dom)
    {
        free(sig_dom);
    }

    if (sig_for)
    {
        free(sig_for);
    }

    if (merge_dates)
    {
        free(merge_dates);
    }

    return err;
}

Err Get_FX_StochRate_TermStructures5F(
    char*    underlying,
    double** sigma_date_dom,
    double** sigma_dom,
    long*    sigma_n_dom,
    double*  fixed_tau_dom,
    double*  fixed_alpha_dom,
    double*  fixed_gamma_dom,
    double*  fixed_rho_dom,
    double** sigma_date_for,
    double** sigma_for,
    long*    sigma_n_for,
    double*  fixed_tau_for,
    double*  fixed_alpha_for,
    double*  fixed_gamma_for,
    double*  fixed_rho_for,
    double** sigma_date_fx,
    double** sigma_fx,
    long*    sigma_n_fx,
    double*  correl_dom_for,
    double*  correl_dom_fx,
    double*  correl_for_fx)
{
    SrtUndPtr und;
    char *    domname, *forname;
    long      corr_date_n;
    double *  corr_date = NULL, *corr = NULL;
    long      today;
    int       i;

    Err err = NULL;

    *sigma_date_dom = NULL;
    *sigma_dom      = NULL;
    *sigma_date_for = NULL;
    *sigma_for      = NULL;
    *sigma_date_fx  = NULL;
    *sigma_fx       = NULL;

    und = lookup_und(underlying);
    if (!und)
    {
        err = serror("Couldn't fin underlying named %s", underlying);
        goto FREE_RETURN;
    }

    if (get_underlying_type(und) != FOREX_UND)
    {
        err = serror("Underlying %s is not of type FX", underlying);
        goto FREE_RETURN;
    }

    if (get_mdltype_from_fxund(und) != FX_STOCH_RATES)
    {
        err = serror("Underlying %s is not of type FX Stoch Rates", underlying);
        goto FREE_RETURN;
    }

    domname = get_domname_from_fxund(und);

    err = Get_LGM2F_TermStructure(
        domname,
        sigma_date_dom,
        sigma_dom,
        sigma_n_dom,
        fixed_tau_dom,
        fixed_alpha_dom,
        fixed_gamma_dom,
        fixed_rho_dom);
    if (err)
    {
        goto FREE_RETURN;
    }

    forname = get_forname_from_fxund(und);

    err = Get_LGM2F_TermStructure(
        forname,
        sigma_date_for,
        sigma_for,
        sigma_n_for,
        fixed_tau_for,
        fixed_alpha_for,
        fixed_gamma_for,
        fixed_rho_for);

    if (err)
    {
        goto FREE_RETURN;
    }

    err = srt_f_display_FX_TermStruct(
        underlying, sigma_n_fx, sigma_date_fx, sigma_fx, &corr_date_n, &corr_date, &corr);
    if (err)
    {
        goto FREE_RETURN;
    }

    today = get_today_from_underlying(und);
    for (i = 0; i < *sigma_n_fx; i++)
    {
        (*sigma_date_fx)[i] = ((*sigma_date_fx)[i] - today) / 365;
    }

    if (corr_date_n != 1)
    {
        err = "no correlation term structure allowed";
        goto FREE_RETURN;
    }

    *correl_dom_for = corr[0];
    *correl_dom_fx  = corr[1];
    *correl_for_fx  = corr[2];

FREE_RETURN:

    if (err)
    {
        if (*sigma_date_dom)
            free(*sigma_date_dom);
        *sigma_date_dom = NULL;
        if (*sigma_dom)
            free(*sigma_dom);
        *sigma_dom = NULL;

        if (*sigma_date_for)
            free(*sigma_date_for);
        *sigma_date_for = NULL;
        if (*sigma_for)
            free(*sigma_for);
        *sigma_for = NULL;

        if (*sigma_date_fx)
            free(*sigma_date_fx);
        *sigma_date_fx = NULL;
        if (*sigma_fx)
            free(*sigma_fx);
        *sigma_fx = NULL;
    }

    if (corr_date)
        free(corr_date);
    if (corr)
        free(corr);

    return err;
}

Err calibrate_correl5D_from_histo(
    double   mat, /* maturity of the long rate */
    double   lam_dom,
    double   alpha_dom,
    double   gamma_dom,
    double   rho_dom,
    int      calib_rho_dom, /* if one we calibrate rho_dom */
    double   lam_for,
    double   alpha_for,
    double   gamma_for,
    double   rho_for,
    int      calib_rho_for, /* if one we calibrate rho_for */
    double** corr_histo,    /* input of the historical correl SD, LD, SF, LF, Fx	*/
    double** corr_model)    /* ouput for correl between the brownian of the model	*/
{
    double   a, b, c, d, ex_dom, ex_for, delta, x1, x2, corr2;
    double **matrix = NULL, **matrix_inv = NULL;

    Err err = NULL;

    matrix = dmatrix(0, 3, 0, 3);

    if (!matrix)
    {
        err = "Memory allocation failure in calibrate_correl_from_histo";
        goto FREE_RETURN;
    }

    /* First step, we calibrate the rho_dom */

    ex_dom = exp(-gamma_dom * mat);

    if (calib_rho_dom)
    {
        corr2 = corr_histo[0][1] * corr_histo[0][1];

        a = alpha_dom * alpha_dom * (4 * corr2 * ex_dom - (1 + ex_dom) * (1 + ex_dom));
        b = alpha_dom * (1 + ex_dom) * (1 + alpha_dom * alpha_dom * ex_dom) * (corr2 - 1.0);
        c = corr2 * (1.0 + alpha_dom * alpha_dom) *
                (1.0 + alpha_dom * alpha_dom * ex_dom * ex_dom) -
            (1.0 + alpha_dom * alpha_dom * ex_dom) * (1.0 + alpha_dom * alpha_dom * ex_dom);

        delta = b * b - a * c;

        if (delta > -1.0E-10)
        {
            if (fabs(delta) < 1.0E-10)
            {
                delta = 0.0;
            }
            else
            {
                delta = sqrt(delta);
            }

            x1 = (-b + delta) / a;
            x2 = (-b - delta) / a;

            if (fabs(x1) < 1.0 && (x1 - x2) * corr_histo[0][1] > -1.0E-10)
            {
                rho_dom = x1;
            }
            else if (fabs(x2) < 1.0 && (x2 - x1) * corr_histo[0][1] > -1.0E-10)
            {
                rho_dom = x2;
            }
            else
            {
                err = "Cannot solve for rho domestic";
                goto FREE_RETURN;
            }
        }
        else
        {
            err = "Cannot solve for rho domestic";
            goto FREE_RETURN;
        }
    }

    corr_model[0][1] = rho_dom;

    /* Second step, solve for correl dom with S */
    a = 1.0 / sqrt(1.0 + alpha_dom * (2.0 * rho_dom + alpha_dom));
    b = alpha_dom * a;

    c = 1.0 / sqrt(1.0 + alpha_dom * ex_dom * (2.0 * rho_dom + alpha_dom * ex_dom));
    d = alpha_dom * ex_dom * c;

    delta = a * d - b * c;

    if (fabs(delta) < 1.0E-10)
    {
        err = "Cannot solve for correls between Dom and S";
        goto FREE_RETURN;
    }
    else
    {
        corr_model[0][4] = 1.0 / delta * (d * corr_histo[0][4] - b * corr_histo[1][4]);
        corr_model[1][4] = 1.0 / delta * (a * corr_histo[1][4] - c * corr_histo[0][4]);
    }

    /* Third step, we calibrate the rho_for */

    ex_for = exp(-gamma_for * mat);

    if (calib_rho_for)
    {
        corr2 = corr_histo[2][3] * corr_histo[2][3];

        a = alpha_for * alpha_for * (4 * corr2 * ex_for - (1 + ex_for) * (1 + ex_for));
        b = alpha_for * (1 + ex_for) * (1 + alpha_for * alpha_for * ex_for) * (corr2 - 1.0);
        c = corr2 * (1.0 + alpha_for * alpha_for) *
                (1.0 + alpha_for * alpha_for * ex_for * ex_for) -
            (1.0 + alpha_for * alpha_for * ex_for) * (1.0 + alpha_for * alpha_for * ex_for);

        delta = b * b - a * c;

        if (delta > -1.0E-10)
        {
            if (fabs(delta) < 1.0E-10)
            {
                delta = 0.0;
            }
            else
            {
                delta = sqrt(delta);
            }

            x1 = (-b + delta) / a;
            x2 = (-b - delta) / a;

            if (fabs(x1) < 1.0 && (x1 - x2) * corr_histo[2][3] > -1.0E-10)
            {
                rho_for = x1;
            }
            else if (fabs(x2) < 1.0 && (x2 - x1) * corr_histo[2][3] > -1.0E-10)
            {
                rho_for = x2;
            }
            else
            {
                err = "Cannot solve for rho foreign";
                goto FREE_RETURN;
            }
        }
        else
        {
            err = "Cannot solve for rho foreign";
            goto FREE_RETURN;
        }
    }

    corr_model[2][3] = rho_for;

    /* Fourth step, solve for correl dom with S */

    a = 1.0 / sqrt(1.0 + alpha_for * (2.0 * rho_for + alpha_for));
    b = alpha_for * a;

    c = 1.0 / sqrt(1.0 + alpha_for * ex_for * (2.0 * rho_for + alpha_for * ex_for));
    d = alpha_for * ex_for * c;

    delta = a * d - b * c;

    if (fabs(delta) < 1.0E-10)
    {
        err = "Cannot solve for correls between For and S";
        goto FREE_RETURN;
    }
    else
    {
        corr_model[2][4] = 1.0 / delta * (d * corr_histo[2][4] - b * corr_histo[3][4]);
        corr_model[3][4] = 1.0 / delta * (a * corr_histo[3][4] - c * corr_histo[2][4]);
    }

    /* Fifth step, solve for the others in a big system */

    a = sqrt(
        (1.0 + alpha_dom * (2.0 * rho_dom + alpha_dom)) *
        (1.0 + alpha_for * (2.0 * rho_for + alpha_for)));
    b = sqrt(
        (1.0 + alpha_dom * ex_dom * (2.0 * rho_dom + alpha_dom * ex_dom)) *
        (1.0 + alpha_for * ex_for * (2.0 * rho_for + alpha_for * ex_for)));
    c = sqrt(
        (1.0 + alpha_dom * (2.0 * rho_dom + alpha_dom)) *
        (1.0 + alpha_for * ex_for * (2.0 * rho_for + alpha_for * ex_for)));
    d = sqrt(
        (1.0 + alpha_dom * ex_dom * (2.0 * rho_dom + alpha_dom * ex_dom)) *
        (1.0 + alpha_for * (2.0 * rho_for + alpha_for)));

    /* Correlation SD and SF */
    matrix[0][0] = 1.0 / a;
    matrix[0][1] = alpha_for / a;
    matrix[0][2] = alpha_dom / a;
    matrix[0][3] = alpha_dom * alpha_for / a;

    /* Correlation LD and LF */
    matrix[1][0] = 1.0 / b;
    matrix[1][1] = alpha_for * ex_for / b;
    matrix[1][2] = alpha_dom * ex_dom / b;
    matrix[1][3] = alpha_dom * alpha_for * ex_for * ex_dom / b;

    /* Correlation SD and LF */
    matrix[2][0] = 1.0 / c;
    matrix[2][1] = alpha_for * ex_for / c;
    matrix[2][2] = alpha_dom / c;
    matrix[2][3] = alpha_dom * alpha_for * ex_for / c;

    /* Correlation LD and SF */
    matrix[3][0] = 1.0 / d;
    matrix[3][1] = alpha_for / d;
    matrix[3][2] = alpha_dom * ex_dom / d;
    matrix[3][3] = alpha_dom * alpha_for * ex_dom / d;

    matrix_inv = inverse_matrix(matrix, 0, 3);

    if (!matrix_inv)
    {
        err = "Cannot solve for cross correlations Dom / For";
        goto FREE_RETURN;
    }

    corr_model[0][2] = matrix_inv[1][1] * corr_histo[0][2] + matrix_inv[1][2] * corr_histo[1][3] +
                       matrix_inv[1][3] * corr_histo[0][3] + matrix_inv[1][4] * corr_histo[1][2];

    corr_model[0][3] = matrix_inv[2][1] * corr_histo[0][2] + matrix_inv[2][2] * corr_histo[1][3] +
                       matrix_inv[2][3] * corr_histo[0][3] + matrix_inv[2][4] * corr_histo[1][2];

    corr_model[1][2] = matrix_inv[3][1] * corr_histo[0][2] + matrix_inv[3][2] * corr_histo[1][3] +
                       matrix_inv[3][3] * corr_histo[0][3] + matrix_inv[3][4] * corr_histo[1][2];

    corr_model[1][3] = matrix_inv[4][1] * corr_histo[0][2] + matrix_inv[4][2] * corr_histo[1][3] +
                       matrix_inv[4][3] * corr_histo[0][3] + matrix_inv[4][4] * corr_histo[1][2];

FREE_RETURN:

    if (matrix)
        free_dmatrix(matrix, 0, 3, 0, 3);
    if (matrix_inv)
        free_dmatrix(matrix_inv, 1, 4, 1, 4);

    return err;
}

Err get_correlSL_from_model(
    double   mat, /* maturity of the long rate */
    double   lam_dom,
    double   alpha_dom,
    double   gamma_dom,
    double   lam_for,
    double   alpha_for,
    double   gamma_for,
    double** corr_model, /* input of the model correl  between brownian	*/
    double** corr_sd)    /* ouput for correl between SD, LD, SF, LF, Fx	*/
{
    double ex_dom, ex_for;
    double a, b, c, d;

    ex_dom = exp(-gamma_dom * mat);
    ex_for = exp(-gamma_for * mat);

    a = sqrt(1.0 + alpha_dom * (2.0 * corr_model[0][1] + alpha_dom));
    b = sqrt(1.0 + alpha_dom * ex_dom * (2.0 * corr_model[0][1] + alpha_dom * ex_dom));
    c = sqrt(1.0 + alpha_for * (2.0 * corr_model[2][3] + alpha_for));
    d = sqrt(1.0 + alpha_for * ex_for * (2.0 * corr_model[2][3] + alpha_for * ex_for));

    /* SD / LD */
    corr_sd[0][1] =
        (1.0 + corr_model[0][1] * alpha_dom * (1.0 + ex_dom) + alpha_dom * alpha_dom * ex_dom) /
        (a * b);
    /* SD / Fx */
    corr_sd[0][4] = (corr_model[0][4] + corr_model[1][4] * alpha_dom) / a;
    /* LD / Fx */
    corr_sd[1][4] = (corr_model[0][4] + corr_model[1][4] * alpha_dom * ex_dom) / b;

    /* SF / LF */
    corr_sd[2][3] =
        (1.0 + corr_model[2][3] * alpha_for * (1.0 + ex_for) + alpha_for * alpha_for * ex_for) /
        (c * d);
    /* SF / Fx */
    corr_sd[2][4] = (corr_model[2][4] + corr_model[3][4] * alpha_for) / c;
    /* LF / Fx */
    corr_sd[3][4] = (corr_model[2][4] + corr_model[3][4] * alpha_for * ex_for) / d;

    /* SD / SF */
    corr_sd[0][2] = (corr_model[0][2] + corr_model[0][3] * alpha_for +
                     corr_model[1][2] * alpha_dom + corr_model[1][3] * alpha_dom * alpha_for) /
                    (a * c);
    /* LD / LF */
    corr_sd[1][3] = (corr_model[0][2] + corr_model[0][3] * alpha_for * ex_for +
                     corr_model[1][2] * alpha_dom * ex_dom +
                     corr_model[1][3] * alpha_dom * alpha_for * ex_dom * ex_for) /
                    (b * d);
    /* SD / LF */
    corr_sd[0][3] =
        (corr_model[0][2] + corr_model[0][3] * alpha_for * ex_for + corr_model[1][2] * alpha_dom +
         corr_model[1][3] * alpha_dom * alpha_for * ex_for) /
        (a * d);
    /* LD / SF */
    corr_sd[1][2] =
        (corr_model[0][2] + corr_model[0][3] * alpha_for + corr_model[1][2] * alpha_dom * ex_dom +
         corr_model[1][3] * alpha_dom * alpha_for * ex_dom) /
        (b * c);

    return NULL;
}
