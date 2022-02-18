/**********************************************************************
 *      Name: FxOptions.c		                                      *
 *--------------------------------------------------------------------*
 *    Author: L.C.				                                      *
 *      Date: 05/25/01                                                *
 *--------------------------------------------------------------------*/

#include "FxOptions.h"

#include "BGMEval.h"
#include "CpdCalib.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "FxSabrAdi.h"
#include "FxSabrGrfn.h"
#include "FxSabrQuadAdi.h"
#include "SrtAccess.h"
#include "math.h"
#include "opfnctns.h"
#include "opsabrcalib.h"
#include "srt_h_all.h"

#define MAX_ITERNEWTON 75
#define DIGIT_SPREAD 0.005
#define GRADIENT_SHIFT 0.01

Err ConvexFx1F(
    char*   underlying,
    long    maturity,
    double  barrier,
    double  strike,
    double  typeBar, /*	0: minimum < Bar   1: maximum > bar */
    long    npaths,
    long    nstp,
    double* res1,
    double* error1,
    double* res2,
    double* error2)
{
    double *time = NULL, *date = NULL, **drift = NULL, **sig_quad = NULL;

    long today, spot_date;

    double mat, dt, dtd, cum_vol;

    double S, lnBar, lnS0, lnS;

    double      spot_fx;
    char *      dom_yc, *for_yc;
    SrtUndPtr * und_ptr = NULL, und = NULL;
    TermStruct* fx_ts;

    int i, j, flag;

    double pay1, pay2, sum1, sum2, sumsq1, sumsq2, df;

    long seed;

    Err err = NULL;

    /* look for the underlying name */
    und = lookup_und(underlying);
    if (!und)
    {
        err = "cannot find the underlying for ConvexFx1F";
        goto FREE_RETURN;
    }

    /* look for the today date */
    today     = get_today_from_underlying(und);
    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    if (get_underlying_type(und) != FOREX_UND)
    {
        err = serror("Underlying %s is not of type FX", underlying);
        goto FREE_RETURN;
    }

    if (get_mdltype_from_fxund(und) != BLACK_SCHOLES)
    {
        err = "Model must be BLACK_SCHOLES";
        goto FREE_RETURN;
    }

    /* Get the term structure */

    dom_yc  = (char*)get_domname_from_fxund(und);
    for_yc  = (char*)get_forname_from_fxund(und);
    spot_fx = get_spot_from_fxund(und) * swp_f_df(today, spot_date, dom_yc) /
              swp_f_df(today, spot_date, for_yc);

    fx_ts = get_ts_from_fxund(und);

    time = dvector(0, nstp - 1);
    date = dvector(0, nstp - 1);

    drift    = dmatrix(0, nstp - 1, 0, 1);
    sig_quad = dmatrix(0, nstp - 1, 0, 1);

    if (!time || !date || !drift || !sig_quad)
    {
        err = "Memory allocation error (1) in ConvexFx1F";
        goto FREE_RETURN;
    }

    /*	Fill the vectors */

    mat = (maturity - today) * YEARS_IN_DAY;

    dt  = mat / nstp;
    dtd = dt * DAYS_IN_YEAR;

    time[0] = 0.0;
    date[0] = today;

    for (i = 1; i < nstp; i++)
    {
        time[i] = time[i - 1] + dt;
        date[i] = date[i - 1] + dtd;

        cum_vol = fx_cum_vol_func(time[i - 1], SRT_TRUE, fx_ts);

        sig_quad[i - 1][0] = sqrt(fx_cum_vol_func(time[i], SRT_TRUE, fx_ts) - cum_vol);
        sig_quad[i - 1][1] = sqrt(fx_cum_vol_func(mat, SRT_TRUE, fx_ts) - cum_vol);

        drift[i - 1][0] =
            (swp_f_zr(date[i - 1], date[i], dom_yc) - swp_f_zr(date[i - 1], date[i], for_yc)) * dt -
            0.5 * sig_quad[i - 1][0] * sig_quad[i - 1][0];
        drift[i - 1][1] =
            (swp_f_zr(date[i - 1], maturity, dom_yc) - swp_f_zr(date[i - 1], maturity, for_yc)) *
                (mat - time[i - 1]) -
            0.5 * sig_quad[i - 1][1] * sig_quad[i - 1][1];
    }

    /*	now launch Monte Carlo Simul */

    lnBar = log(barrier);
    lnS0  = log(spot_fx);

    df = swp_f_df(today, maturity, dom_yc);

    seed = -123456789;

    sum1   = 0.0;
    sum2   = 0.0;
    sumsq1 = 0.0;
    sumsq2 = 0.0;

    if (typeBar)
    {
        for (i = 0; i < npaths; i++)
        {
            lnS  = lnS0;
            flag = 0;

            for (j = 0; j < nstp; j++)
            {
                if (lnS < lnBar)
                {
                    /* go to next time step */
                    lnS += drift[j][0] + sig_quad[j][0] * gauss_sample(&seed);
                }
                else
                {
                    /* jump to maturity */

                    lnS += drift[j][1] + sig_quad[j][1] * gauss_sample(&seed);
                    j    = nstp;
                    flag = 1;
                }
            }

            if (flag)
            {
                S = exp(lnS);

                pay1 = 1.0 / S;
                sum1 += pay1 / npaths;
                sumsq1 += pay1 * pay1 / npaths;

                if (S < strike)
                {
                    pay2 = strike - S;
                    sum2 += pay2 / npaths;
                    sumsq2 += pay2 * pay2 / npaths;
                }
            }
        }
    }
    else
    {
        for (i = 0; i < npaths; i++)
        {
            lnS  = lnS0;
            flag = 0;

            for (j = 0; j < nstp; j++)
            {
                if (lnS > lnBar)
                {
                    /* go to next time step */
                    lnS += drift[j][0] + sig_quad[j][0] * gauss_sample(&seed);
                }
                else
                {
                    /* jump to maturity */

                    lnS += drift[j][1] + sig_quad[j][1] * gauss_sample(&seed);
                    j    = nstp;
                    flag = 1;
                }
            }

            if (flag)
            {
                S = exp(lnS);

                pay1 = 1.0 / S;
                sum1 += pay1 / npaths;
                sumsq1 += pay1 * pay1 / npaths;

                if (S > strike)
                {
                    pay2 = S - strike;
                    sum2 += pay2 / npaths;
                    sumsq2 += pay2 * pay2 / npaths;
                }
            }
        }
    }

    /* save results */

    *res1   = sum1 * df * barrier;
    *error1 = sqrt((sumsq1 - sum1 * sum1) / npaths) * df * barrier;

    *res2   = sum2 * df;
    *error2 = sqrt((sumsq2 - sum2 * sum2) / npaths) * df;

FREE_RETURN:

    if (time)
        free_dvector(time, 0, nstp - 1);
    if (date)
        free_dvector(date, 0, nstp - 1);

    if (drift)
        free_dmatrix(drift, 0, nstp - 1, 0, 1);
    if (sig_quad)
        free_dmatrix(sig_quad, 0, nstp - 1, 0, 1);

    return err;
}

double DeltaPIPS(double F, double T, double DfDom, double DfFor, double K, double Vol, int IsCall)
{
    double d1, std;

    std = Vol * sqrt(T);
    d1  = log(F / K) / std + 0.5 * std;

    if (IsCall)
    {
        return DfFor * norm(d1);
    }
    else
    {
        return DfFor * (norm(d1) - 1.0);
    }
}

double DeltaPercent(
    double F, double T, double DfDom, double DfFor, double K, double Vol, int IsCall)
{
    double d2, std;

    std = Vol * sqrt(T);
    d2  = log(F / K) / std - 0.5 * std;

    if (IsCall)
    {
        return DfFor * K / F * norm(d2);
    }
    else
    {
        return DfFor * K / F * (norm(d2) - 1.0);
    }
}

double DeltaSimple(double F, double T, double DfDom, double DfFor, double K, double Vol, int IsCall)
{
    double d1, std;

    std = Vol * sqrt(T);
    d1  = log(F / K) / std;

    if (IsCall)
    {
        return DfDom * norm(d1);
    }
    else
    {
        return DfDom * (norm(d1) - 1.0);
    }
}

Err SABR_to_Strikes(
    double           Fwd,
    double           T,
    double           DfDom,
    double           DfFor,
    double           Sigma,
    double           Alpha,
    double           Beta,
    double           Rho,
    SrtDiffusionType VolType,
    int              DeltaType,
    int              IsForeign,
    double           Precision,
    double*          Strike1,
    double*          Strike2,
    double*          Strike3,
    double*          Vol1,
    double*          Vol2,
    double*          Vol3,
    double*          VolATM)
{
    Err err = NULL;

    switch (DeltaType)
    {
    case 0:
    {
        err = SABR_SolveDeltaNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            0.25,
            1,
            Sigma,
            Alpha,
            Beta,
            Rho,
            VolType,
            DeltaSimple,
            IsForeign,
            Precision,
            Strike1,
            Vol1);
        if (err)
        {
            return err;
        }

        err = SABR_SolveStraddleNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            Sigma,
            Alpha,
            Beta,
            Rho,
            VolType,
            DeltaSimple,
            IsForeign,
            Precision,
            Strike2,
            Vol2);
        if (err)
        {
            return err;
        }

        err = SABR_SolveDeltaNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            -0.25,
            0,
            Sigma,
            Alpha,
            Beta,
            Rho,
            VolType,
            DeltaSimple,
            IsForeign,
            Precision,
            Strike3,
            Vol3);
        if (err)
        {
            return err;
        }
        break;
    }
    case 1:
    {
        err = SABR_SolveDeltaNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            0.25,
            1,
            Sigma,
            Alpha,
            Beta,
            Rho,
            VolType,
            DeltaPercent,
            IsForeign,
            Precision,
            Strike1,
            Vol1);
        if (err)
        {
            return err;
        }

        err = SABR_SolveStraddleNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            Sigma,
            Alpha,
            Beta,
            Rho,
            VolType,
            DeltaPercent,
            IsForeign,
            Precision,
            Strike2,
            Vol2);
        if (err)
        {
            return err;
        }

        err = SABR_SolveDeltaNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            -0.25,
            0,
            Sigma,
            Alpha,
            Beta,
            Rho,
            VolType,
            DeltaPercent,
            IsForeign,
            Precision,
            Strike3,
            Vol3);
        if (err)
        {
            return err;
        }
        break;
    }
    case 2:
    {
        err = SABR_SolveDeltaNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            0.25,
            1,
            Sigma,
            Alpha,
            Beta,
            Rho,
            VolType,
            DeltaPIPS,
            IsForeign,
            Precision,
            Strike1,
            Vol1);
        if (err)
        {
            return err;
        }

        err = SABR_SolveStraddleNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            Sigma,
            Alpha,
            Beta,
            Rho,
            VolType,
            DeltaPIPS,
            IsForeign,
            Precision,
            Strike2,
            Vol2);
        if (err)
        {
            return err;
        }

        err = SABR_SolveDeltaNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            -0.25,
            0,
            Sigma,
            Alpha,
            Beta,
            Rho,
            VolType,
            DeltaPIPS,
            IsForeign,
            Precision,
            Strike3,
            Vol3);
        if (err)
        {
            return err;
        }
        break;
    }
    default:
    {
        return "DeltaType is 0, 1 or 2";
    }
    }

    err = srt_f_optsarbvol(Fwd, Fwd, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, VolATM);

    return err;
}

Err SABR_SolveDeltaNewton(
    double           F,
    double           T,
    double           DfDom,
    double           DfFor,
    double           Delta,
    int              IsCall,
    double           Sigma,
    double           Alpha,
    double           Beta,
    double           Rho,
    SrtDiffusionType VolType,
    double (*GetDelta)(
        double F, double T, double DfDom, double DfFor, double K, double Vol, int IsCall),
    int     IsForeign,
    double  Precision,
    double* K,
    double* Sig)
{
    int    i;
    double K1, K2, Sig1, Sig2, Delta1, Delta2, deriv, error, Finv;
    double volATM;
    Err    err = NULL;

    /* First Guess in the Simple case */

    err = srt_f_optsarbvol(F, F, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &volATM);
    if (err)
    {
        return err;
    }

    if (IsForeign)
    {
        Finv = 1.0 / F;

        if (IsCall)
        {
            K1  = Finv * exp(-inv_cumnorm_fast(Delta) * volATM * sqrt(T));
            err = srt_f_optsarbvol(
                F, 1.0 / K1, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig1);
            if (err)
                return err;
            Delta1 = GetDelta(Finv, T, DfFor, DfDom, K1, Sig1, IsCall);

            if (Delta1 > Delta)
            {
                K2 = K1 * exp(-0.25 * volATM * sqrt(T));
            }
            else
            {
                K2 = K1 * exp(0.25 * volATM * sqrt(T));
            }
        }
        else
        {
            K1  = Finv * exp(-inv_cumnorm_fast(Delta + 1.0) * volATM * sqrt(T));
            err = srt_f_optsarbvol(
                F, 1.0 / K1, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig1);
            if (err)
                return err;
            Delta1 = GetDelta(Finv, T, DfFor, DfDom, K1, Sig1, IsCall);

            if (Delta1 > Delta)
            {
                K2 = K1 * exp(0.25 * volATM * sqrt(T));
            }
            else
            {
                K2 = K1 * exp(-0.25 * volATM * sqrt(T));
            }
        }

        err = srt_f_optsarbvol(
            F, 1.0 / K2, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig2);
        if (err)
            return err;
        Delta2 = GetDelta(Finv, T, DfFor, DfDom, K2, Sig2, IsCall);

        i     = 0;
        error = fabs(Delta2 - Delta);

        while (i < MAX_ITERNEWTON && error > Precision)
        {
            deriv = (Delta2 - Delta1) / (K2 - K1);

            K1     = K2;
            Delta1 = Delta2;

            K2 += (Delta - Delta2) / deriv;
            err = srt_f_optsarbvol(
                F, 1.0 / K2, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig2);
            if (err)
                return err;
            Delta2 = GetDelta(Finv, T, DfFor, DfDom, K2, Sig2, IsCall);

            error = fabs(Delta2 - Delta);

            i++;
        }
    }
    else
    {
        if (IsCall)
        {
            K1 = F * exp(-inv_cumnorm_fast(Delta) * volATM * sqrt(T));
            err =
                srt_f_optsarbvol(F, K1, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig1);
            if (err)
                return err;
            Delta1 = GetDelta(F, T, DfDom, DfFor, K1, Sig1, IsCall);

            if (Delta1 > Delta)
            {
                K2 = K1 * exp(-0.25 * volATM * sqrt(T));
            }
            else
            {
                K2 = K1 * exp(0.25 * volATM * sqrt(T));
            }
        }
        else
        {
            K1 = F * exp(-inv_cumnorm_fast(Delta + 1.0) * volATM * sqrt(T));
            err =
                srt_f_optsarbvol(F, K1, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig1);
            if (err)
                return err;
            Delta1 = GetDelta(F, T, DfDom, DfFor, K1, Sig1, IsCall);

            if (Delta1 > Delta)
            {
                K2 = K1 * exp(0.25 * volATM * sqrt(T));
            }
            else
            {
                K2 = K1 * exp(-0.25 * volATM * sqrt(T));
            }
        }

        err = srt_f_optsarbvol(F, K2, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig2);
        if (err)
            return err;
        Delta2 = GetDelta(F, T, DfDom, DfFor, K2, Sig2, IsCall);

        i     = 0;
        error = fabs(Delta2 - Delta);

        while (i < MAX_ITERNEWTON && error > Precision)
        {
            deriv = (Delta2 - Delta1) / (K2 - K1);

            K1     = K2;
            Delta1 = Delta2;

            K2 += (Delta - Delta2) / deriv;
            err =
                srt_f_optsarbvol(F, K2, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig2);
            if (err)
                return err;
            Delta2 = GetDelta(F, T, DfDom, DfFor, K2, Sig2, IsCall);

            error = fabs(Delta2 - Delta);

            i++;
        }
    }

    if (error > Precision)
    {
        *K   = 0.0;
        *Sig = 0.0;
        err  = "Could not find solution using Newton";
    }
    else
    {
        if (IsForeign)
        {
            *K   = 1.0 / K2;
            *Sig = Sig2;
        }
        else
        {
            *K   = K2;
            *Sig = Sig2;
        }
    }

    return err;
}

Err SABR_SolveStraddleNewton(
    double           F,
    double           T,
    double           DfDom,
    double           DfFor,
    double           Sigma,
    double           Alpha,
    double           Beta,
    double           Rho,
    SrtDiffusionType VolType,
    double (*GetDelta)(
        double F, double T, double DfDom, double DfFor, double K, double Vol, int IsCall),
    int     IsForeign,
    double  Precision,
    double* K,
    double* Sig)
{
    int    i;
    double K1, K2, Sig1, Sig2, Delta1, Delta2, deriv, error, Finv;
    double volATM;
    Err    err = NULL;

    /* First Guess in the Simple case */

    err = srt_f_optsarbvol(F, F, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &volATM);
    if (err)
    {
        return err;
    }

    if (IsForeign)
    {
        Finv = 1.0 / F;

        K1  = Finv;
        err = srt_f_optsarbvol(
            F, 1.0 / K1, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig1);
        if (err)
            return err;
        Delta1 = GetDelta(Finv, T, DfFor, DfDom, K1, Sig1, 1);
        Delta1 += GetDelta(Finv, T, DfFor, DfDom, K1, Sig1, 0);
        K2 = K1 * exp(0.1 * volATM * sqrt(T));

        err = srt_f_optsarbvol(
            F, 1.0 / K2, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig2);
        if (err)
            return err;
        Delta2 = GetDelta(Finv, T, DfFor, DfDom, K2, Sig2, 1);
        Delta2 += GetDelta(Finv, T, DfFor, DfDom, K2, Sig2, 0);

        i     = 0;
        error = fabs(Delta2);

        while (i < MAX_ITERNEWTON && error > Precision)
        {
            deriv = (Delta2 - Delta1) / (K2 - K1);

            K1     = K2;
            Delta1 = Delta2;

            K2 += -Delta2 / deriv;
            err = srt_f_optsarbvol(
                F, 1.0 / K2, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig2);
            if (err)
                return err;
            Delta2 = GetDelta(Finv, T, DfFor, DfDom, K2, Sig2, 1);
            Delta2 += GetDelta(Finv, T, DfFor, DfDom, K2, Sig2, 0);

            error = fabs(Delta2);

            i++;
        }
    }
    else
    {
        K1  = F;
        err = srt_f_optsarbvol(F, K1, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig1);
        if (err)
            return err;
        Delta1 = GetDelta(F, T, DfDom, DfFor, K1, Sig1, 1);
        Delta1 += GetDelta(F, T, DfDom, DfFor, K1, Sig1, 0);
        K2 = K1 * exp(0.1 * volATM * sqrt(T));

        err = srt_f_optsarbvol(F, K2, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig2);
        if (err)
            return err;
        Delta2 = GetDelta(F, T, DfDom, DfFor, K2, Sig2, 1);
        Delta2 += GetDelta(F, T, DfDom, DfFor, K2, Sig2, 0);

        i     = 0;
        error = fabs(Delta2);

        while (i < MAX_ITERNEWTON && error > Precision)
        {
            deriv = (Delta2 - Delta1) / (K2 - K1);

            K1     = K2;
            Delta1 = Delta2;

            K2 += -Delta2 / deriv;
            err =
                srt_f_optsarbvol(F, K2, T, Sigma, Alpha, Beta, Rho, VolType, SRT_LOGNORMAL, &Sig2);
            if (err)
                return err;
            Delta2 = GetDelta(F, T, DfDom, DfFor, K2, Sig2, 1);
            Delta2 += GetDelta(F, T, DfDom, DfFor, K2, Sig2, 0);

            error = fabs(Delta2);

            i++;
        }
    }

    if (error > Precision)
    {
        *K   = 0.0;
        *Sig = 0.0;
        err  = "Could not find solution using Newton";
    }
    else
    {
        if (IsForeign)
        {
            *K   = 1.0 / K2;
            *Sig = Sig2;
        }
        else
        {
            *K   = K2;
            *Sig = Sig2;
        }
    }

    return err;
}

Err SolveDeltaNewton(
    double F,
    double T,
    double DfDom,
    double DfFor,
    double Delta,
    int    IsCall,
    double Sigma,
    double (*GetDelta)(
        double F, double T, double DfDom, double DfFor, double K, double Vol, int IsCall),
    int     IsForeign,
    double  Precision,
    double* K)
{
    int    i;
    double K1, K2, Delta1, Delta2, deriv, error, DfTemp;
    Err    err = NULL;

    if (IsForeign)
    {
        F      = 1.0 / F;
        DfTemp = DfFor;
        DfFor  = DfDom;
        DfDom  = DfTemp;
    }

    if (IsCall)
    {
        K1     = F * exp(-inv_cumnorm_fast(Delta) * Sigma * sqrt(T));
        Delta1 = GetDelta(F, T, DfDom, DfFor, K1, Sigma, IsCall);

        if (Delta1 > Delta)
        {
            K2 = K1 * exp(-0.25 * Sigma * sqrt(T));
        }
        else
        {
            K2 = K1 * exp(0.25 * Sigma * sqrt(T));
        }
    }
    else
    {
        K1     = F * exp(-inv_cumnorm_fast(Delta + 1.0) * Sigma * sqrt(T));
        Delta1 = GetDelta(F, T, DfDom, DfFor, K1, Sigma, IsCall);

        if (Delta1 > Delta)
        {
            K2 = K1 * exp(0.25 * Sigma * sqrt(T));
        }
        else
        {
            K2 = K1 * exp(-0.25 * Sigma * sqrt(T));
        }
    }

    Delta2 = GetDelta(F, T, DfDom, DfFor, K2, Sigma, IsCall);

    i     = 0;
    error = fabs(Delta2 - Delta);

    while (i < MAX_ITERNEWTON && error > Precision)
    {
        deriv = (Delta2 - Delta1) / (K2 - K1);

        K1     = K2;
        Delta1 = Delta2;

        K2 += (Delta - Delta2) / deriv;
        Delta2 = GetDelta(F, T, DfDom, DfFor, K2, Sigma, IsCall);

        error = fabs(Delta2 - Delta);

        i++;
    }

    if (error > Precision)
    {
        *K  = 0.0;
        err = "Could not find solution using Newton";
    }
    else
    {
        if (IsForeign)
        {
            *K = 1.0 / K2;
        }
        else
        {
            *K = K2;
        }
    }

    return err;
}

Err SolveStraddleNewton(
    double F,
    double T,
    double DfDom,
    double DfFor,
    double Sigma,
    double (*GetDelta)(
        double F, double T, double DfDom, double DfFor, double K, double Vol, int IsCall),
    int     IsForeign,
    double  Precision,
    double* K)
{
    int    i;
    double K1, K2, Delta1, Delta2, deriv, error, DfTemp;
    Err    err = NULL;

    if (IsForeign)
    {
        F      = 1.0 / F;
        DfTemp = DfDom;
        DfDom  = DfFor;
        DfFor  = DfTemp;
    }

    K1 = F;
    K2 = K1 * exp(-0.1 * Sigma * sqrt(T));

    Delta1 = GetDelta(F, T, DfDom, DfFor, K1, Sigma, 1);
    Delta1 += GetDelta(F, T, DfDom, DfFor, K1, Sigma, 0);
    Delta2 = GetDelta(F, T, DfDom, DfFor, K2, Sigma, 1);
    Delta2 += GetDelta(F, T, DfDom, DfFor, K2, Sigma, 0);

    i     = 0;
    error = fabs(Delta2);

    while (i < MAX_ITERNEWTON && error > Precision)
    {
        deriv = (Delta2 - Delta1) / (K2 - K1);

        K1     = K2;
        Delta1 = Delta2;

        K2 += -Delta2 / deriv;
        Delta2 = GetDelta(F, T, DfDom, DfFor, K2, Sigma, 1);
        Delta2 += GetDelta(F, T, DfDom, DfFor, K2, Sigma, 0);

        error = fabs(Delta2);

        i++;
    }

    if (error > Precision)
    {
        *K  = 0.0;
        err = "Could not find solution using Newton";
    }
    else
    {
        if (IsForeign)
        {
            *K = 1.0 / K2;
        }
        else
        {
            *K = K2;
        }
    }

    return err;
}

Err FindStrikes(
    double  Fwd,
    double  T,
    double  DfDom,
    double  DfFor,
    int     DeltaType,
    int     IsForeign,
    double  Precision,
    double  Vol25Put,
    double  VolStraddle,
    double  Vol25Call,
    double* Strike25Put,
    double* StrikeStraddle,
    double* Strike25Call)
{
    Err err = NULL;

    switch (DeltaType)
    {
    case 0:
    {
        err = SolveStraddleNewton(
            Fwd, T, DfDom, DfFor, VolStraddle, DeltaSimple, IsForeign, Precision, StrikeStraddle);
        if (err)
        {
            return err;
        }
        err = SolveDeltaNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            0.25,
            1,
            Vol25Call,
            DeltaSimple,
            IsForeign,
            Precision,
            Strike25Call);
        if (err)
        {
            return err;
        }
        err = SolveDeltaNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            -0.25,
            0,
            Vol25Put,
            DeltaSimple,
            IsForeign,
            Precision,
            Strike25Put);
        if (err)
        {
            return err;
        }
        break;
    }
    case 1:
    {
        err = SolveStraddleNewton(
            Fwd, T, DfDom, DfFor, VolStraddle, DeltaPercent, IsForeign, Precision, StrikeStraddle);
        if (err)
        {
            return err;
        }
        err = SolveDeltaNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            0.25,
            1,
            Vol25Call,
            DeltaPercent,
            IsForeign,
            Precision,
            Strike25Call);
        if (err)
        {
            return err;
        }
        err = SolveDeltaNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            -0.25,
            0,
            Vol25Put,
            DeltaPercent,
            IsForeign,
            Precision,
            Strike25Put);
        if (err)
        {
            return err;
        }
        break;
    }
    case 2:
    {
        err = SolveStraddleNewton(
            Fwd, T, DfDom, DfFor, VolStraddle, DeltaPIPS, IsForeign, Precision, StrikeStraddle);
        if (err)
        {
            return err;
        }
        err = SolveDeltaNewton(
            Fwd,
            T,
            DfDom,
            DfFor,
            0.25,
            1,
            Vol25Call,
            DeltaPIPS,
            IsForeign,
            Precision,
            Strike25Call);
        if (err)
        {
            return err;
        }
        err = SolveDeltaNewton(
            Fwd, T, DfDom, DfFor, -0.25, 0, Vol25Put, DeltaPIPS, IsForeign, Precision, Strike25Put);
        if (err)
        {
            return err;
        }
        break;
    }
    default:
    {
        return "DeltaType is 0, 1 or 2";
    }
    }

    return err;
}

Err Strikes_to_SABR(
    double  Fwd,
    double  T,
    double  DfDom,
    double  DfFor,
    double  Vol25Put,
    double  VolStraddle,
    double  Vol25Call,
    double* Strike25Put,
    double* StrikeStraddle,
    double* Strike25Call,
    double* volATM,
    double* Alpha,
    double* Beta,
    double* Rho,
    int     FreezeAlpha,
    int     FreezeBeta,
    int     FreezeRho,
    int     DeltaType,
    int     IsForeign,
    double  Precision,
    double* fitting_error)
{
    Err     err         = NULL;
    double *strike_sabr = NULL, *vol_sabr = NULL;

    err = FindStrikes(
        Fwd,
        T,
        DfDom,
        DfFor,
        DeltaType,
        IsForeign,
        Precision,
        Vol25Put,
        VolStraddle,
        Vol25Call,
        Strike25Put,
        StrikeStraddle,
        Strike25Call);

    if (err)
    {
        return err;
    }

    strike_sabr = dvector(1, 3);
    vol_sabr    = dvector(1, 3);

    if (!strike_sabr || !vol_sabr)
    {
        if (strike_sabr)
            free_dvector(strike_sabr, 1, 3);
        if (vol_sabr)
            free_dvector(vol_sabr, 1, 3);
        return "Memory allocation failure";
    }

    if (*Strike25Call < *StrikeStraddle)
    {
        if (*StrikeStraddle < *Strike25Put)
        {
            strike_sabr[1] = *Strike25Call;
            strike_sabr[2] = *StrikeStraddle;
            strike_sabr[3] = *Strike25Put;

            vol_sabr[1] = Vol25Call;
            vol_sabr[2] = VolStraddle;
            vol_sabr[3] = Vol25Put;
        }
        else if (*Strike25Call < *Strike25Put)
        {
            strike_sabr[1] = *Strike25Call;
            strike_sabr[2] = *Strike25Put;
            strike_sabr[3] = *StrikeStraddle;

            vol_sabr[1] = Vol25Call;
            vol_sabr[2] = Vol25Put;
            vol_sabr[3] = VolStraddle;
        }
        else
        {
            strike_sabr[1] = *Strike25Put;
            strike_sabr[2] = *Strike25Call;
            strike_sabr[3] = *StrikeStraddle;

            vol_sabr[1] = Vol25Put;
            vol_sabr[2] = Vol25Call;
            vol_sabr[3] = VolStraddle;
        }
    }
    else
    {
        if (*Strike25Call < *Strike25Put)
        {
            strike_sabr[1] = *StrikeStraddle;
            strike_sabr[2] = *Strike25Call;
            strike_sabr[3] = *Strike25Put;

            vol_sabr[1] = VolStraddle;
            vol_sabr[2] = Vol25Call;
            vol_sabr[3] = Vol25Put;
        }
        else if (*StrikeStraddle < *Strike25Put)
        {
            strike_sabr[1] = *StrikeStraddle;
            strike_sabr[2] = *Strike25Put;
            strike_sabr[3] = *Strike25Call;

            vol_sabr[1] = VolStraddle;
            vol_sabr[2] = Vol25Put;
            vol_sabr[3] = Vol25Call;
        }
        else
        {
            strike_sabr[1] = *Strike25Put;
            strike_sabr[2] = *StrikeStraddle;
            strike_sabr[3] = *Strike25Call;

            vol_sabr[1] = Vol25Put;
            vol_sabr[2] = VolStraddle;
            vol_sabr[3] = Vol25Call;
        }
    }

    *volATM = VolStraddle;

    err = opsabrcalib(
        Fwd,
        T,
        3,
        strike_sabr,
        vol_sabr,
        volATM,
        Alpha,
        FreezeAlpha,
        Beta,
        FreezeBeta,
        Rho,
        FreezeRho,
        fitting_error);

    return err;
}

double DeltaSoho(
    double Fwd,
    double T,
    double DfDom,
    double DfFor,
    double K,
    double Vol,
    int    IsCall,
    int    DeltaType)
{
    switch (DeltaType)
    {
    case 0:
    {
        return DeltaSimple(Fwd, T, DfDom, DfFor, K, Vol, IsCall);
    }
    case 1:
    {
        return DeltaPercent(Fwd, T, DfDom, DfFor, K, Vol, IsCall);
    }
    case 2:
    {
        return DeltaPIPS(Fwd, T, DfDom, DfFor, K, Vol, IsCall);
    }
    default:
    {
        return 0.0;
    }
    }
}

Err BarrierOption_autocal(
    /*	Market parameters */
    long   today,
    char*  dom_yc,
    char*  for_yc,
    double spot_fx, /* as quoted in the market */

    /*	Product parameters */
    double notional,
    long   exercise_date,
    long   settlmt_date,
    double strike,
    double barrier_up,
    double barrier_down,
    double rebate_up,
    double rebate_down,
    int    is_call,     /* 1: Call, 0: Put */
    int    is_ko,       /* 1: KO, 0: KI */
    int    is_cvx,      /* 1: use 1/Fx, 0: use Fx */
    int    is_digital,  /* 1: digital payoff, 0, regular option payoff */
    int    is_american, /* 1: American, 0: European */

    /*	Model parameters */
    double mdl_alpha,
    double mdl_beta,
    double mdl_rho,
    double mdl_lambda,
    double mdl_gamma,
    double mdl_floor,

    /*	Calib parameters */
    long*   opt_matdates, /* mkt options maturity */
    double* opt_atmvols,  /* mkt options ATM vols */
    int     nb_opt,
    int     calib_all_ts, /* 1: calib all the TS, 0: calib only the relevant */

    int    do_smile_calib, /* 1: calibrates alpha and rho, 0: uses model input */
    int    cal_stochvol,   /* 1: calib alpha and rho, 0: calib beta and cvxty */
    long   cal_maturity,   /* option maturity, exercise is 2bd before */
    double cal_alpha,      /* SABR alpha */
    double cal_beta,       /* SABR beta */
    double cal_rho,        /* SABR rho */

    /*	Numerical parameters */
    int pr_nstept,   /* time steps for pricing */
    int pr_nstepfx,  /* Fx steps for pricing */
    int pr_nstepvol, /* Vol steps for pricing */

    int    cal_nstept,       /* time steps for calibration */
    int    cal_nstepfx,      /* Fx steps for calibration */
    int    cal_nstepvol,     /* Vol steps for calibration */
    int    cal_atmmaxiter,   /* Maximum number of iterations for ATM calibration */
    double cal_atmprec,      /* Precision on ATM vols */
    int    cal_smilemaxiter, /* Maximum number of iterations for ATM calibration */
    double cal_smileprec,    /* Precision on ATM vols */
    int    cal_usetotalts,   /* 1: when calibrating smile, consider total TS, 0: consider only the
                                corresponding ATM vol */

    /* IOD / EOD flags */
    int eod_fix_flag, /*	EOD Fixing Flag 0: I, 1: E */
    int eod_pay_flag, /*	EOD Payment Flag 0: I, 1: E */
    int eod_ex_flag,  /*	EOD Exercise Flag 0: I, 1: E */

    /*	Exercised flag */
    int    exercised,      /*	Is exercised Flag */
    int    knocked,        /*	Has knocked flag */
    int    has_knocked_up, /*	Knocked on the up barrier:1, on the down one: 0 */
    long   ex_date_ex,     /*	Date when exercised */
    long   ex_date_set,    /*	Corresponding settlement date */
    double spot_fx_fixing, /*	Fixing of the Fx */

    /*	Outputs */
    double*        price,
    double*        greeks,
    int            export_ts,
    int            use_old_ts,
    BARAUTOCAL_UND calib_und)
{
    Err  err = NULL;
    long spot_date;

    double *opt_mat = NULL, *opt_exe = NULL, *calib_ts = NULL;
    double *time = NULL, *sig = NULL, *drift = NULL, *date = NULL;
    double *bar_lvl_up = NULL, *bar_lvl_down = NULL, *prod_val = NULL, *prod_val_ki = NULL,
           *greeks_ki = NULL;

    long   exercise_date_temp, date_temp;
    double cal_fwd, cal_mat, cal_dfdom, cal_dffor, cal_strike1, cal_strike2, cal_strike3, cal_sigma,
        cal_vol1, cal_vol2, cal_vol3, cal_volATM;
    long   cal_exercise, cal_index;
    double cash_fx, fwd_fx;
    double bar_mat;
    int    i, ns, index, nstp, num_col;
    double coef, a, b, c;

    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
    cash_fx   = spot_fx * swp_f_df(today, spot_date, dom_yc) / swp_f_df(today, spot_date, for_yc);

    if (export_ts && !use_old_ts)
    {
        calib_und->nb_sig = 0;
    }

    /* step 0: check exercised flag */
    if (exercised)
    {
        if (today + eod_pay_flag > ex_date_set)
        {
            /* already settled */
            *price = 0.0;
        }
        else
        {
            /* not yet settled */
            if (today + eod_pay_flag > ex_date_set)
            {
                /* already fixed */
                if (is_cvx)
                {
                    if (is_call)
                    {
                        *price = max(spot_fx_fixing - strike, 0.0) *
                                 swp_f_df(today, ex_date_set, dom_yc) * notional;
                    }
                    else
                    {
                        *price = max(strike - spot_fx_fixing, 0.0) *
                                 swp_f_df(today, ex_date_set, dom_yc) * notional;
                    }
                }
                else
                {
                    if (is_call)
                    {
                        *price = max(1.0 / spot_fx_fixing - strike, 0.0) *
                                 swp_f_df(today, ex_date_set, dom_yc) * notional;
                    }
                    else
                    {
                        *price = max(strike - 1.0 / spot_fx_fixing, 0.0) *
                                 swp_f_df(today, ex_date_set, dom_yc) * notional;
                    }
                }
            }
            else
            {
                /* fx not yet fixed */
                fwd_fx = cash_fx * swp_f_df(today, ex_date_set, for_yc) /
                         swp_f_df(today, ex_date_set, dom_yc);

                if (is_cvx)
                {
                    /*	choto wrong */
                    if (is_call)
                    {
                        *price = (1.0 / fwd_fx - strike) * swp_f_df(today, ex_date_set, dom_yc) *
                                 notional;
                    }
                    else
                    {
                        *price = (strike - 1.0 / fwd_fx) * swp_f_df(today, ex_date_set, dom_yc) *
                                 notional;
                    }
                }
                else
                {
                    if (is_call)
                    {
                        *price =
                            (fwd_fx - strike) * swp_f_df(today, ex_date_set, dom_yc) * notional;
                    }
                    else
                    {
                        *price =
                            (strike - fwd_fx) * swp_f_df(today, ex_date_set, dom_yc) * notional;
                    }
                }
            }
        }

        goto FREE_RETURN;
    }

    if (knocked)
    {
        if (is_ko)
        {
            if (today + eod_pay_flag > ex_date_set)
            {
                /* already settled */
                *price = 0.0;
            }
            else
            {
                if (has_knocked_up)
                {
                    *price = rebate_up * swp_f_df(today, ex_date_set, dom_yc) * notional;
                }
                else
                {
                    *price = rebate_down * swp_f_df(today, ex_date_set, dom_yc) * notional;
                }
            }
        }
        else
        {
            err = BarrierOption_autocal(
                today,
                dom_yc,
                for_yc,
                spot_fx,
                notional,
                exercise_date,
                settlmt_date,
                strike,
                1000000 * spot_fx,
                -1000000 * spot_fx,
                rebate_up,
                rebate_down,
                is_call,
                1,
                is_cvx,
                is_digital,
                is_american,
                mdl_alpha,
                mdl_beta,
                mdl_rho,
                mdl_lambda,
                mdl_gamma,
                mdl_floor,
                opt_matdates,
                opt_atmvols,
                nb_opt,
                calib_all_ts,
                do_smile_calib,
                cal_stochvol,
                cal_maturity,
                cal_alpha,
                cal_beta,
                cal_rho,
                pr_nstept,
                pr_nstepfx,
                pr_nstepvol,
                cal_nstept,
                cal_nstepfx,
                cal_nstepvol,
                cal_atmmaxiter,
                cal_atmprec,
                cal_smilemaxiter,
                cal_smileprec,
                cal_usetotalts,
                eod_fix_flag,
                eod_pay_flag,
                eod_ex_flag,
                exercised,
                0,
                has_knocked_up,
                ex_date_ex,
                ex_date_set,
                spot_fx_fixing,
                price,
                greeks,
                export_ts,
                use_old_ts,
                calib_und);
        }

        goto FREE_RETURN;
    }

    /* step 1: model calibration */

    if (!use_old_ts)
    {
        /* get the relevant number of options */
        i = 0;
        while (i < nb_opt && opt_matdates[i] < settlmt_date)
        {
            i++;
        }
        if (i < nb_opt)
        {
            nb_opt = i + 1;
        }

        if (!calib_all_ts && nb_opt > 1)
        {
            if (opt_matdates[nb_opt - 1] > settlmt_date)
            {
                opt_matdates = &(opt_matdates[nb_opt - 2]);
                opt_atmvols  = &(opt_atmvols[nb_opt - 2]);
                nb_opt       = 2;
            }
            else
            {
                opt_matdates = &(opt_matdates[nb_opt - 1]);
                opt_atmvols  = &(opt_atmvols[nb_opt - 1]);
                nb_opt       = 1;
            }
        }

        opt_mat  = dvector(0, nb_opt - 1);
        opt_exe  = dvector(0, nb_opt - 1);
        calib_ts = dvector(0, nb_opt - 1);

        if (!opt_mat || !opt_exe || !calib_ts)
        {
            err = "Memory allocation faillure (1) in BarrierOption_autocal";
            goto FREE_RETURN;
        }

        for (i = 0; i < nb_opt; i++)
        {
            exercise_date_temp = add_unit(opt_matdates[i], -2, SRT_BDAY, MODIFIED_SUCCEEDING);
            opt_mat[i]         = (opt_matdates[i] - today) * YEARS_IN_DAY;
            opt_exe[i]         = (exercise_date_temp - today) * YEARS_IN_DAY;
        }

        if (do_smile_calib)
        {
            /* get two strikes and two vols from SABR */

            cal_exercise = add_unit(cal_maturity, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
            cal_mat      = (cal_exercise - today) * YEARS_IN_DAY;
            cal_dfdom    = swp_f_df(today, cal_maturity, dom_yc);
            cal_dffor    = swp_f_df(today, cal_maturity, for_yc);
            cal_fwd      = cash_fx * cal_dffor / cal_dfdom;

            cal_index = 0;
            while (fabs(opt_matdates[cal_index] - cal_maturity) > 1.0E-08 && cal_index < nb_opt)
            {
                cal_index++;
            }

            if (cal_index == nb_opt)
            {
                err = serror(
                    "Maturity %d for calibration is not present in the market ATM term struct",
                    cal_maturity);
                goto FREE_RETURN;
            }

            cal_sigma = opt_atmvols[cal_index];

            err = SABR_to_Strikes(
                cal_fwd,
                cal_mat,
                1.0,
                1.0,
                cal_sigma,
                cal_alpha,
                cal_beta,
                cal_rho,
                SRT_LOGNORMAL,
                0,
                0,
                0.0001,
                &cal_strike1,
                &cal_strike2,
                &cal_strike3,
                &cal_vol1,
                &cal_vol2,
                &cal_vol3,
                &cal_volATM);

            if (err)
            {
                goto FREE_RETURN;
            }
        }

        if ((fabs(mdl_alpha) < 1.0E-04 &&
             (do_smile_calib == 0 || (do_smile_calib == 1 && cal_stochvol == 0))) ||
            (do_smile_calib == 1 && fabs(cal_alpha) < 1.0E-04))
        {
            mdl_alpha    = 1.0E-08;
            pr_nstepvol  = 1;
            cal_nstepvol = 1;
        }

        if (do_smile_calib == 1 && cal_stochvol == 0 && mdl_gamma < 1.0E-08)
        {
            mdl_gamma = max(cal_alpha - mdl_alpha, 0.0001) / 100;
        }

        if (do_smile_calib == 1 && cal_stochvol == 0 && mdl_beta < 1.0E-08)
        {
            mdl_beta = cal_beta;
        }

        if (do_smile_calib == 1 && cal_stochvol == 1 && mdl_alpha < 1.0E-04)
        {
            mdl_alpha = cal_alpha;
            mdl_rho   = cal_rho;
        }

        err = FxSabrQuadCalibSmile2(
            dom_yc,
            for_yc,
            today,
            spot_fx,
            &mdl_alpha,
            &mdl_gamma,
            &mdl_beta,
            &mdl_rho,
            mdl_lambda,
            mdl_floor,
            do_smile_calib,
            1,
            cal_usetotalts,
            cal_stochvol,
            opt_exe,
            opt_mat,
            opt_atmvols,
            nb_opt,
            0,
            cal_maturity,
            0.0,
            0.0,
            cal_nstept,
            cal_nstepfx,
            cal_nstepvol,
            cal_atmmaxiter,
            cal_atmprec,
            cal_smilemaxiter,
            cal_smileprec,
            &cal_strike1,
            &cal_vol1,
            &cal_strike3,
            &cal_vol3,
            calib_ts);
        if (err)
        {
            goto FREE_RETURN;
        }
    }
    else
    {
        mdl_alpha  = calib_und->alpha;
        mdl_gamma  = calib_und->cvxty;
        mdl_beta   = calib_und->beta;
        mdl_rho    = calib_und->rho;
        mdl_lambda = calib_und->lambda;
        mdl_floor  = calib_und->floormu;

        nb_opt   = calib_und->nb_sig;
        opt_exe  = dvector(0, nb_opt - 1);
        calib_ts = dvector(0, nb_opt - 1);

        if (!opt_exe || !calib_ts)
        {
            err = "Memory allocation faillure (1) in BarrierOption_autocal";
            goto FREE_RETURN;
        }

        memcpy(calib_ts, calib_und->sig, calib_und->nb_sig * sizeof(double));
        memcpy(opt_exe, calib_und->sig_time, calib_und->nb_sig * sizeof(double));
    }

    if (fabs(mdl_alpha) < 1.0E-04)
    {
        mdl_alpha   = 1.0E-08;
        pr_nstepvol = 1;
    }

    /* Step2: price the deal */

    bar_mat = (exercise_date - today) * YEARS_IN_DAY;

    /* discretise in time */
    ns   = 2;
    time = (double*)calloc(ns, sizeof(double));
    if (!time)
    {
        err = "Memory allocation error (2) in BarrierOption_autocal";
        goto FREE_RETURN;
    }

    time[0] = 0.0;
    time[1] = bar_mat;

    /*	Add revelant vol times	*/

    for (i = 0; i < nb_opt; i++)
    {
        if (opt_exe[i] < bar_mat)
        {
            num_f_add_number(&ns, &time, opt_exe[i]);
        }
        else
        {
            break;
        }
    }

    num_f_sort_vector(ns, time);
    num_f_unique_vector(&ns, time);

    /*	Fill the time vector */
    err = fill_time_vector(&time, &ns, 0, NULL, 0, NULL, pr_nstept);
    if (err)
    {
        goto FREE_RETURN;
    }

    nstp = ns;

    date  = (double*)calloc(nstp, sizeof(double));
    sig   = (double*)calloc(nstp, sizeof(double));
    drift = (double*)calloc(nstp, sizeof(double));

    bar_lvl_up   = (double*)calloc(nstp, sizeof(double));
    bar_lvl_down = (double*)calloc(nstp, sizeof(double));

    if (!sig || !date || !drift || !bar_lvl_up || !bar_lvl_down)
    {
        err = "Memory allocation failure (3) in BarrierOption_autocal";
        goto FREE_RETURN;
    }

    /* first part */

    i     = 0;
    index = 0;

    while (i < nstp && time[i] < opt_exe[index] - 1.0E-08)
    {
        drift[i] = 0.0;
        i++;
    }

    /* middle and end part */
    while (i < nstp)
    {
        if (index < nb_opt - 1)
        {
            /* middle part */
            coef =
                log(calib_ts[index + 1] / calib_ts[index]) / (opt_exe[index + 1] - opt_exe[index]);

            index++;

            while (i < nstp && time[i] < opt_exe[index] - 1.0E-08)
            {
                drift[i] = coef;
                i++;
            }
        }
        else
        {
            /* end part */
            coef = 0.0;

            while (i < nstp)
            {
                drift[i] = coef;
                i++;
            }
        }
    }

    for (i = nstp - 1; i >= 0; i--)
    {
        date[i]   = today + time[i] * DAYS_IN_YEAR;
        date_temp = add_unit((long)(date[i] + 1.0E-08), 2, SRT_BDAY, MODIFIED_SUCCEEDING);
        bar_lvl_up[i] =
            barrier_up * swp_f_df(today, date_temp, dom_yc) / swp_f_df(today, date_temp, for_yc);
        bar_lvl_down[i] =
            barrier_down * swp_f_df(today, date_temp, dom_yc) / swp_f_df(today, date_temp, for_yc);
    }

    /*	Eventually! call to function */

    /* Transformation from beta and convexity to b and c */
    SabrQuadGetParams(cash_fx, mdl_gamma, mdl_beta, &a, &b, &c);

    num_col  = 1;
    prod_val = dvector(0, num_col - 1);

    if (!prod_val)
    {
        err = "Memory allocation failure (4) in BarrierOption_autocal";
        goto FREE_RETURN;
    }

    if (((cash_fx > barrier_up) || (cash_fx < barrier_down)) && eod_fix_flag == 0)
    {
        prod_val[0] = 0.0;
    }
    else
    {
        err = FxSabrQuad_KOOption(
            nstp,
            time,
            date,
            pr_nstepfx,
            pr_nstepvol,
            calib_ts[0],
            drift,
            mdl_alpha,
            a,
            b,
            c,
            mdl_rho,
            mdl_lambda,
            mdl_floor,
            settlmt_date,
            strike,
            is_call,
            is_american,
            is_cvx,
            is_digital,
            bar_lvl_up,
            bar_lvl_down,
            rebate_up,
            rebate_down,
            cash_fx,
            dom_yc,
            for_yc,
            eod_fix_flag,
            eod_pay_flag,
            eod_ex_flag,
            prod_val,
            1,
            greeks);
    }

    if (err)
    {
        goto FREE_RETURN;
    }

    if (!is_ko)
    {
        /* now we need to price the regular option */

        greeks_ki   = dvector(0, 5);
        prod_val_ki = dvector(0, num_col - 1);

        if (!greeks_ki || !prod_val_ki)
        {
            err = "Memory allocation failure (3) in BarrierOption_autocal";
            goto FREE_RETURN;
        }

        for (i = nstp - 1; i >= 0; i--)
        {
            bar_lvl_up[i]   = 1.0E09;
            bar_lvl_down[i] = -1.0E09;
        }

        err = FxSabrQuad_KOOption(
            nstp,
            time,
            date,
            pr_nstepfx,
            pr_nstepvol,
            calib_ts[0],
            drift,
            mdl_alpha,
            a,
            b,
            c,
            mdl_rho,
            mdl_lambda,
            mdl_floor,
            settlmt_date,
            strike,
            is_call,
            is_american,
            is_cvx,
            is_digital,
            bar_lvl_up,
            bar_lvl_down,
            0.0,
            0.0,
            cash_fx,
            dom_yc,
            for_yc,
            eod_fix_flag,
            eod_pay_flag,
            eod_ex_flag,
            prod_val_ki,
            1,
            greeks_ki);

        if (err)
        {
            goto FREE_RETURN;
        }

        *price = (prod_val_ki[0] - prod_val[0]) * notional;

        greeks[0] = (greeks_ki[0] - greeks[0]) * notional;
        greeks[1] = (greeks_ki[1] - greeks[1]) * notional;
        greeks[2] = (greeks_ki[2] - greeks[2]) * notional;
        greeks[3] = (greeks_ki[3] - greeks[3]) * notional;
        greeks[4] = (greeks_ki[4] - greeks[4]) * notional;
        greeks[5] = (greeks_ki[5] - greeks[5]) * notional;
    }
    else
    {
        *price = prod_val[0] * notional;
        greeks[0] *= notional;
        greeks[1] *= notional;
        greeks[2] *= notional;
        greeks[3] *= notional;
        greeks[4] *= notional;
        greeks[5] *= notional;
    }

    /* Step 3: export the term structure */

    if (export_ts && !use_old_ts)
    {
        calib_und->today    = today;
        calib_und->sig_time = dvector(0, nb_opt - 1);
        calib_und->sig      = dvector(0, nb_opt - 1);

        if (!calib_und->sig || !calib_und->sig_time)
        {
            err = "Memory allocation failure (4) in BarrierOption_autocal";
            goto FREE_RETURN;
        }

        for (i = 0; i < nb_opt; i++)
        {
            calib_und->sig_time[i] = opt_exe[i];
            calib_und->sig[i]      = calib_ts[i];
        }

        calib_und->spot_fx = spot_fx;
        calib_und->nb_sig  = nb_opt;
        calib_und->alpha   = mdl_alpha;
        calib_und->beta    = mdl_beta;
        calib_und->rho     = mdl_rho;
        calib_und->lambda  = mdl_lambda;
        calib_und->cvxty   = mdl_gamma;
        calib_und->floormu = mdl_floor;
    }

FREE_RETURN:

    if (opt_mat)
        free_dvector(opt_mat, 0, nb_opt - 1);
    if (opt_exe)
        free_dvector(opt_exe, 0, nb_opt - 1);
    if (calib_ts)
        free_dvector(calib_ts, 0, nb_opt - 1);

    if (date)
        free(date);
    if (time)
        free(time);
    if (sig)
        free(sig);
    if (drift)
        free(drift);

    if (bar_lvl_up)
        free(bar_lvl_up);
    if (bar_lvl_down)
        free(bar_lvl_down);
    if (prod_val)
        free_dvector(prod_val, 0, num_col - 1);

    if (greeks_ki)
        free_dvector(greeks_ki, 0, 5);
    if (prod_val_ki)
        free_dvector(prod_val_ki, 0, num_col - 1);

    return err;
}

void barautocal_free_und(BARAUTOCAL_UND und)
{
    if (und->sig_time)
    {
        free_dvector(und->sig_time, 0, und->nb_sig - 1);
        und->sig_time = NULL;
    }

    if (und->sig)
    {
        free_dvector(und->sig, 0, und->nb_sig - 1);
        und->sig = NULL;
    }

    und->nb_sig = 0;

    und = NULL;
}

Err FxSabrQuadSmilets(
    long    today,
    double* sigma_time_fx,
    double* sigma_fx,
    long    sigma_n_fx,
    char*   dom_yc,
    char*   for_yc,
    double  spot_fx,
    double  alpha,
    double  gamma,
    double  beta,
    double  rho,
    double  lambda,
    double  floormu,
    /*	Product data */
    long     mat_date,
    double*  strike,
    int      nb_strike,
    int      nstp,
    int      nstpfx,
    int      nstpvol,
    double** res)
{
    double *time = NULL, *sig = NULL, *drift = NULL, *date = NULL;

    long spot_date, exe_date;

    int num_col, ns;

    double* prod_val = NULL;

    long index;

    double coef, fwd_fx, df, maturity, bs_vol;

    SrtUndPtr *und_ptr = NULL, und = NULL, dom_und = NULL, for_und = NULL;

    double** func_parm = NULL;

    int* is_event = NULL;

    int i;

    double a, b, c;

    Err err = NULL;

    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    spot_fx *= swp_f_df(today, spot_date, dom_yc) / swp_f_df(today, spot_date, for_yc);

    exe_date = add_unit(mat_date, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
    maturity = (exe_date - today) * YEARS_IN_DAY;

    fwd_fx = get_spot_from_fxund(und) * swp_f_df(spot_date, mat_date, for_yc) /
             swp_f_df(spot_date, mat_date, dom_yc);

    df = swp_f_df(today, exe_date, dom_yc);

    /* discretise in time			*/

    ns = 2;

    time = (double*)calloc(ns, sizeof(double));
    if (!time)
    {
        err = "Memory allocation error (1) in SrtGrfnFxSabrAdiBar";
        goto FREE_RETURN;
    }

    time[0] = 0.0;
    time[1] = maturity;

    /*	Add revelant vol times	*/

    for (i = 0; i < sigma_n_fx; i++)
    {
        if (sigma_time_fx[i] < maturity)
        {
            num_f_add_number(&ns, &time, sigma_time_fx[i]);
        }
        else
        {
            break;
        }
    }

    num_f_sort_vector(ns, time);
    num_f_unique_vector(&ns, time);

    /*	Fill the time vector */
    err = fill_time_vector(&time, &ns, 0, NULL, 0, NULL, nstp);
    if (err)
    {
        goto FREE_RETURN;
    }

    nstp = ns;

    date  = (double*)calloc(nstp, sizeof(double));
    sig   = (double*)calloc(nstp, sizeof(double));
    drift = (double*)calloc(nstp, sizeof(double));

    if (!sig || !date || !drift)
    {
        err = "Memory allocation failure in SrtGrfnFxSabrAdiBar";
        goto FREE_RETURN;
    }

    /* first part */

    i     = 0;
    index = 0;

    while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
    {
        drift[i] = 0.0;
        i++;
    }

    /* middle and end part */
    while (i < nstp)
    {
        if (index < sigma_n_fx - 1)
        {
            /* middle part */
            coef = log(sigma_fx[index + 1] / sigma_fx[index]) /
                   (sigma_time_fx[index + 1] - sigma_time_fx[index]);

            index++;

            while (i < nstp && time[i] < sigma_time_fx[index] - 1.0E-08)
            {
                drift[i] = coef;
                i++;
            }
        }
        else
        {
            /* end part */
            coef = 0.0;

            while (i < nstp)
            {
                drift[i] = coef;
                i++;
            }
        }
    }

    for (i = nstp - 1; i >= 0; i--)
    {
        date[i] = today + time[i] * DAYS_IN_YEAR;
    }

    /*	Eventually! call to function */

    num_col = nb_strike;

    /* now we need to price the regular option */

    func_parm = dmatrix(0, nstp - 1, 0, num_col - 1);
    *res      = dvector(0, nb_strike - 1);
    is_event  = calloc(nstp, sizeof(int));

    if (!func_parm || !is_event || !(*res))
    {
        err = "Memory allocation faillure in FxSabrAdiKO";
        goto FREE_RETURN;
    }

    for (i = 0; i < num_col; i++)
    {
        func_parm[nstp - 1][i] = strike[i];
    }

    is_event[nstp - 1] = 1;

    /* Transformation from beta and convexity to b and c */
    SabrQuadGetParams(spot_fx, gamma, beta, &a, &b, &c);

    err = FxSabrQuad_adi(
        nstp,
        time,
        date,
        nstpfx,
        nstpvol,
        sigma_fx[0],
        drift,
        alpha,
        a,
        b,
        c,
        rho,
        lambda,
        floormu,
        func_parm,
        is_event,
        NULL,
        NULL,
        NULL,
        spot_fx,
        dom_yc,
        for_yc,
        payoff_fx_sabr_adi_opt,
        num_col,
        *res,
        0,
        NULL,
        0,
        0,
        0,
        NULL,
        NULL);

    if (err)
    {
        goto FREE_RETURN;
    }

    for (i = 0; i < nb_strike; i++)
    {
        err = srt_f_optimpvol(
            (*res)[i], fwd_fx, strike[i], maturity, df, SRT_CALL, SRT_LOGNORMAL, &bs_vol);

        (*res)[i] = bs_vol;
    }

FREE_RETURN:

    if (date)
        free(date);
    if (time)
        free(time);
    if (sig)
        free(sig);
    if (drift)
        free(drift);

    if (func_parm)
        free_dmatrix(func_parm, 0, nstp - 1, 0, num_col - 1);
    if (is_event)
        free(is_event);

    return err;
}

/* Price in the simple lognormal model (1 / FX(T) - K)+ * 1{Hd < FX(T) < Hu} */
Err ConvexFwdSimple(
    /*	Market parameters */
    long   today,
    char*  dom_yc,
    char*  for_yc,
    double spot_fx, /* as quoted in the market */

    /*	Product parameters */
    double notional,
    long   exercise_date,
    long   settlmt_date,
    double strike,
    int    is_call,
    double barrier_up,
    double barrier_down,

    /*	Calib parameters */
    long*   opt_matdates, /* mkt options maturity */
    double* opt_atmvols,  /* mkt options ATM vols */
    int     nb_opt,

    /* IOD / EOD flags */
    int    eod_fix_flag,   /*	EOD Fixing Flag 0: I, 1: E */
    int    eod_pay_flag,   /*	EOD Payment Flag 0: I, 1: E */
    int    eod_ex_flag,    /*	EOD Exercise Flag 0: I, 1: E */
    double spot_fx_fixing, /*	fixing of the spot fx if today > exercise */
    int    exercised,

    double* price,
    double* vol)
{
    double du, dd;
    double cash_fx, fwd_fx, mat;
    double cvx_vol, cvx_std;
    double df;
    long   spot_date, settlmt_fx_date;
    int    i;
    int    remove_option = 0;
    Err    err           = NULL;

    /* First check that the inputs make sense */
    if (exercise_date > settlmt_date)
    {
        err = "settlmt date should be greater than exercise date";
        goto FREE_RETURN;
    }

    if (barrier_up < -1.0E-10)
    {
        err = "barrier up should be positive";
        goto FREE_RETURN;
    }

    if (barrier_up < 1.0E-10)
    {
        barrier_up = spot_fx * 1.0E20;
    }

    if (barrier_up < barrier_down)
    {
        err = "barrier up should be greater than barrier down";
        goto FREE_RETURN;
    }

    if (strike < -1.0E-10)
    {
        err = "strike should be positive";
        goto FREE_RETURN;
    }

    /* check all the flags */
    if (today + eod_pay_flag > settlmt_date)
    {
        *price = 0.0;
        *vol   = 0.0;
        goto FREE_RETURN;
    }

    df = swp_f_df(today, settlmt_date, dom_yc);

    if (strike < 1.0E-10)
    {
        /* no option */
        if (today + eod_fix_flag > exercise_date)
        {
            if (spot_fx_fixing <= barrier_up && spot_fx_fixing >= barrier_down)
            {
                *price = (1.0 / spot_fx_fixing) * df * notional;
                *vol   = 0.0;
                goto FREE_RETURN;
            }
            else
            {
                *price = 0.0;
                *vol   = 0.0;
                goto FREE_RETURN;
            }
        }
    }
    else
    {
        if (today + eod_fix_flag > exercise_date && today + eod_ex_flag > exercise_date)
        {
            /* Fx has been fixed and decision has been taken but deal is not paid yet */
            if (!exercised || spot_fx_fixing > barrier_up || spot_fx_fixing < barrier_down)
            {
                *price = 0.0;
                *vol   = 0.0;
                goto FREE_RETURN;
            }
            else
            {
                if (is_call)
                {
                    *price = (1.0 / spot_fx_fixing - strike) * df * notional;
                    *vol   = 0.0;
                    goto FREE_RETURN;
                }
                else
                {
                    *price = (strike - 1.0 / spot_fx_fixing) * df * notional;
                    *vol   = 0.0;
                    goto FREE_RETURN;
                }
            }
        }
        else if (today + eod_fix_flag > exercise_date)
        {
            /* Fixed but not exercised */
            if (spot_fx_fixing <= barrier_up && spot_fx_fixing >= barrier_down)
            {
                if (is_call)
                {
                    *price = max(1.0 / spot_fx_fixing - strike, 0.0) * df * notional;
                    *vol   = 0.0;
                    goto FREE_RETURN;
                }
                else
                {
                    *price = max(strike - 1.0 / spot_fx_fixing, 0.0) * df * notional;
                    *vol   = 0.0;
                    goto FREE_RETURN;
                }
            }
            else
            {
                *price = 0.0;
                *vol   = 0.0;
                goto FREE_RETURN;
            }
        }
        else if (today + eod_ex_flag > exercise_date)
        {
            /* Exercised but not yet fixed */
            if (!exercised)
            {
                *price = 0.0;
                *vol   = 0.0;
                goto FREE_RETURN;
            }
            else
            {
                /* exercised but fx not yet fixed */
                remove_option = 1;
            }
        }
    }

    spot_date       = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
    settlmt_fx_date = add_unit(exercise_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
    cash_fx = spot_fx * swp_f_df(today, spot_date, dom_yc) / swp_f_df(today, spot_date, for_yc);
    fwd_fx =
        cash_fx * swp_f_df(today, settlmt_date, for_yc) / swp_f_df(today, settlmt_date, dom_yc);
    mat = (exercise_date - today) * YEARS_IN_DAY;

    if (mat < 1.0E-08)
    {
        mat = YEARS_IN_DAY / 10.0;
    }

    /* get the volatility */

    if (nb_opt > 1 && settlmt_fx_date > opt_matdates[0])
    {
        i = 0;
        while (i < nb_opt && opt_matdates[i] < settlmt_fx_date)
        {
            i++;
        }
        if (i < nb_opt)
        {
            i++;
        }

        if (opt_matdates[i - 1] > settlmt_fx_date)
        {
            /* linear interpolation */
            cvx_vol = (opt_atmvols[i - 1] - opt_atmvols[i - 2]) /
                          (opt_matdates[i - 1] - opt_matdates[i - 2]) *
                          (settlmt_fx_date - opt_matdates[i - 2]) +
                      opt_atmvols[i - 2];
        }
        else
        {
            cvx_vol = opt_atmvols[i - 1];
        }
    }
    else
    {
        cvx_vol = opt_atmvols[0];
    }

    *vol = cvx_vol;

    cvx_std = cvx_vol * sqrt(mat);

    if (fabs(strike) < 1.0E-10)
    {
        /* this is a simple convex forward */
        du = log(fwd_fx / barrier_up) / cvx_std + 0.5 * cvx_std;

        if (fabs(barrier_down) > 1.0E-08)
        {
            dd = log(fwd_fx / barrier_down) / cvx_std + 0.5 * cvx_std;
        }
        else
        {
            dd = 1.0E20;
        }

        *price =
            exp(cvx_std * cvx_std) / fwd_fx * (norm(dd - 2.0 * cvx_std) - norm(du - 2.0 * cvx_std));
    }
    else
    {
        if (is_call)
        {
            if (!remove_option)
            {
                barrier_up = min(barrier_up, 1.0 / strike);
            }
            du = log(fwd_fx / barrier_up) / cvx_std + 0.5 * cvx_std;

            if (fabs(barrier_down) > 1.0E-08)
            {
                dd = log(fwd_fx / barrier_down) / cvx_std + 0.5 * cvx_std;
            }
            else
            {
                dd = 1.0E20;
            }

            *price = exp(cvx_std * cvx_std) / fwd_fx *
                     (norm(dd - 2.0 * cvx_std) - norm(du - 2.0 * cvx_std));

            *price -= strike * (norm(dd - cvx_std) - norm(du - cvx_std));
        }
        else
        {
            if (!remove_option)
            {
                barrier_down = max(barrier_down, 1.0 / strike);
            }
            du = log(fwd_fx / barrier_up) / cvx_std + 0.5 * cvx_std;

            if (fabs(barrier_down) > 1.0E-08)
            {
                dd = log(fwd_fx / barrier_down) / cvx_std + 0.5 * cvx_std;
            }
            else
            {
                dd = 1.0E20;
            }

            *price = strike * (norm(dd - cvx_std) - norm(du - cvx_std));
            *price -= exp(cvx_std * cvx_std) / fwd_fx *
                      (norm(dd - 2.0 * cvx_std) - norm(du - 2.0 * cvx_std));
        }
    }

    *price *= df * notional;

FREE_RETURN:

    return err;
}

Err cvx_static_replication(double* strikes, long nb_strike, double* coefJ, double* res)
{
    int i;
    Err err = NULL;

    if (nb_strike == 1)
    {
        res[0] = 4.0 / strikes[nb_strike - 1] / strikes[nb_strike - 1];
    }
    else
    {
        res[nb_strike - 1] = 0.0;
        res[nb_strike - 2] = 4.0 * (strikes[nb_strike - 1] - strikes[nb_strike - 2]) /
                             strikes[nb_strike - 1] / strikes[nb_strike - 1];

        for (i = nb_strike - 3; i >= 0; i--)
        {
            res[i] = pow(res[i + 1] / (1.0 - sqrt(1 - res[i + 1] * strikes[i + 1])), 2) *
                         (strikes[i + 1] - strikes[i]) +
                     res[i + 1];
        }

        for (i = nb_strike - 1; i >= 1; i--)
        {
            coefJ[i] = (res[i] - res[i - 1]) / (strikes[i] - strikes[i - 1]);
        }

        coefJ[0] = -pow(res[0] / (1.0 - sqrt(1.0 - res[0] * strikes[0])), 2);

        res[nb_strike - 1] = -coefJ[nb_strike - 1];

        for (i = nb_strike - 2; i >= 0; i--)
        {
            res[i] = -coefJ[i] + coefJ[i + 1];
        }
    }

    return err;
}

static LOGCONVEX_PARAM logcvx_param;

void init_static_logconvex(
    double          mat,
    double          df,
    double          fwd,
    double          barrier_down,
    double          barrier_up,
    double          log_vol,
    long            nb_strike,
    double          dig_spread,
    LOGCONVEX_PARAM cvx_param)
{
    cvx_param->mat          = mat;
    cvx_param->df           = df;
    cvx_param->fwd          = fwd;
    cvx_param->barrier_down = barrier_down;
    cvx_param->barrier_up   = barrier_up;
    cvx_param->log_vol      = log_vol;
    cvx_param->log_std      = log_vol * sqrt(mat);
    cvx_param->nb_strike    = nb_strike;

    if (dig_spread < 1.0E-10)
    {
        cvx_param->dig_spread = DIGIT_SPREAD * fwd;
    }
    else
    {
        cvx_param->dig_spread = dig_spread;
    }

    cvx_param->coef    = dvector(0, nb_strike - 1);
    cvx_param->coefJ   = dvector(0, nb_strike - 1);
    cvx_param->strikes = dvector(0, nb_strike - 1);
}

void free_static_logconvex(LOGCONVEX_PARAM cvx_param)
{
    if (cvx_param)
    {
        if (cvx_param->coef)
            free_dvector(cvx_param->coef, 0, cvx_param->nb_strike - 1);
        if (cvx_param->coefJ)
            free_dvector(cvx_param->coefJ, 0, cvx_param->nb_strike - 1);
        if (cvx_param->strikes)
            free_dvector(cvx_param->strikes, 0, cvx_param->nb_strike - 1);

        free(cvx_param);
        cvx_param = NULL;
    }
}

Err static_logconvex(double param, double std[], double* price, double* gradient, int nb_strike)
{
    Err    err = NULL;
    int    i, j, index_up, index_down;
    double sum;
    double slope_down, const_down, strike_down, strike_up, slope_up, const_up;
    double strike_down_down, strike_down_up, coef_down, digit_down;
    double strike_up_down, strike_up_up, coef_up, digit_up;

    /* check that strikes are increasing */
    for (i = 2; i <= nb_strike; i++)
    {
        if (std[i] < 0.0)
        {
            std[i] = 2.0 * GRADIENT_SHIFT;
        }
    }

    /* check input are not too bad */
    if (std[1] > 30)
    {
        std[1] = -10;
    }

    for (i = 2; i <= nb_strike; i++)
    {
        if (std[i] > 30)
        {
            std[i] = 5;
        }
    }

    /* first convert std in strikes */
    sum = 0.0;
    for (i = 0; i < nb_strike; i++)
    {
        sum += std[i + 1];
        logcvx_param->strikes[i] = logcvx_param->fwd * exp(sum * logcvx_param->log_std);
    }

    sum = 0.0;

    /* Static replication */
    /* take into account the up barrier */
    if (logcvx_param->barrier_up > 1.0E-10)
    {
        index_up =
            Get_Index(logcvx_param->barrier_up, logcvx_param->strikes, logcvx_param->nb_strike);
    }
    else
    {
        index_up = logcvx_param->nb_strike - 1;
    }

    err = cvx_static_replication(
        logcvx_param->strikes, index_up + 1, logcvx_param->coefJ, logcvx_param->coef);

    if (err)
    {
        goto FREE_RETURN;
    }

    /* take into account the down barrier */
    index_down = Get_Index(logcvx_param->barrier_down, logcvx_param->strikes, index_up + 1);

    for (i = 0; i < index_down; i++)
    {
        logcvx_param->coef[i] = 0.0;
    }

    if (logcvx_param->barrier_down > 1.0E-10)
    {
        slope_down = 0.0;
        const_down = 0.0;

        for (i = index_down; i < index_up + 1; i++)
        {
            slope_down += logcvx_param->coef[i];
            const_down += logcvx_param->coef[i] * logcvx_param->strikes[i];
        }
    }

    if (logcvx_param->barrier_up > 1.0E-10)
    {
        /* remove barrier up */
        slope_up = logcvx_param->coef[index_up];
        const_up = logcvx_param->coef[index_up] * logcvx_param->strikes[index_up];
        logcvx_param->coef[index_up] = 0.0;
    }

    /* Price computation */
    sum = 0.0;

    for (i = 0; i < index_up + 1; i++)
    {
        sum += logcvx_param->coef[i] * srt_f_optblksch(
                                           logcvx_param->fwd,
                                           logcvx_param->strikes[i],
                                           logcvx_param->log_std,
                                           1.0,
                                           logcvx_param->df,
                                           SRT_PUT,
                                           PREMIUM);
    }

    /* remove barrier up */
    if (logcvx_param->barrier_up > 1.0E-10)
    {
        strike_up = logcvx_param->barrier_up;
        digit_up  = const_up - slope_up * strike_up;

        /* add option part */

        sum += slope_up * srt_f_optblksch(
                              logcvx_param->fwd,
                              strike_up,
                              logcvx_param->log_std,
                              1.0,
                              logcvx_param->df,
                              SRT_PUT,
                              PREMIUM);

        /* add digital part */
        strike_up_up   = strike_up + logcvx_param->dig_spread;
        strike_up_down = strike_up - logcvx_param->dig_spread;
        coef_up        = digit_up / (2.0 * logcvx_param->dig_spread);

        sum += coef_up * srt_f_optblksch(
                             logcvx_param->fwd,
                             strike_up_up,
                             logcvx_param->log_std,
                             1.0,
                             logcvx_param->df,
                             SRT_PUT,
                             PREMIUM);

        sum -= coef_up * srt_f_optblksch(
                             logcvx_param->fwd,
                             strike_up_down,
                             logcvx_param->log_std,
                             1.0,
                             logcvx_param->df,
                             SRT_PUT,
                             PREMIUM);
    }

    /* remove barrier down */
    if (logcvx_param->barrier_down > 1.0E-10)
    {
        strike_down = logcvx_param->barrier_down;
        digit_down  = const_down - slope_down * strike_down;

        /* remove option part */
        sum -= slope_down * srt_f_optblksch(
                                logcvx_param->fwd,
                                strike_down,
                                logcvx_param->log_std,
                                1.0,
                                logcvx_param->df,
                                SRT_PUT,
                                PREMIUM);

        /* remove digital part */
        strike_down_up   = strike_down + logcvx_param->dig_spread;
        strike_down_down = strike_down - logcvx_param->dig_spread;
        coef_down        = digit_down / (2.0 * logcvx_param->dig_spread);

        sum -= coef_down * srt_f_optblksch(
                               logcvx_param->fwd,
                               strike_down_up,
                               logcvx_param->log_std,
                               1.0,
                               logcvx_param->df,
                               SRT_PUT,
                               PREMIUM);

        sum += coef_down * srt_f_optblksch(
                               logcvx_param->fwd,
                               strike_down_down,
                               logcvx_param->log_std,
                               1.0,
                               logcvx_param->df,
                               SRT_PUT,
                               PREMIUM);
    }

    *price = sum;

    /* calculates the gradient */

    for (j = 1; j <= nb_strike; j++)
    {
        /* shift the strike */
        std[j] += GRADIENT_SHIFT;
        sum = 0.0;
        for (i = 0; i < nb_strike; i++)
        {
            sum += std[i + 1];
            logcvx_param->strikes[i] = logcvx_param->fwd * exp(sum * logcvx_param->log_std);
        }

        sum = 0.0;

        /* Static replication */
        /* take into account the up barrier */
        if (logcvx_param->barrier_up > 1.0E-10)
        {
            index_up =
                Get_Index(logcvx_param->barrier_up, logcvx_param->strikes, logcvx_param->nb_strike);
        }
        else
        {
            index_up = logcvx_param->nb_strike - 1;
        }

        err = cvx_static_replication(
            logcvx_param->strikes, index_up + 1, logcvx_param->coefJ, logcvx_param->coef);

        if (err)
        {
            goto FREE_RETURN;
        }

        /* take into account the down barrier */
        index_down = Get_Index(logcvx_param->barrier_down, logcvx_param->strikes, index_up + 1);

        for (i = 0; i < index_down; i++)
        {
            logcvx_param->coef[i] = 0.0;
        }

        if (logcvx_param->barrier_down > 1.0E-10)
        {
            slope_down = 0.0;
            const_down = 0.0;

            for (i = index_down; i < index_up + 1; i++)
            {
                slope_down += logcvx_param->coef[i];
                const_down += logcvx_param->coef[i] * logcvx_param->strikes[i];
            }
        }

        if (logcvx_param->barrier_up > 1.0E-10)
        {
            /* remove barrier up */
            slope_up = logcvx_param->coef[index_up];
            const_up = logcvx_param->coef[index_up] * logcvx_param->strikes[index_up];
            logcvx_param->coef[index_up] = 0.0;
        }

        /* Price computation */
        sum = 0.0;

        for (i = 0; i < index_up + 1; i++)
        {
            sum += logcvx_param->coef[i] * srt_f_optblksch(
                                               logcvx_param->fwd,
                                               logcvx_param->strikes[i],
                                               logcvx_param->log_std,
                                               1.0,
                                               logcvx_param->df,
                                               SRT_PUT,
                                               PREMIUM);
        }

        /* remove barrier up */
        if (logcvx_param->barrier_up > 1.0E-10)
        {
            strike_up = logcvx_param->barrier_up;
            digit_up  = const_up - slope_up * strike_up;

            /* add option part */

            sum += slope_up * srt_f_optblksch(
                                  logcvx_param->fwd,
                                  strike_up,
                                  logcvx_param->log_std,
                                  1.0,
                                  logcvx_param->df,
                                  SRT_PUT,
                                  PREMIUM);

            /* add digital part */
            strike_up_up   = strike_up + logcvx_param->dig_spread;
            strike_up_down = strike_up - logcvx_param->dig_spread;
            coef_up        = digit_up / (2.0 * logcvx_param->dig_spread);

            sum += coef_up * srt_f_optblksch(
                                 logcvx_param->fwd,
                                 strike_up_up,
                                 logcvx_param->log_std,
                                 1.0,
                                 logcvx_param->df,
                                 SRT_PUT,
                                 PREMIUM);

            sum -= coef_up * srt_f_optblksch(
                                 logcvx_param->fwd,
                                 strike_up_down,
                                 logcvx_param->log_std,
                                 1.0,
                                 logcvx_param->df,
                                 SRT_PUT,
                                 PREMIUM);
        }

        /* remove barrier down */
        if (logcvx_param->barrier_down > 1.0E-10)
        {
            strike_down = logcvx_param->barrier_down;
            digit_down  = const_down - slope_down * strike_down;

            /* remove option part */
            sum -= slope_down * srt_f_optblksch(
                                    logcvx_param->fwd,
                                    strike_down,
                                    logcvx_param->log_std,
                                    1.0,
                                    logcvx_param->df,
                                    SRT_PUT,
                                    PREMIUM);

            /* remove digital part */
            strike_down_up   = strike_down + logcvx_param->dig_spread;
            strike_down_down = strike_down - logcvx_param->dig_spread;
            coef_down        = digit_down / (2.0 * logcvx_param->dig_spread);

            sum -= coef_down * srt_f_optblksch(
                                   logcvx_param->fwd,
                                   strike_down_up,
                                   logcvx_param->log_std,
                                   1.0,
                                   logcvx_param->df,
                                   SRT_PUT,
                                   PREMIUM);

            sum += coef_down * srt_f_optblksch(
                                   logcvx_param->fwd,
                                   strike_down_down,
                                   logcvx_param->log_std,
                                   1.0,
                                   logcvx_param->df,
                                   SRT_PUT,
                                   PREMIUM);
        }

        /* shift the strike */
        std[j] -= GRADIENT_SHIFT;

        gradient[j] = (sum - (*price)) / (GRADIENT_SHIFT);
    }

FREE_RETURN:

    return err;
}

Err ConvexFwdSABR_static_replication(
    /*	Market parameters */
    long   today,
    char*  dom_yc,
    char*  for_yc,
    double spot_fx, /* as quoted in the market */

    /*	Product parameters */
    double notional,
    long   exercise_date,
    long   settlmt_date,
    double barrier_up,
    double barrier_down,

    /*	Calib parameters */
    long*   opt_matdates, /* mkt options maturity */
    double* opt_atmvols,  /* mkt options ATM vols */
    int     nb_opt,

    double alpha,
    double beta,
    double rho,

    /*	Replication parameters */
    double* strikes,
    long    nb_strike,
    double  floorvol_strike,
    long    is_std,
    double  dig_spread,
    int     calib_strikes,

    /* IOD / EOD flags */
    int    eod_fix_flag,   /*	EOD Fixing Flag 0: I, 1: E */
    int    eod_pay_flag,   /*	EOD Payment Flag 0: I, 1: E */
    int    eod_ex_flag,    /*	EOD Exercise Flag 0: I, 1: E */
    double spot_fx_fixing, /*	fixing of the spot fx if today > exercise */
    int    exercised,

    double*  price,
    double*  coef,
    double*  vols,
    double** bar_adjust) /*	bar adjust is an array 5 * 3 */
{
    long    spot_date, settlmt_fx_date;
    double  cash_fx, fwd_fx, mat;
    double  atm_vol, atm_std, sigma_beta, floor_vol, sum, df;
    double  slope_down, const_down, strike_down, vol_down, strike_up, vol_up, slope_up, const_up;
    double  strike_down_up, vol_down_up, digit_down;
    double  strike_up_down, vol_up_down, digit_up;
    double* coefJ = NULL;
    int     i, index_up, index_down;
    double  log_price, log_vol;
    logconvex_param logconvex;
    double *        data = NULL, *weight = NULL, *target = NULL, *param = NULL, error;
    long            nb_iter = 20;

    Err err = NULL;

    if (barrier_up < barrier_down && barrier_up > 1.0E10)
    {
        err = "barrier up should be greater than barrier down";
        goto FREE_RETURN;
    }

    if (barrier_up < -1.0E-10)
    {
        err = "barrier up should be positive";
        goto FREE_RETURN;
    }

    /* check all the flags */
    if (today + eod_pay_flag > settlmt_date)
    {
        *price = 0.0;
        goto FREE_RETURN;
    }

    df = swp_f_df(today, settlmt_date, dom_yc);

    /* no option */
    if (today + eod_fix_flag > exercise_date)
    {
        if (spot_fx_fixing <= barrier_up && spot_fx_fixing >= barrier_down)
        {
            *price = (1.0 / spot_fx_fixing) * df * notional;
            goto FREE_RETURN;
        }
        else
        {
            *price = 0.0;
            goto FREE_RETURN;
        }
    }

    spot_date       = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
    settlmt_fx_date = add_unit(exercise_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
    cash_fx = spot_fx * swp_f_df(today, spot_date, dom_yc) / swp_f_df(today, spot_date, for_yc);
    fwd_fx =
        cash_fx * swp_f_df(today, settlmt_date, for_yc) / swp_f_df(today, settlmt_date, dom_yc);
    mat = (exercise_date - today) * YEARS_IN_DAY;

    /* get the volatility */
    if (nb_opt > 1 && settlmt_fx_date > opt_matdates[0])
    {
        i = 0;
        while (i < nb_opt && opt_matdates[i] < settlmt_fx_date)
        {
            i++;
        }
        if (i < nb_opt)
        {
            i++;
        }

        if (opt_matdates[i - 1] > settlmt_fx_date)
        {
            /* linear interpolation */
            atm_vol = (opt_atmvols[i - 1] - opt_atmvols[i - 2]) /
                          (opt_matdates[i - 1] - opt_matdates[i - 2]) *
                          (settlmt_fx_date - opt_matdates[i - 2]) +
                      opt_atmvols[i - 2];
        }
        else
        {
            atm_vol = opt_atmvols[i - 1];
        }
    }
    else
    {
        atm_vol = opt_atmvols[0];
    }

    if (calib_strikes)
    {
        data   = dvector(0, 1);
        weight = dvector(0, 1);
        target = dvector(0, 1);
        param  = dvector(0, nb_strike);

        if (!data || !weight || !target || !param)
        {
            err = "Memory allocation faillure in ConvexFwdSABR_static_replication";
            goto FREE_RETURN;
        }

        /* first we try to find the strikes which match as much as possible the log price */

        /* First guess */

        sum = 0.0;
        for (i = 0; i < nb_strike; i++)
        {
            sum += strikes[i] * strikes[i];
        }

        if (sum < 1.0E-10 && nb_strike > 1)
        {
            /* we provide the first guess */
            sum      = 10.0 / nb_strike;
            param[1] = -5.0;
            for (i = 2; i <= nb_strike; i++)
            {
                param[i] = sum;
            }
        }
        else if (!is_std)
        {
            /* convert into std */
            param[1] = log(fwd_fx / strikes[0]) / atm_vol / sqrt(mat);

            for (i = 2; i <= nb_strike; i++)
            {
                param[i] = log(fwd_fx / strikes[i - 1]) / atm_vol / sqrt(mat) - strikes[i - 2];
            }
        }
        else
        {
            for (i = nb_strike; i >= 2; i--)
            {
                param[i] = strikes[i - 1] - strikes[i - 2];
            }

            param[1] = strikes[0];
        }

        logcvx_param = &logconvex;

        err = ConvexFwdSimple(
            today,
            dom_yc,
            for_yc,
            spot_fx,
            1.0,
            exercise_date,
            settlmt_date,
            0.0,
            0,
            barrier_up,
            barrier_down,
            opt_matdates,
            opt_atmvols,
            nb_opt,
            0,
            0,
            0,
            0.0,
            0,
            &log_price,
            &log_vol);

        if (err)
            goto FREE_RETURN;

        init_static_logconvex(
            mat,
            df,
            fwd_fx,
            barrier_down,
            barrier_up,
            log_vol,
            nb_strike,
            dig_spread,
            logcvx_param);

        /* optimise */
        data[1]   = 0.0;
        weight[1] = 1.0;
        target[1] = log_price;

        err = levenberg_marquardt(
            data, target, weight, 1, param, nb_strike, nb_iter, static_logconvex, &error);

        /* convert into strikes */
        sum = 0.0;

        for (i = 0; i < nb_strike; i++)
        {
            sum += param[i + 1];
            strikes[i] = sum;
        }

        is_std = 1;
    }
    else
    {
        logcvx_param = NULL;
    }

    /* transform the strikes */
    if (is_std)
    {
        atm_std = atm_vol * sqrt(mat);

        for (i = 0; i < nb_strike; i++)
        {
            strikes[i] = fwd_fx * exp(strikes[i] * atm_std);
        }

        floorvol_strike = fwd_fx * exp(floorvol_strike * atm_std);
    }

    /* Static replication */

    /* take into account the up barrier */
    if (barrier_up > 1.0E-10)
    {
        index_up = Get_Index(barrier_up, strikes, nb_strike);

        if (strikes[index_up] - barrier_up > 1.0E-08)
        {
            index_up--;
        }
    }
    else
    {
        index_up = nb_strike - 1;
    }

    coefJ = dvector(0, nb_strike - 1);

    if (!coefJ)
    {
        err = "Memory allocation faillure in ConvexFwdSABR_static_replication";
        goto FREE_RETURN;
    }

    /* first replicate the convex forward */
    err = cvx_static_replication(strikes, nb_strike, coefJ, coef);

    if (err)
        goto FREE_RETURN;

    /* take into account the up barrier */
    slope_up = 0.0;
    const_up = 0.0;

    if (barrier_up > 1.0E-10)
    {
        for (i = index_up + 1; i < nb_strike; i++)
        {
            slope_up += coef[i];
            const_up += coef[i] * (strikes[i] - strikes[index_up]);
            coef[i] = 0.0;
        }

        coef[index_up] += slope_up;
    }

    /* take into account the down barrier */
    slope_down = 0.0;
    const_down = const_up;

    if (barrier_down > 1.0E-10)
    {
        index_down = Get_Index(barrier_down, strikes, index_up + 1);

        if (barrier_down - strikes[index_down] > 1.0E-08)
        {
            index_down++;
        }

        for (i = 0; i < index_down; i++)
        {
            coef[i] = 0.0;
        }

        for (i = index_down; i < index_up + 1; i++)
        {
            slope_down += coef[i];
            const_down += coef[i] * (strikes[i] - strikes[index_down]);
        }
    }

    /* Price computation */
    /* ***************** */

    /* Get the sigma beta */
    err = srt_f_optsarbvol(
        fwd_fx, fwd_fx, mat, atm_vol, alpha, beta, rho, SRT_LOGNORMAL, SRT_BETAVOL, &sigma_beta);

    if (err)
        goto FREE_RETURN;

    /* Calculate vol at floor */
    err = srt_f_optsarbvol(
        fwd_fx,
        floorvol_strike,
        mat,
        sigma_beta,
        alpha,
        beta,
        rho,
        SRT_BETAVOL,
        SRT_LOGNORMAL,
        &floor_vol);

    if (err)
        goto FREE_RETURN;

    sum = 0.0;
    df  = swp_f_df(today, settlmt_date, dom_yc);

    for (i = 0; i < nb_strike; i++)
    {
        if (fabs(coef[i]) > 1.0E-10)
        {
            if (strikes[i] > floorvol_strike)
            {
                err = srt_f_optsarbvol(
                    fwd_fx,
                    strikes[i],
                    mat,
                    sigma_beta,
                    alpha,
                    beta,
                    rho,
                    SRT_BETAVOL,
                    SRT_LOGNORMAL,
                    &(vols[i]));

                if (err)
                    goto FREE_RETURN;
            }
            else
            {
                vols[i] = floor_vol;
            }

            sum +=
                coef[i] * srt_f_optblksch(fwd_fx, strikes[i], vols[i], mat, df, SRT_PUT, PREMIUM);
        }
    }

    if (dig_spread < 1.0E-10)
    {
        dig_spread = DIGIT_SPREAD * fwd_fx;
    }

    /* Remove upper part of the Convex */
    if (barrier_up > 1.0E-10)
    {
        strike_up = strikes[index_up];
        digit_up  = const_up / dig_spread;

        /* substract digital */
        if (strike_up > floorvol_strike)
        {
            err = srt_f_optsarbvol(
                fwd_fx,
                strike_up,
                mat,
                sigma_beta,
                alpha,
                beta,
                rho,
                SRT_BETAVOL,
                SRT_LOGNORMAL,
                &vol_up);
        }
        else
        {
            vol_up = floor_vol;
        }

        sum += digit_up * srt_f_optblksch(fwd_fx, strike_up, vol_up, mat, df, SRT_PUT, PREMIUM);

        strike_up_down = strike_up - dig_spread;

        if (strike_up_down > floorvol_strike)
        {
            err = srt_f_optsarbvol(
                fwd_fx,
                strike_up_down,
                mat,
                sigma_beta,
                alpha,
                beta,
                rho,
                SRT_BETAVOL,
                SRT_LOGNORMAL,
                &vol_up_down);
        }
        else
        {
            vol_up_down = floor_vol;
        }

        sum -= digit_up *
               srt_f_optblksch(fwd_fx, strike_up_down, vol_up_down, mat, df, SRT_PUT, PREMIUM);

        bar_adjust[2][0] = digit_up;
        bar_adjust[2][1] = vol_up;
        bar_adjust[2][2] = strike_up;
        bar_adjust[3][0] = -digit_up;
        bar_adjust[3][1] = vol_up_down;
        bar_adjust[3][2] = strike_up_down;
    }

    /* Remove lower part of the digital */
    if (barrier_down > 1.0E-10)
    {
        strike_down = strikes[index_down];
        digit_down  = const_down / dig_spread;

        /* remove option part */
        if (strike_down > floorvol_strike)
        {
            err = srt_f_optsarbvol(
                fwd_fx,
                strike_down,
                mat,
                sigma_beta,
                alpha,
                beta,
                rho,
                SRT_BETAVOL,
                SRT_LOGNORMAL,
                &vol_down);
        }
        else
        {
            vol_down = floor_vol;
        }

        sum -=
            slope_down * srt_f_optblksch(fwd_fx, strike_down, vol_down, mat, df, SRT_PUT, PREMIUM);

        /* remove digital part */
        strike_down_up = strike_down + dig_spread;

        sum +=
            digit_down * srt_f_optblksch(fwd_fx, strike_down, vol_down, mat, df, SRT_PUT, PREMIUM);

        if (strike_down_up > floorvol_strike)
        {
            err = srt_f_optsarbvol(
                fwd_fx,
                strike_down_up,
                mat,
                sigma_beta,
                alpha,
                beta,
                rho,
                SRT_BETAVOL,
                SRT_LOGNORMAL,
                &vol_down_up);
        }
        else
        {
            vol_down_up = floor_vol;
        }

        sum -= digit_down *
               srt_f_optblksch(fwd_fx, strike_down_up, vol_down_up, mat, df, SRT_PUT, PREMIUM);

        bar_adjust[0][0] = -slope_down + digit_down;
        bar_adjust[0][1] = vol_down;
        bar_adjust[0][2] = strike_down;
        bar_adjust[1][0] = -digit_down;
        bar_adjust[1][1] = vol_down_up;
        bar_adjust[1][2] = strike_down_up;
    }

    *price = sum * notional;

FREE_RETURN:

    if (coefJ)
        free_dvector(coefJ, 0, nb_strike - 1);
    if (data)
        free_dvector(data, 0, 1);
    if (weight)
        free_dvector(weight, 0, 1);
    if (target)
        free_dvector(target, 0, 1);
    if (param)
        free_dvector(param, 0, nb_strike);

    return err;
}

Err fwd_fxoption(
    long    today,
    long    spot_date,
    double  spot_fx, /*	2bd fwd */
    char*   dom_yc,
    char*   dom_vc,
    char*   dom_ref,
    char*   dom_swap_freq,
    char*   dom_swap_basis,
    double  dom_lam,
    char*   for_yc,
    char*   for_vc,
    char*   for_ref,
    char*   for_swap_freq,
    char*   for_swap_basis,
    double  for_lam,
    double* corr_times,
    double* correl_dom_for, /*	Correlations */
    double* correl_dom_fx,
    double* correl_for_fx,
    long    corr_n_times,
    Err (*get_cash_vol)(/*	Function to get IR cash vol from the markets */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),

    /*	Fx vol from the market */
    long*   fx_mkt_vol_date,
    double* fx_mkt_vol,
    int     num_fx_mkt_vol,

    /*	The structure */
    double notional,
    long   exercise_date,
    long   fx_fixing_date,
    long   settlmt_date,
    double strike,
    int    callput,

    /* IOD / EOD flags */
    int eod_fix_flag, /*	EOD Fixing Flag 0: I, 1: E */
    int eod_pay_flag, /*	EOD Payment Flag 0: I, 1: E */
    int eod_ex_flag,  /*	EOD Exercise Flag 0: I, 1: E */

    int    exercised,      /*	Is exercised Flag */
    double spot_fx_fixing, /*	Fixing of the Fx */

    /* Extra params */
    char*  calib_freq,
    long   calib_dom_tau,
    long   calib_for_tau,
    double dom_vol_shift,
    double for_vol_shift,
    double fx_vol_shift,

    /* Outputs */
    double* price,
    double* vol,
    int     export_ts, /*	1: Export TS, 0: don't */
    CPD_UND und_exp)
{
    double *fx_vol_curve = NULL, *sigma_date_dom = NULL, *sigma_dom = NULL, *sigma_date_for = NULL,
           *sigma_for = NULL, *merge_dates = NULL, *sig_dom = NULL, *sig_for = NULL,
           *opt_exe = NULL, *opt_mat = NULL;

    long* ex_date = NULL;

    long           sigma_n_dom, sigma_n_for, nb_merge_dates, num_calib, num_calib_rates;
    SrtCompounding call_comp, dom_comp, for_comp;
    int            i, call_nb, nb_exe, call_nb_max;
    long           cur_date, cur_date_adj;
    double         opt_maturity, exe_maturity;
    long           opt_date_lim, maturity_calib;
    double         df, cash_fx, fwd_fx;
    SrtCallPutType CallPutFlag;

    Err err = NULL;

    /* check all the flags */
    if (today + eod_pay_flag > settlmt_date)
    {
        *price = 0.0;
        *vol   = 0.0;
        goto FREE_RETURN;
    }

    df      = swp_f_df(today, settlmt_date, dom_yc);
    cash_fx = spot_fx * swp_f_df(today, spot_date, dom_yc) / swp_f_df(today, spot_date, for_yc);
    fwd_fx =
        cash_fx * swp_f_df(today, settlmt_date, for_yc) / swp_f_df(today, settlmt_date, dom_yc);

    if (today + eod_ex_flag > exercise_date && today + eod_fix_flag > fx_fixing_date)
    {
        /* Fx has been fixed and decision has been taken but deal is not paid yet */
        if (!exercised)
        {
            *price = 0.0;
            *vol   = 0.0;
            goto FREE_RETURN;
        }
        else
        {
            if (callput)
            {
                *price = (spot_fx_fixing - strike) * df;
                *vol   = 0.0;
                goto FREE_RETURN;
            }
            else
            {
                *price = (strike - spot_fx_fixing) * df;
                *vol   = 0.0;
                goto FREE_RETURN;
            }
        }
    }
    else if (today + eod_fix_flag > fx_fixing_date)
    {
        /* Fixed but not exercised */

        if (callput)
        {
            *price = max(spot_fx_fixing - strike, 0.0) * df;
            *vol   = 0.0;
            goto FREE_RETURN;
        }
        else
        {
            *price = max(strike - spot_fx_fixing, 0.0) * df;
            *vol   = 0.0;
            goto FREE_RETURN;
        }
    }
    else if (today + eod_ex_flag > exercise_date)
    {
        /* Exercised but not yet fixed */
        if (!exercised)
        {
            *price = 0.0;
            *vol   = 0.0;
            goto FREE_RETURN;
        }
        else
        {
            /* exercised but fx not yet fixed */
            if (callput)
            {
                *price = (fwd_fx - strike) * df;
                *vol   = 0.0;
                goto FREE_RETURN;
            }
            else
            {
                *price = (strike - fwd_fx) * df;
                *vol   = 0.0;
                goto FREE_RETURN;
            }
        }
    }

    /* first calibrate the underlying */
    opt_maturity = (settlmt_date - today) * YEARS_IN_DAY;
    exe_maturity = (exercise_date - today) * YEARS_IN_DAY;

    err = interp_compounding(calib_freq, &call_comp);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = interp_compounding(dom_swap_freq, &dom_comp);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = interp_compounding(for_swap_freq, &for_comp);
    if (err)
    {
        goto FREE_RETURN;
    }

    call_nb_max = 12 / (min(min(call_comp, dom_comp), for_comp));

    opt_date_lim = add_unit(spot_date, 2 * call_nb_max, SRT_MONTH, MODIFIED_SUCCEEDING);

    if (settlmt_date < opt_date_lim)
    {
        maturity_calib = opt_date_lim;
    }
    else
    {
        maturity_calib = settlmt_date;
    }

    /* calibration dates */

    call_nb = 12 / call_comp;

    cur_date     = maturity_calib;
    cur_date_adj = add_unit(cur_date, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

    nb_exe = 0;

    while (cur_date_adj > spot_date)
    {
        cur_date     = add_unit(cur_date, -call_nb, SRT_MONTH, NO_BUSDAY_CONVENTION);
        cur_date_adj = add_unit(cur_date, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

        nb_exe++;
    }

    nb_exe--;

    ex_date = calloc(nb_exe, sizeof(long));

    if (!ex_date)
    {
        err = "Memory Allocation faillure in FX3DVolDef2";
        goto FREE_RETURN;
    }

    cur_date = maturity_calib;

    for (i = nb_exe - 1; i >= 0; i--)
    {
        cur_date     = add_unit(cur_date, -call_nb, SRT_MONTH, NO_BUSDAY_CONVENTION);
        cur_date_adj = add_unit(cur_date, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

        ex_date[i] = add_unit(cur_date_adj, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
    }

    for (num_calib_rates = 0;
         (num_calib_rates < nb_exe) && (ex_date[num_calib_rates] < exercise_date);
         num_calib_rates++)
        ;

    if (num_calib_rates >= nb_exe)
        num_calib_rates--;

    num_calib_rates = max(num_calib_rates, 1);

    err = cpd_calib_diagonal(
        dom_yc,
        dom_vc,
        dom_ref,
        get_cash_vol,
        dom_vol_shift,
        1,
        num_calib_rates,
        ex_date,
        maturity_calib,
        NULL,
        NULL,
        0,
        1.0,
        1.0,
        dom_swap_freq,
        dom_swap_basis,
        !calib_dom_tau,
        0,
        0,
        CALPRES,
        CALPRES,
        0,
        0,
        0,
        0,
        NULL,
        &dom_lam,
        1,
        0.0,
        0.0,
        0.0,
        &sigma_n_dom,
        &sigma_date_dom,
        &sigma_dom,
        NULL);

    if (err)
    {
        goto FREE_RETURN;
    }

    err = cpd_calib_diagonal(
        for_yc,
        for_vc,
        for_ref,
        get_cash_vol,
        for_vol_shift,
        1,
        num_calib_rates,
        ex_date,
        maturity_calib,
        NULL,
        NULL,
        0,
        1.0,
        1.0,
        for_swap_freq,
        for_swap_basis,
        !calib_for_tau,
        0,
        0,
        CALPRES,
        CALPRES,
        0,
        0,
        0,
        0,
        NULL,
        &for_lam,
        1,
        0.0,
        0.0,
        0.0,
        &sigma_n_for,
        &sigma_date_for,
        &sigma_for,
        NULL);

    if (err)
    {
        goto FREE_RETURN;
    }

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

    /* transform dates into maturities */

    opt_exe = dvector(0, num_fx_mkt_vol - 1);
    opt_mat = dvector(0, num_fx_mkt_vol - 1);

    if (!opt_exe || !opt_mat)
    {
        err = "Memory allocation faillure in fwdfx_option (1)";
        goto FREE_RETURN;
    }

    for (i = 0; i < num_fx_mkt_vol; i++)
    {
        opt_exe[i] = add_unit(fx_mkt_vol_date[i], -2, SRT_BDAY, MODIFIED_SUCCEEDING);
        opt_mat[i] = (fx_mkt_vol_date[i] - today) * YEARS_IN_DAY;
        opt_exe[i] = (opt_exe[i] - today) * YEARS_IN_DAY;
    }

    num_calib = Get_Index(exe_maturity, opt_exe, num_fx_mkt_vol);
    num_calib = max(min(num_calib + 1, num_fx_mkt_vol), 1);

    err = Fx3DtsCalibration_corr(
        opt_exe,
        opt_mat,
        fx_mkt_vol,
        num_calib,
        merge_dates,
        nb_merge_dates,
        sig_dom,
        dom_lam,
        sig_for,
        for_lam,
        corr_times,
        correl_dom_for,
        correl_dom_fx,
        correl_for_fx,
        corr_n_times,
        &fx_vol_curve);

    if (err)
    {
        goto FREE_RETURN;
    }

    /* calculates the volatility */

    err = Fx3DtsImpliedVol_corr(
        opt_maturity,
        0.0,
        exe_maturity,
        merge_dates,
        nb_merge_dates,
        sig_dom,
        dom_lam,
        sig_for,
        for_lam,
        opt_exe,
        fx_vol_curve,
        num_calib,
        corr_times,
        correl_dom_for,
        correl_dom_fx,
        correl_for_fx,
        corr_n_times,
        vol);

    if (err)
    {
        goto FREE_RETURN;
    }

    if (callput)
    {
        CallPutFlag = SRT_CALL;
    }
    else
    {
        CallPutFlag = SRT_PUT;
    }

    *price =
        notional * srt_f_optblksch(fwd_fx, strike, *vol, exe_maturity, df, CallPutFlag, PREMIUM);

FREE_RETURN:

    if (export_ts)
    {
        if (merge_dates)
        {
            for (i = 0; i < nb_merge_dates; i++)
            {
                merge_dates[i] = today + merge_dates[i] * DAYS_IN_YEAR;
            }

            und_exp->sigma_date_rates = merge_dates;
            und_exp->sigma_n_rates    = nb_merge_dates;
        }
        else
        {
            und_exp->sigma_n_rates    = 0;
            und_exp->sigma_date_rates = NULL;
        }

        if (sig_dom)
        {
            und_exp->sigma_dom = sig_dom;
            und_exp->lda_dom   = dom_lam;
        }
        else
        {
            und_exp->sigma_dom = NULL;
            und_exp->lda_dom   = dom_lam;
        }

        if (sig_for)
        {
            und_exp->sigma_for = sig_for;
            und_exp->lda_for   = for_lam;
        }
        else
        {
            und_exp->sigma_for = NULL;
            und_exp->lda_for   = for_lam;
        }

        if (fx_vol_curve)
        {
            und_exp->sigma_date_fx = calloc(num_calib, sizeof(double));

            for (i = 0; i < num_calib; i++)
            {
                und_exp->sigma_date_fx[i] =
                    add_unit(fx_mkt_vol_date[i], -2, SRT_BDAY, MODIFIED_SUCCEEDING);
            }
            und_exp->sigma_fx   = fx_vol_curve;
            und_exp->sigma_n_fx = num_calib;
        }
        else
        {
            und_exp->sigma_date_fx = NULL;
            und_exp->sigma_fx      = NULL;
            und_exp->sigma_n_fx    = 0;
        }
    }
    else
    {
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

        if (fx_vol_curve)
        {
            free(fx_vol_curve);
        }
    }

    if (ex_date)
    {
        free(ex_date);
    }

    if (opt_exe)
    {
        free_dvector(opt_exe, 0, num_fx_mkt_vol - 1);
    }

    if (opt_mat)
    {
        free_dvector(opt_mat, 0, num_fx_mkt_vol - 1);
    }

    if (sigma_date_dom)
    {
        free(sigma_date_dom);
    }

    if (sigma_date_for)
    {
        free(sigma_date_for);
    }

    if (sigma_dom)
    {
        free(sigma_dom);
    }

    if (sigma_for)
    {
        free(sigma_for);
    }

    return NULL;
}