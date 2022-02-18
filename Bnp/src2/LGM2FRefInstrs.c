
#include "LGM2FRefInstrs.h"

#include "opfnctns.h"
#include "math.h"
#include "num_h_minbfgs.h"
#include "num_h_newton.h"
#include "srt_h_all.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"

#define MAX_ITERATION 50
#define FIXED_TAU 10
#define MIN_MEAN_REVERSION -5
#define MAX_MEAN_REVERSION 5.0
#define MIN_ITER_TAU_DIFF 1e-2
#define MAX_NR_STDEV_FROM_PREV_STRIKE 2.0
#define DAYS_IN_QUARTER 90
#define MAX_STDEV 2
#define FIRST_EX_SIGMA_MIN_SQ 2.0e-5
#define SIGMA_MIN_SQ 1.0e-8
#define SIGMA_MAX_SQ 1.0e-2
#define ONE_HALF 0.5
#define SWAPTION_SQ_VOL_TOL 1e-8
#define CAP_VOL_TOL 1e-5
#define WEAK_CAP_VOL_TOL 1e-4
#define CEV_VOL_SHIFT 1e-2
#define NBR_HERMITE_POINTS 25

typedef struct
{
    Date date_today, *ex_dates, *start_dates, *caplet_end_dates, *swp_pay_dates, last_swp_pay_date,
        cap_theo_end_date, *float_start_dates, *float_end_dates, *float_fixing_dates;

    double *cvg_cap, *cvg_pay, *cvg_first, *df_start, *df_pay, *df_cap, **long_swps_payts,
        *long_swps_strikes, *long_swps_fwd, *long_swps_beta, *long_swps_cev, *long_swps_stdev,
        **short_swps_payts, *short_swps_strikes, *short_swps_fwd, *short_swps_beta, *short_swps_cev,
        *short_swps_pv, *short_swps_vega, *long_swps_pv, *long_swps_vega, *float_coverages,
        *caplets_strikes, *caplets_pv, *caplets_vega, *lgm_cap_vol,
        *float_start_dates_discount_factors, *float_end_dates_discount_factors, hybrid_swpts_pv,
        hybrid_swpts_vega, pv_full_cap, vega_full_cap;

    long *         first_cpn_index, nPay, nEx, ex_index, num_float_dates;
    SrtCompounding float_freq;
    SrtBasisCode   float_basis;

    double *Zeta1, *Zeta12, *Zeta2, *H1, *H2, alpha, gamma, rho;
    double *gi, *Wi, v;

    double      zeta1_min, zeta1_max;
    SRT_Boolean apply_criteria;

    long*      num_short_swp_cpn;
    LGMCalMeth cal_method;

    double *    start_zeta1s, *start_taus;
    SRT_Boolean use_start_point_for_taus, use_start_point_for_zeta1s;

} STATICPARAMS;

static STATICPARAMS SParms;

static double *s_PVs, *s_WeighedPVs, *s_zeta1_ptr, *s_du, *s_dv, *s_optim_lambda,
    **s_long_swpts_cond_prob, **s_short_swpts_cond_prob, *s_levenberg_marquard_data_for_swpcalboot,
    *s_levenberg_marquard_target_for_swpcalboot, *s_levenberg_marquard_weight_for_swpcalboot,
    *s_levenberg_marquard_data_for_calibrate, *s_levenberg_marquard_target_for_calibrate,
    *s_levenberg_marquard_weight_for_calibrate, *s_levenberg_marquard_data_for_fullcap_diag,
    *s_levenberg_marquard_target_for_fullcap_diag, *s_levenberg_marquard_weight_for_fullcap_diag;

static double *s_yr_to_ex_dates_ptr, *s_yr_to_start_ptr, s_yr_to_last_pay_date;

double Zeta1_Func(Date tRef, Date* ex_dates, long nEx, Date date_today, double* Zeta1)
{
    double tRefZeta, W;
    int    i;

    if (tRef > ex_dates[nEx])
    {
        W        = (double)(tRef - ex_dates[nEx]) / (ex_dates[nEx] - ex_dates[nEx - 1]);
        tRefZeta = Zeta1[nEx - 1] + W * (Zeta1[nEx - 1] - Zeta1[nEx - 2]);
    }
    else if (tRef <= ex_dates[nEx])
    {
        i = BasicDichotomie(ex_dates, 1, nEx, tRef);
        /* ex_dates[i]=<tRef<ex_dates[i+1] */
        W = (double)(tRef - ex_dates[i]) / (ex_dates[i + 1] - ex_dates[i]);
        tRefZeta =
            Zeta1[i - 1] + W * (Zeta1[i] - Zeta1[i - 1]); /* Zeta1[i-1] = Zeta1(ex_dates[i]) */
    }

    return tRefZeta;
}

double XL_New_Zeta2_Func(
    Date    tRef,
    Date*   ex_dates,
    long    nEx,
    Date    date_today,
    double* Zeta1,
    double  alpha,
    double  gamma,
    double  rho)
{
    long   i;
    double tRefZeta, *Zeta2;

    Zeta2                = dvector(0, nEx - 1);
    s_yr_to_ex_dates_ptr = dvector(1, nEx);

    for (i = 1; i <= nEx; i++)
        s_yr_to_ex_dates_ptr[i] = (double)(ex_dates[i] - date_today) * YEARS_IN_DAY;

    tRefZeta = New_Zeta2_Func(tRef, ex_dates, nEx, date_today, Zeta1, Zeta2, alpha, gamma, rho);

    if (Zeta2)
        free_dvector(Zeta2, 0, nEx - 1);
    Zeta2 = NULL;
    if (s_yr_to_ex_dates_ptr)
        free_dvector(s_yr_to_ex_dates_ptr, 1, nEx);
    s_yr_to_ex_dates_ptr = NULL;

    return tRefZeta;
}

double XL_New_Zeta12_Func(
    Date    tRef,
    Date*   ex_dates,
    long    nEx,
    Date    date_today,
    double* Zeta1,
    double  alpha,
    double  gamma,
    double  rho)
{
    long   i;
    double tRefZeta, *Zeta12;

    Zeta12               = dvector(0, nEx - 1);
    s_yr_to_ex_dates_ptr = dvector(1, nEx);

    for (i = 1; i <= nEx; i++)
        s_yr_to_ex_dates_ptr[i] = (double)(ex_dates[i] - date_today) * YEARS_IN_DAY;

    tRefZeta = New_Zeta12_Func(tRef, ex_dates, nEx, date_today, Zeta1, Zeta12, alpha, gamma, rho);

    if (Zeta12)
        free_dvector(Zeta12, 0, nEx - 1);
    Zeta12 = NULL;
    if (s_yr_to_ex_dates_ptr)
        free_dvector(s_yr_to_ex_dates_ptr, 1, nEx);
    s_yr_to_ex_dates_ptr = NULL;

    return tRefZeta;
}

double New_Zeta2_Func(
    Date    tRef,
    Date*   ex_dates,
    long    nEx,
    Date    date_today,
    double* Zeta1,
    double* Zeta2,
    double  alpha,
    double  gamma,
    double  rho)
{
    int    i;
    double AlphaSq, Weight, tRefZeta, Nextt, Prevt, t, dt, W, tRefZeta1;

    i = BasicDichotomie(ex_dates, 1, nEx, tRef);
    /* ex_dates[i] =< tRef < ex_dates[i+1] */

    if (i == 1) /* CASE WHERE ex_dates[1] =< tRef < ex_dates[2] */
    {
        AlphaSq = alpha * alpha;
        /* COMPUTE ZETA2(ex_dates[0]) FIRST */
        t        = s_yr_to_ex_dates_ptr[1];
        Weight   = AlphaSq / (2 * gamma) * (exp(2 * gamma * t) - 1) / t;
        tRefZeta = Weight * Zeta1[0];

        if (tRef > ex_dates[1]) /* CASE WHERE tRef > ex_dates[1] */
        {
            Nextt  = s_yr_to_ex_dates_ptr[i + 1];
            Prevt  = s_yr_to_ex_dates_ptr[i];
            dt     = (Nextt - Prevt);
            t      = (double)(tRef - date_today) * YEARS_IN_DAY;
            Weight = AlphaSq / (2 * gamma) * (exp(2 * gamma * t) - exp(2 * gamma * Prevt)) / dt;
            tRefZeta += Weight * (Zeta1[i] - Zeta1[i - 1]);

            return tRefZeta;
        }

        else
            return tRefZeta;
    }
    else /* CASE WHERE tRef >= ex_dates[2] */
    {
        AlphaSq = alpha * alpha;
        /* COMPUTATION OF Zeta2(ex_dates[i]) = Zeta2[i-1]; */
        tRefZeta = Zeta2[i - 2];
        t        = s_yr_to_ex_dates_ptr[i];
        Prevt    = s_yr_to_ex_dates_ptr[i - 1];
        dt       = t - Prevt;
        Weight   = AlphaSq / (2 * gamma) * (exp(2 * gamma * t) - exp(2 * gamma * Prevt)) / dt;
        tRefZeta += Weight * (Zeta1[i - 1] - Zeta1[i - 2]);

        if ((tRef > ex_dates[i]) && (i < nEx)) /* CASE WHERE tRef > ex_dates[i] , i < nEx */
        {
            Nextt  = s_yr_to_ex_dates_ptr[i + 1];
            Prevt  = s_yr_to_ex_dates_ptr[i];
            dt     = (Nextt - Prevt);
            t      = (double)(tRef - date_today) * YEARS_IN_DAY;
            Weight = AlphaSq / (2 * gamma) * (exp(2 * gamma * t) - exp(2 * gamma * Prevt)) / dt;
            tRefZeta += Weight * (Zeta1[i] - Zeta1[i - 1]);

            return tRefZeta;
        }
        else if (tRef > ex_dates[nEx]) /* CASE WHERE tRef > ex_dates[nEx] */
        {
            /* COMPUTE Zeta1(tRef) */
            W         = (double)(tRef - ex_dates[nEx]) / (ex_dates[nEx] - ex_dates[nEx - 1]);
            tRefZeta1 = Zeta1[nEx - 1] + W * (Zeta1[nEx - 1] - Zeta1[nEx - 2]);
            /* COMPUTE Zeta2(tRef) */
            dt     = (double)(tRef - ex_dates[nEx]) * YEARS_IN_DAY;
            Prevt  = (double)(ex_dates[nEx] - date_today) * YEARS_IN_DAY;
            t      = (double)(tRef - date_today) * YEARS_IN_DAY;
            Weight = AlphaSq / (2 * gamma) * (exp(2 * gamma * t) - exp(2 * gamma * Prevt)) / dt;
            tRefZeta += Weight * (tRefZeta1 - Zeta1[nEx - 1]);

            return tRefZeta;
        }
        else
            return tRefZeta;
    }
}

double New_Zeta12_Func(
    Date    tRef,
    Date*   ex_dates,
    long    nEx,
    Date    date_today,
    double* Zeta1,
    double* Zeta12,
    double  alpha,
    double  gamma,
    double  rho)
{
    int    i;
    double AlphaRho, Weight, tRefZeta, Nextt, Prevt, t, dt, W, tRefZeta1;

    i = BasicDichotomie(ex_dates, 1, nEx, tRef);
    /* ex_dates[i] =< tRef < ex_dates[i+1] */

    if (i == 1) /* CASE WHERE ex_dates[1] =< tRef < ex_dates[2] */
    {
        AlphaRho = alpha * rho;
        /* COMPUTE ZETA12(ex_dates[0]) FIRST */
        t        = s_yr_to_ex_dates_ptr[1];
        Weight   = (AlphaRho / gamma) * (exp(gamma * t) - 1) / t;
        tRefZeta = Weight * Zeta1[0];

        if (tRef > ex_dates[1]) /* CASE WHERE tRef > ex_dates[1] */
        {
            Nextt  = s_yr_to_ex_dates_ptr[i + 1];
            Prevt  = s_yr_to_ex_dates_ptr[i];
            dt     = (Nextt - Prevt);
            t      = (double)(tRef - date_today) * YEARS_IN_DAY;
            Weight = (AlphaRho / gamma) * (exp(gamma * t) - exp(gamma * Prevt)) / dt;
            tRefZeta += Weight * (Zeta1[i] - Zeta1[i - 1]);

            return tRefZeta;
        }

        else
            return tRefZeta;
    }
    else /* CASE WHERE tRef >= ex_dates[2] */
    {
        AlphaRho = alpha * rho;
        /* COMPUTATION OF Zeta12(ex_dates[i]) = Zeta12[i-1]; */
        tRefZeta = Zeta12[i - 2];
        t        = s_yr_to_ex_dates_ptr[i];
        Prevt    = s_yr_to_ex_dates_ptr[i - 1];
        dt       = t - Prevt;
        Weight   = (AlphaRho / gamma) * (exp(gamma * t) - exp(gamma * Prevt)) / dt;
        tRefZeta += Weight * (Zeta1[i - 1] - Zeta1[i - 2]);

        if (tRef > ex_dates[i]) /* CASE WHERE tRef > ex_dates[i] , i < nEx */
        {
            Nextt  = s_yr_to_ex_dates_ptr[i + 1];
            Prevt  = s_yr_to_ex_dates_ptr[i];
            dt     = (Nextt - Prevt);
            t      = (double)(tRef - date_today) * YEARS_IN_DAY;
            Weight = (AlphaRho / gamma) * (exp(gamma * t) - exp(gamma * Prevt)) / dt;
            tRefZeta += Weight * (Zeta1[i] - Zeta1[i - 1]);

            return tRefZeta;
        }
        else if (tRef > ex_dates[nEx]) /* CASE WHERE tRef > ex_dates[nEx] */
        {
            /* COMPUTE Zeta1(tRef) */
            W         = (double)(tRef - ex_dates[nEx]) / (ex_dates[nEx] - ex_dates[nEx - 1]);
            tRefZeta1 = Zeta1[nEx - 1] + W * (Zeta1[nEx - 1] - Zeta1[nEx - 2]);
            /* COMPUTE Zeta2(tRef) */
            dt     = (double)(tRef - ex_dates[nEx]) * YEARS_IN_DAY;
            Prevt  = (double)(ex_dates[nEx] - date_today) * YEARS_IN_DAY;
            t      = (double)(tRef - date_today) * YEARS_IN_DAY;
            Weight = (AlphaRho / gamma) * (exp(gamma * t) - exp(gamma * Prevt)) / dt;
            tRefZeta += Weight * (tRefZeta1 - Zeta1[nEx - 1]);

            return tRefZeta;
        }
        else
            return tRefZeta;
    }
}

/* BETWEEN TWO PAY DATES H1 IS LINEARLY INTERPOLATED */
double H1_Func(Date tRef, Date* start_dates, Date date_today, long nEx, double* H1)
{
    double tRefH, W;
    long   i;
    Date   last_swp_pay_date;

    last_swp_pay_date = SParms.last_swp_pay_date;

    /* start_dates[nEx] IS THE START DATE OF THE LAST SWAPTION && last_swp_pay_date IS THE LAST
     * PAYMENT DATE */

    if ((tRef >= start_dates[nEx]) && (tRef <= last_swp_pay_date))
    {
        W     = (double)(tRef - start_dates[nEx]) / (last_swp_pay_date - start_dates[nEx]);
        tRefH = (1 - W) * H1[nEx - 1];

        return tRefH;
    }
    else if ((tRef < start_dates[nEx]) && (tRef >= start_dates[1]))
    {
        i = BasicDichotomie(start_dates, 1, nEx, tRef);
        /* start_dates[i] =< tRef < start_dates[i+1] */
        W     = (double)(tRef - start_dates[i]) / (start_dates[i + 1] - start_dates[i]);
        tRefH = H1[i - 1] + W * (H1[i] - H1[i - 1]);
    }
    else if (tRef < start_dates[1])
    {
        W     = (double)(start_dates[1] - tRef) / (start_dates[2] - start_dates[1]);
        tRefH = H1[0] + W * (H1[0] - H1[1]);
    }

    return tRefH;
}

double XLH1_Func(
    Date tRef, Date* start_dates, Date date_today, Date last_swp_pay_date, long nEx, double* H1)
{
    double tRefH1;

    SParms.last_swp_pay_date = last_swp_pay_date;

    tRefH1 = H1_Func(tRef, start_dates, date_today, nEx, H1);

    return tRefH1;
}

double XLH2_Func(
    Date    tRef,
    Date*   start_dates,
    Date    date_today,
    Date    last_swp_pay_date,
    long    nEx,
    double* H1,
    double  gamma)
{
    double tRefH2;

    SParms.last_swp_pay_date = last_swp_pay_date;

    tRefH2 = H2_Func(tRef, start_dates, date_today, nEx, H1, gamma);

    return tRefH2;
}

double New_H2_Func(
    Date tRef, Date* start_dates, Date date_today, long nEx, double* H1, double* H2, double gamma)
{
    double tRefH, last_swp_pay_date = SParms.last_swp_pay_date;
    double Prevt, Nextt, t, dt, Weight;
    long   j;

    if (tRef == last_swp_pay_date)
    {
        tRefH = 0.;
        return tRefH;
    }
    else if ((tRef < last_swp_pay_date) && (tRef > start_dates[nEx]))
    {
        Prevt  = s_yr_to_start_ptr[nEx];
        Nextt  = s_yr_to_last_pay_date;
        dt     = (Nextt - Prevt);
        t      = (double)(tRef - date_today) * YEARS_IN_DAY;
        Weight = (exp(-gamma * Nextt) - exp(-gamma * t)) / (gamma * dt);
        tRefH  = -Weight * H1[nEx - 1];

        return tRefH;
    }
    else if ((tRef < last_swp_pay_date) && (tRef <= start_dates[nEx]))
    {
        tRefH = 0.;
        if (tRef == start_dates[nEx])
        {
            Prevt  = s_yr_to_start_ptr[nEx];
            Nextt  = s_yr_to_last_pay_date;
            dt     = (Nextt - Prevt);
            Weight = (exp(-gamma * Nextt) - exp(-gamma * Prevt)) / (gamma * dt);
            tRefH += -Weight * H1[nEx - 1];

            return tRefH;
        }
        else if (tRef < start_dates[nEx])
        {
            j = BasicDichotomie(start_dates, 1, nEx, tRef);
            /* start_dates[j] =< tRef < start_dates[j+1]  */

            tRefH  = H2[nEx - 1];
            Prevt  = s_yr_to_start_ptr[j];
            t      = s_yr_to_start_ptr[j + 1];
            dt     = (t - Prevt);
            Weight = (exp(-gamma * t) - exp(-gamma * Prevt)) / (gamma * dt);
            tRefH += Weight * (H1[j] - H1[j - 1]);

            if (tRef > start_dates[j])
            {
                Prevt  = s_yr_to_start_ptr[j];
                Nextt  = s_yr_to_start_ptr[j + 1];
                dt     = (Nextt - Prevt);
                t      = (double)(tRef - date_today) * YEARS_IN_DAY;
                Weight = (exp(-gamma * Prevt) - exp(-gamma * t)) / (gamma * dt);
                tRefH += Weight * (H1[j] - H1[j - 1]);

                return tRefH;
            }
            else
                return tRefH;
        }
    }
    return tRefH;
}

double H2_Func(Date tRef, Date* start_dates, Date date_today, long nEx, double* H1, double gamma)
{
    double tRefH; /* RETURN */
    double t, Prevt, Nextt, dt;
    long   j, i;
    double Weight;
    Date   last_swp_pay_date;

    last_swp_pay_date = SParms.last_swp_pay_date;
    tRefH             = 0.;

    if (tRef == last_swp_pay_date)
    {
        tRefH = 0;
        return tRefH;
    }
    else if ((tRef < last_swp_pay_date) && (tRef > start_dates[nEx]))
    {
        Prevt  = (double)(start_dates[nEx] - date_today) * YEARS_IN_DAY;
        Nextt  = (double)(last_swp_pay_date - date_today) * YEARS_IN_DAY;
        dt     = (Nextt - Prevt);
        t      = (double)(tRef - date_today) * YEARS_IN_DAY;             /* Prevt < t < Nextt */
        Weight = (exp(-gamma * Nextt) - exp(-gamma * t)) / (gamma * dt); /* Weight > 0  */
        tRefH  = -Weight * H1[nEx - 1];

        return tRefH;
    }
    else if ((tRef < last_swp_pay_date) && (tRef <= start_dates[nEx]))
    {
        j = BasicDichotomie(start_dates, 1, nEx, tRef);

        /* start_dates[j] =< tRef < start_dates[j+1]  */

        tRefH = 0.;

        if (j == nEx)
        {
            Prevt  = (double)(start_dates[nEx] - date_today) * YEARS_IN_DAY;
            Nextt  = (double)(last_swp_pay_date - date_today) * YEARS_IN_DAY;
            dt     = (Nextt - Prevt);
            Weight = (exp(-gamma * Nextt) - exp(-gamma * Prevt)) / (gamma * dt);
            tRefH += -Weight * H1[nEx - 1];
        }

        else if (j < nEx)
        {
            Prevt  = (double)(start_dates[nEx] - date_today) * YEARS_IN_DAY;
            Nextt  = (double)(last_swp_pay_date - date_today) * YEARS_IN_DAY;
            dt     = (Nextt - Prevt);
            Weight = (exp(-gamma * Nextt) - exp(-gamma * Prevt)) / (gamma * dt);
            tRefH += -Weight * H1[nEx - 1];

            for (i = (j + 1); i <= nEx; i++) /* compute H2(start_dates[j]) */
            {
                Prevt  = (double)(start_dates[i - 1] - date_today) * YEARS_IN_DAY;
                t      = (double)(start_dates[i] - date_today) * YEARS_IN_DAY;
                dt     = (t - Prevt);
                Weight = (exp(-gamma * t) - exp(-gamma * Prevt)) / (gamma * dt);
                tRefH += Weight * (H1[i - 1] - H1[i - 2]);
            }

            if (tRef > start_dates[j])
            {
                Prevt  = (double)(start_dates[j] - date_today) * YEARS_IN_DAY;
                Nextt  = (double)(start_dates[j + 1] - date_today) * YEARS_IN_DAY;
                dt     = (Nextt - Prevt);
                t      = (double)(tRef - date_today) * YEARS_IN_DAY; /* Prevt < t < Nextt  */
                Weight = (exp(-gamma * Prevt) - exp(-gamma * t)) / (gamma * dt); /* Weight > 0 */
                tRefH += Weight * (H1[j] - H1[j - 1]);
            }
        }
    }

    return tRefH;
}

Err XLLGMZeta1Func(double t, char* UndName, double* Zeta)
{
    Err         err = NULL;
    TermStruct* ts;
    SrtTFTSMat  G;
    SrtUndPtr   sUndPtr;
    SrtMdlType  mdl_type;
    long        mdl_dim;
    double      H;

    sUndPtr = lookup_und(UndName);

    err = get_underlying_mdldim(sUndPtr, &mdl_dim);
    if (err)
        return err;

    err = get_underlying_mdltype(sUndPtr, &mdl_type);
    if (err)
        return err;

    err = get_underlying_ts(sUndPtr, &ts);
    if (err)
        return err;

    if (mdl_dim == TWO_FAC)
    {
        err = get_2f_G_funcs(t, ts, &G);
        if (err)
            return err;

        *Zeta = G[0][0];
    }
    else
    {
        G_H_func(t, ts, Zeta, &H);
    }

    return err;
}

Err XLLGMZeta2Func(double t, char* UndName, double* Zeta)
{
    Err         err = NULL;
    TermStruct* ts;
    SrtTFTSMat  G;
    SrtUndPtr   sUndPtr;

    sUndPtr = lookup_und(UndName);
    err     = get_underlying_ts(sUndPtr, &ts);
    if (err)
        return err;

    err = get_2f_G_funcs(t, ts, &G);
    if (err)
        return err;

    *Zeta = G[1][1];

    return err;
}

static Err srt_f_autocal_df(
    Date    swp_pay_dates,
    double  df_pay,
    double  U,
    double  V,
    double  zeta1_ex,
    double  zeta2_ex,
    double  zeta12_ex,
    double  H1tPay,
    double  H2tPay,
    double* ZC)
{
    double X1tEx, X2tEx;
    double RootZeta1tEx, RootZeta2tEx;
    double RHO;
    Err    err = NULL;

    if ((zeta1_ex > 0.0) && (zeta2_ex > 0.0))
    {
        RootZeta1tEx = sqrt(zeta1_ex);
        RootZeta2tEx = sqrt(zeta2_ex);

        RHO = zeta12_ex / (RootZeta1tEx * RootZeta2tEx);
        if (fabs(RHO) > 1.0)
            return serror("Fatal error in srt_f_autocal_df");

        X1tEx = RootZeta1tEx * U;
        X2tEx = RootZeta2tEx * (RHO * U + sqrt(1 - RHO * RHO) * V);

        (*ZC) = df_pay * exp(H1tPay * X1tEx + H2tPay * X2tEx - H1tPay * H2tPay * zeta12_ex -
                             0.5 * (H1tPay * H1tPay * zeta1_ex + H2tPay * H2tPay * zeta2_ex));
    }
    else
        return serror("Fatal error in srt_f_autocal_df");

    return err;
}

static Err srt_f_autocal_swap_pv(
    long    ex_index,
    long    num_cpn,
    double  Rf,
    double  U,
    double  V,
    double  zeta1_ex,
    double  zeta2_ex,
    double  zeta12_ex,
    double* H1,
    double* SWAP)
{
    Date *start_dates = SParms.start_dates, *swp_pay_dates = SParms.swp_pay_dates,
         date_today = SParms.date_today;
    long    first_cpn_index, nEx = SParms.nEx, k;
    double *cvg_pay, ak, cvg_first, *df_start = SParms.df_start, *df_pay = SParms.df_pay, H1tStart,
                                    H2tStart, H1tPay, H2tPay, ZCtStart, ZCtPayk,
                                    gamma = SParms.gamma;
    Err err                               = NULL;

    (*SWAP) = 0.0;

    H1tStart = H1_Func(start_dates[ex_index], start_dates, date_today, nEx, H1);
    H2tStart = H2_Func(start_dates[ex_index], start_dates, date_today, nEx, H1, gamma);

    err = srt_f_autocal_df(
        start_dates[ex_index],
        df_start[ex_index],
        U,
        V,
        zeta1_ex,
        zeta2_ex,
        zeta12_ex,
        H1tStart,
        H2tStart,
        &ZCtStart);
    if (err)
        return err;

    (*SWAP) -= ZCtStart;

    first_cpn_index = SParms.first_cpn_index[ex_index];
    cvg_first       = SParms.cvg_first[ex_index];
    cvg_pay         = SParms.cvg_pay;

    for (k = 0; k <= num_cpn; k++)
    {
        ak     = 0.;
        H1tPay = H1_Func(swp_pay_dates[first_cpn_index + k], start_dates, date_today, nEx, H1);
        H2tPay =
            H2_Func(swp_pay_dates[first_cpn_index + k], start_dates, date_today, nEx, H1, gamma);

        err = srt_f_autocal_df(
            swp_pay_dates[first_cpn_index + k],
            df_pay[first_cpn_index + k],
            U,
            V,
            zeta1_ex,
            zeta2_ex,
            zeta12_ex,
            H1tPay,
            H2tPay,
            &ZCtPayk);
        if (err)
            return err;

        if ((k == 0) && (num_cpn > 0))
            ak = cvg_first * Rf;
        if ((k > 0) && (k < num_cpn))
            ak = cvg_pay[first_cpn_index + k - 1] * Rf;

        if (k == num_cpn)
        {
            if (num_cpn > 0)
                ak = (1 + cvg_pay[first_cpn_index + k - 1] * Rf);
            else if (num_cpn == 0)
                ak = (1 + cvg_first * Rf);
        }

        (*SWAP) += ak * ZCtPayk;
    }

    return err;
}

Err make_autocal_cond_prob(
    double*   gi,
    double*   Wi,
    double    du_s,
    double*   s_du,
    long      jEx,
    long      num_cpn,
    double    Rf,
    double*   Zeta1,
    double*   Zeta2,
    double*   Zeta12,
    double*   H1,
    long      nher,
    double*** cond_prob)
{
    Err         err = NULL;
    double      a[3], b[3], nstop, d;
    long        l, k, i;
    SRT_Boolean is_decreasing;

    for (l = 1; l <= nher; l++)
    {
        nstop = 0.0;
        a[0]  = 1;
        srt_f_autocal_swap_pv(
            jEx,
            num_cpn,
            Rf,
            a[0],
            gi[l],
            Zeta1[jEx - 1],
            Zeta2[jEx - 1],
            Zeta12[jEx - 1],
            H1,
            &b[0]);
        a[1] = 1.01;
        srt_f_autocal_swap_pv(
            jEx,
            num_cpn,
            Rf,
            a[1],
            gi[l],
            Zeta1[jEx - 1],
            Zeta2[jEx - 1],
            Zeta12[jEx - 1],
            H1,
            &b[1]);

        if (b[0] > b[1])
            is_decreasing = SRT_YES;
        else
            is_decreasing = SRT_NO;

        a[2] = 1.02;

        k = 0;
        while (nstop < 1.0 && k < MAX_ITERATION)
        {
            srt_f_autocal_swap_pv(
                jEx,
                num_cpn,
                Rf,
                a[2],
                gi[l],
                Zeta1[jEx - 1],
                Zeta2[jEx - 1],
                Zeta12[jEx - 1],
                H1,
                &b[2]);
            newton(0.0, 2.0, a, b, &nstop);
            k++;
            /* sets an arbitary bound to prevent overflow */
            if (a[2] - du_s > 30.)
                a[2] = 30 + du_s;
            else if (a[2] - du_s < -30.)
                a[2] = -30 + du_s;
        }

        d = a[2];

        if (is_decreasing == SRT_YES)
        {
            (*cond_prob)[0][l] = norm(d - du_s);
            for (i = 1; i <= (num_cpn + 1); i++)
                (*cond_prob)[i][l] = norm(d - s_du[i - 1]);
        }
        else
        {
            (*cond_prob)[0][l] = (1 - norm(d - du_s));
            for (i = 1; i <= (num_cpn + 1); i++)
                (*cond_prob)[i][l] = (1 - norm(d - s_du[i - 1]));
        }
    }
    return (err);
}

static Err make_static_allocations(long nEx, long nPay)
{
    Err err = NULL;

    s_PVs                   = dvector(1, nEx + nEx);
    s_WeighedPVs            = dvector(1, 2 * nEx);
    s_zeta1_ptr             = dvector(1, 1);
    s_yr_to_ex_dates_ptr    = dvector(1, nEx);
    s_yr_to_start_ptr       = dvector(1, nEx);
    s_du                    = dvector(0, nPay - 1);
    s_dv                    = dvector(0, nPay - 1);
    s_long_swpts_cond_prob  = dmatrix(0, nPay, 1, NBR_HERMITE_POINTS);
    s_short_swpts_cond_prob = dmatrix(0, nPay, 1, NBR_HERMITE_POINTS);
    s_optim_lambda          = dvector(1, 1);

    if ((s_PVs == NULL) || (s_WeighedPVs == NULL) || (s_zeta1_ptr == NULL) ||
        (s_yr_to_ex_dates_ptr == NULL) || (s_yr_to_start_ptr == NULL) || (s_du == NULL) ||
        (s_dv == NULL) || (s_long_swpts_cond_prob == NULL) || (s_short_swpts_cond_prob == NULL))
        return serror("Allocation Memory Failure");

    return err;
}

static Err make_static_desallocations(long nEx, long nPay)
{
    Err err = NULL;

    if (nEx == 0)
        return serror("No Exercise in make_static_desallocations");

    if (s_PVs)
        free_dvector(s_PVs, 1, 2 * nEx);
    s_PVs = NULL;
    if (s_WeighedPVs)
        free_dvector(s_WeighedPVs, 1, 2 * nEx);
    s_WeighedPVs = NULL;
    if (s_zeta1_ptr)
        free_dvector(s_zeta1_ptr, 1, 1);
    s_zeta1_ptr = NULL;
    if (s_yr_to_ex_dates_ptr)
        free_dvector(s_yr_to_ex_dates_ptr, 1, nEx);
    s_yr_to_ex_dates_ptr = NULL;
    if (s_yr_to_start_ptr)
        free_dvector(s_yr_to_start_ptr, 1, nEx);
    s_yr_to_start_ptr = NULL;

    if (s_du)
        free_dvector(s_du, 0, nPay - 1);
    s_du = NULL;
    if (s_dv)
        free_dvector(s_dv, 0, nPay - 1);
    s_dv = NULL;
    if (s_long_swpts_cond_prob)
        free_dmatrix(s_long_swpts_cond_prob, 0, nPay, 1, NBR_HERMITE_POINTS);
    s_long_swpts_cond_prob = NULL;
    if (s_long_swpts_cond_prob)
        free_dmatrix(s_short_swpts_cond_prob, 0, nPay, 1, NBR_HERMITE_POINTS);
    s_short_swpts_cond_prob = NULL;
    if (s_optim_lambda)
        free_dvector(s_optim_lambda, 1, 1);
    s_optim_lambda = NULL;

    return err;
}

static Err make_levenberg_struct(long nEx)
{
    Err  err = NULL;
    long i;

    /*FOR THE BOOTSTRAP ON DIAG SWAPTIONS */
    s_levenberg_marquard_data_for_swpcalboot   = dvector(1, 1);
    s_levenberg_marquard_target_for_swpcalboot = dvector(1, 1);
    s_levenberg_marquard_weight_for_swpcalboot = dvector(1, 1);

    /* FOR THE GLOBAL CALIBRATION ON CAPLETS AND DIAG SWAPTIONS */
    s_levenberg_marquard_data_for_calibrate   = dvector(1, nEx + nEx);
    s_levenberg_marquard_target_for_calibrate = dvector(1, nEx + nEx);
    s_levenberg_marquard_weight_for_calibrate = dvector(1, nEx + nEx);

    /* FOR THE GLOBAL CALIBRATION ON FULL CAP AND DIAG SWAPTIONS */
    s_levenberg_marquard_data_for_fullcap_diag   = dvector(1, nEx + 1);
    s_levenberg_marquard_target_for_fullcap_diag = dvector(1, nEx + 1);
    s_levenberg_marquard_weight_for_fullcap_diag = dvector(1, nEx + 1);

    if ((s_levenberg_marquard_data_for_swpcalboot == NULL) ||
        (s_levenberg_marquard_target_for_swpcalboot == NULL) ||
        (s_levenberg_marquard_weight_for_swpcalboot == NULL) ||
        (s_levenberg_marquard_data_for_calibrate == NULL) ||
        (s_levenberg_marquard_target_for_calibrate == NULL) ||
        (s_levenberg_marquard_weight_for_calibrate == NULL) ||
        (s_levenberg_marquard_data_for_fullcap_diag == NULL) ||
        (s_levenberg_marquard_target_for_fullcap_diag == NULL) ||
        (s_levenberg_marquard_weight_for_fullcap_diag == NULL))

        return serror("Memory Allocation Failure");

    s_levenberg_marquard_data_for_swpcalboot[1] = 1;

    for (i = 1; i <= nEx; i++)
    {
        s_levenberg_marquard_data_for_calibrate[i]   = (double)i;
        s_levenberg_marquard_target_for_calibrate[i] = SParms.short_swps_pv[i];
        s_levenberg_marquard_weight_for_calibrate[i] = 1.0;
    }

    for (i = nEx + 1; i <= nEx + nEx; i++)
    {
        s_levenberg_marquard_data_for_calibrate[i]   = (double)i;
        s_levenberg_marquard_target_for_calibrate[i] = SParms.long_swps_pv[i - nEx];
        s_levenberg_marquard_weight_for_calibrate[i] = SParms.long_swps_vega[i - nEx];
    }

    for (i = 1; i <= nEx; i++)
    {
        s_levenberg_marquard_data_for_fullcap_diag[i]   = (double)i;
        s_levenberg_marquard_target_for_fullcap_diag[i] = SParms.long_swps_pv[i];
        s_levenberg_marquard_weight_for_fullcap_diag[i] = SParms.long_swps_vega[i];
    }

    s_levenberg_marquard_data_for_fullcap_diag[nEx + 1]   = (double)(nEx + 1);
    s_levenberg_marquard_target_for_fullcap_diag[nEx + 1] = SParms.pv_full_cap;
    s_levenberg_marquard_weight_for_fullcap_diag[nEx + 1] = SParms.vega_full_cap;

    return err;
}

static Err delete_levenberg_struct(long nEx)
{
    Err err = NULL;

    if (s_levenberg_marquard_data_for_swpcalboot)
        free_dvector(s_levenberg_marquard_data_for_swpcalboot, 1, 1);
    s_levenberg_marquard_data_for_swpcalboot = NULL;
    if (s_levenberg_marquard_target_for_swpcalboot)
        free_dvector(s_levenberg_marquard_target_for_swpcalboot, 1, 1);
    s_levenberg_marquard_target_for_swpcalboot = NULL;
    if (s_levenberg_marquard_weight_for_swpcalboot)
        free_dvector(s_levenberg_marquard_target_for_swpcalboot, 1, 1);
    s_levenberg_marquard_weight_for_swpcalboot = NULL;

    if (s_levenberg_marquard_data_for_calibrate)
        free_dvector(s_levenberg_marquard_data_for_calibrate, 1, nEx + nEx);
    s_levenberg_marquard_data_for_calibrate = NULL;
    if (s_levenberg_marquard_target_for_calibrate)
        free_dvector(s_levenberg_marquard_target_for_calibrate, 1, nEx + nEx);
    s_levenberg_marquard_target_for_calibrate = NULL;
    if (s_levenberg_marquard_weight_for_calibrate)
        free_dvector(s_levenberg_marquard_weight_for_calibrate, 1, nEx + nEx);
    s_levenberg_marquard_weight_for_calibrate = NULL;

    if (s_levenberg_marquard_data_for_fullcap_diag)
        free_dvector(s_levenberg_marquard_data_for_fullcap_diag, 1, nEx + 1);
    s_levenberg_marquard_data_for_fullcap_diag = NULL;
    if (s_levenberg_marquard_target_for_fullcap_diag)
        free_dvector(s_levenberg_marquard_target_for_fullcap_diag, 1, nEx + 1);
    s_levenberg_marquard_target_for_fullcap_diag = NULL;
    if (s_levenberg_marquard_weight_for_fullcap_diag)
        free_dvector(s_levenberg_marquard_weight_for_fullcap_diag, 1, nEx + 1);
    s_levenberg_marquard_weight_for_fullcap_diag = NULL;

    return err;
}

static Err srt_f_autocal_caplet_greeks(long jEx, double zeta1_ex, double* caplet_pv)
{
    Err     err = NULL;
    double *H1 = SParms.H1, *H2 = SParms.H2, *Zeta1 = SParms.Zeta1, *Zeta2 = SParms.Zeta2,
           *Zeta12 = SParms.Zeta12, *df_start = SParms.df_start, *df_cap = SParms.df_cap,
           *cvg_cap = SParms.cvg_cap;
    Date *ex_dates = SParms.ex_dates, date_today = SParms.date_today,
         *start_dates = SParms.start_dates, *caplet_end_dates = SParms.caplet_end_dates;
    double zeta2_ex, zeta12_ex, H1tStart, H1tEnd, H2tStart, H2tEnd, cap_vol, sq_cap_vol,
        gamma = SParms.gamma, t;
    long nEx  = SParms.nEx;

    t = (double)(ex_dates[jEx] - date_today) * YEARS_IN_DAY;

    Zeta1[jEx - 1] = zeta1_ex;
    zeta2_ex       = Zeta2[jEx - 1];
    zeta12_ex      = Zeta12[jEx - 1];

    H1tStart = H1[jEx - 1];
    H1tEnd   = H1_Func(caplet_end_dates[jEx], start_dates, date_today, nEx, H1);

    H2tStart = H2[jEx - 1];
    H2tEnd   = H2_Func(caplet_end_dates[jEx], start_dates, date_today, nEx, H1, gamma);

    sq_cap_vol = ((H1tEnd - H1tStart) * (H1tEnd - H1tStart) * zeta1_ex +
                  (H2tEnd - H2tStart) * (H2tEnd - H2tStart) * zeta2_ex +
                  2 * (H1tEnd - H1tStart) * (H2tEnd - H2tStart) * zeta12_ex) /
                 t;

    if (sq_cap_vol == 0)
        return serror("Error while computing caplet volatility");

    cap_vol = sqrt(sq_cap_vol);

    (*caplet_pv) = 0.;
    (*caplet_pv) = srt_f_optblksch(
        df_start[jEx] / df_cap[jEx],
        (1 + cvg_cap[jEx] * SParms.short_swps_strikes[jEx]),
        cap_vol,
        t,
        df_cap[jEx],
        SRT_PUT,
        SRT_PREMIUM);

    return err;
}

static Err srt_f_autocal_swpts_greeks(
    long jEx, long num_cpn, double Rf, double** SwapPayts, double zeta1_ex, double* swption_pv)
{
    Err     err   = NULL;
    double *Zeta1 = SParms.Zeta1, *Zeta2 = SParms.Zeta2, *Zeta12 = SParms.Zeta12, *H1 = SParms.H1,
           *H2 = SParms.H2, *gi = SParms.gi, *Wi = SParms.Wi, alpha = SParms.alpha,
           gamma = SParms.gamma, rho = SParms.rho, *df_start = SParms.df_start, RootZeta1tEx,
           zeta2_ex, RootZeta2tEx, zeta12_ex, RHO, H1tStart, H2tStart, duStart, dvStart, dvStartSq,
           H1tk, H2tk, SwpPayoff;

    Date *ex_dates = SParms.ex_dates, *start_dates = SParms.start_dates,
         date_today = SParms.date_today, *swp_pay_dates = SParms.swp_pay_dates;
    long nPay = SParms.nPay, nEx = SParms.nEx, *first_cpn_index = SParms.first_cpn_index, k, l,
         PayIndex;

    Zeta1[jEx - 1] = zeta1_ex;
    zeta2_ex =
        New_Zeta2_Func(ex_dates[jEx], ex_dates, nEx, date_today, Zeta1, Zeta2, alpha, gamma, rho);
    Zeta2[jEx - 1] = zeta2_ex;
    zeta12_ex =
        New_Zeta12_Func(ex_dates[jEx], ex_dates, nEx, date_today, Zeta1, Zeta12, alpha, gamma, rho);
    Zeta12[jEx - 1] = zeta12_ex;
    H1tStart        = H1[jEx - 1];
    H2tStart        = H2[jEx - 1];

    (*swption_pv) = 0.0;

    RootZeta1tEx = sqrt(zeta1_ex);
    RootZeta2tEx = sqrt(zeta2_ex);
    RHO          = zeta12_ex / (RootZeta1tEx * RootZeta2tEx);

    /* START DATE CASH FLOW */
    duStart   = H1tStart * RootZeta1tEx + RHO * H2tStart * RootZeta2tEx;
    dvStart   = sqrt(1 - RHO * RHO) * H2tStart * RootZeta2tEx;
    dvStartSq = dvStart * dvStart;

    for (k = 0; k <= num_cpn; k++)
    {
        PayIndex = first_cpn_index[jEx] + k;
        H1tk     = H1_Func(swp_pay_dates[PayIndex], start_dates, date_today, nEx, H1);
        H2tk     = H2_Func(swp_pay_dates[PayIndex], start_dates, date_today, nEx, H1, gamma);

        s_du[k] = H1tk * RootZeta1tEx + RHO * H2tk * RootZeta2tEx;
        s_dv[k] = sqrt(1 - RHO * RHO) * H2tk * RootZeta2tEx;
    }

    err = make_autocal_cond_prob(
        gi,
        Wi,
        duStart,
        s_du,
        jEx,
        num_cpn,
        Rf,
        Zeta1,
        Zeta2,
        Zeta12,
        H1,
        NBR_HERMITE_POINTS,
        &s_long_swpts_cond_prob);
    if (err)
        return err;

    for (l = 1; l <= NBR_HERMITE_POINTS; l++)
    {
        SwpPayoff = 0.0; /* RESET OF PAYOFF VALUES */
        SwpPayoff -=
            df_start[jEx] * s_long_swpts_cond_prob[0][l] * exp(gi[l] * dvStart - 0.5 * dvStartSq);

        for (k = 0; k <= num_cpn; k++)
            SwpPayoff += SwapPayts[jEx][k] * s_long_swpts_cond_prob[k + 1][l] *
                         exp(gi[l] * s_dv[k] - 0.5 * s_dv[k] * s_dv[k]);

        (*swption_pv) += Wi[l] * SwpPayoff;

    } /* END OF LOOP ON l: HERMITE POINTS NBR */

    return err;
}

static Err srt_f_autocal_cap_greeks(double* cap_pv)
{
    Err   err               = NULL;
    Date *float_start_dates = SParms.float_start_dates, *float_end_dates = SParms.float_end_dates,
         *float_fixing_dates = SParms.float_fixing_dates, *start_dates = SParms.start_dates,
         *ex_dates = SParms.ex_dates, date_today = SParms.date_today;
    double *H1 = SParms.H1, *Zeta1 = SParms.Zeta1, *Zeta2 = SParms.Zeta2, *Zeta12 = SParms.Zeta12,
           alpha = SParms.alpha, gamma = SParms.gamma, rho = SParms.rho,
           *float_start_dates_discount_factors = SParms.float_start_dates_discount_factors,
           *float_end_dates_discount_factors   = SParms.float_end_dates_discount_factors,
           *float_coverages = SParms.float_coverages, zeta1_ex, zeta2_ex, zeta12_ex, H1tStart,
           H1tEnd, H2tStart, H2tEnd, DtStart, DtEnd, cap_vol, sq_cap_vol, t, caplet_pv, CvgFloat;
    long i, j, nEx = SParms.nEx, num_float_dates = SParms.num_float_dates;

    (*cap_pv) = 0.;

    for (j = 1; j <= num_float_dates; j++)
    {
        i = BasicDichotomie(start_dates, 1, nEx, float_start_dates[j - 1]);

        caplet_pv = 0.;
        H1tStart  = H1_Func(float_start_dates[j - 1], start_dates, date_today, nEx, H1);
        H1tEnd    = H1_Func(float_end_dates[j - 1], start_dates, date_today, nEx, H1);

        H2tStart = H2_Func(float_start_dates[j - 1], start_dates, date_today, nEx, H1, gamma);
        H2tEnd   = H2_Func(float_end_dates[j - 1], start_dates, date_today, nEx, H1, gamma);

        t = (double)(float_fixing_dates[j - 1] - date_today) * YEARS_IN_DAY;

        zeta1_ex = Zeta1_Func(float_fixing_dates[j - 1], ex_dates, nEx, date_today, Zeta1);
        zeta2_ex = New_Zeta2_Func(
            float_fixing_dates[j - 1], ex_dates, nEx, date_today, Zeta1, Zeta2, alpha, gamma, rho);
        zeta12_ex = New_Zeta12_Func(
            float_fixing_dates[j - 1], ex_dates, nEx, date_today, Zeta1, Zeta12, alpha, gamma, rho);

        sq_cap_vol = ((H1tEnd - H1tStart) * (H1tEnd - H1tStart) * zeta1_ex +
                      (H2tEnd - H2tStart) * (H2tEnd - H2tStart) * zeta2_ex +
                      2 * (H1tEnd - H1tStart) * (H2tEnd - H2tStart) * zeta12_ex) /
                     t;

        cap_vol = sqrt(sq_cap_vol);

        DtStart = float_start_dates_discount_factors[j - 1];
        DtEnd   = float_end_dates_discount_factors[j - 1];

        if (num_float_dates > 1)
            CvgFloat = float_coverages[j - 1];
        else
            CvgFloat =
                coverage(float_start_dates[j - 1], float_end_dates[j - 1], SParms.float_basis);

        caplet_pv = srt_f_optblksch(
            DtStart / DtEnd,
            (1 + CvgFloat * SParms.short_swps_strikes[i]),
            cap_vol,
            t,
            DtEnd,
            SRT_PUT,
            SRT_PREMIUM);

        (*cap_pv) += caplet_pv;
    }

    return err;
}

static Err srt_f_autocal_hybrid_swpts_greeks(double* hybrid_swpts_pv)
{
    Err      err = NULL;
    long     nEx = SParms.nEx, i, *num_short_swp_cpn = SParms.num_short_swp_cpn;
    double **short_swps_payts = SParms.short_swps_payts, *Zeta1 = SParms.Zeta1, hybrid_swaplet_pv;
    double   Rf;

    (*hybrid_swpts_pv) = 0.0;
    for (i = 1; i < nEx; i++)
    {
        Rf  = SParms.short_swps_strikes[i];
        err = srt_f_autocal_swpts_greeks(
            i, num_short_swp_cpn[i], Rf, short_swps_payts, Zeta1[i - 1], &hybrid_swaplet_pv);
        if (err)
            return err;

        (*hybrid_swpts_pv) += hybrid_swaplet_pv;
    }

    return err;
}

static Err levenberg_func_for_tau(
    double  diag_swp_index,
    double* zeta1_param,
    double* levenberg_value,
    double* levenberg_value_der,
    int     nparm)
{
    Err  err = NULL;
    long jEx = SParms.ex_index, nEx = SParms.nEx, nPay = SParms.nPay,
         *first_cpn_index = SParms.first_cpn_index, num_cpn;
    double *Zeta1 = SParms.Zeta1, *Zeta2 = SParms.Zeta2, *Zeta12 = SParms.Zeta12,
           alpha = SParms.alpha, gamma = SParms.gamma, rho = SParms.rho, zeta1_ex, diag_swaption_pv,
           param_shift, zero_shift = 1e-6, prop_shift = 1e-3, zeta1_min = SParms.zeta1_min,
           zeta1_max = SParms.zeta1_max;
    Date *ex_dates = SParms.ex_dates, *start_dates = SParms.start_dates,
         date_today          = SParms.date_today;
    double** long_swps_payts = SParms.long_swps_payts;
    double   Rf;

    zeta1_ex       = zeta1_param[1];
    Zeta1[jEx - 1] = zeta1_ex;

    (*levenberg_value) = 0.0;
    (*levenberg_value) += max(-Zeta1[0], 0);
    if (jEx > 1)
        (*levenberg_value) += max(-(Zeta1[jEx - 1] - Zeta1[jEx - 2]), 0.0);

    if ((Zeta1[jEx - 1] < zeta1_min) || (Zeta1[jEx - 1] > zeta1_max))
        (*levenberg_value) += 10.0;

    if ((*levenberg_value) == 0)
    {
        Zeta2[jEx - 1] = New_Zeta2_Func(
            ex_dates[jEx], ex_dates, nEx, date_today, Zeta1, Zeta2, alpha, gamma, rho);
        Zeta12[jEx - 1] = New_Zeta12_Func(
            ex_dates[jEx], ex_dates, nEx, date_today, Zeta1, Zeta12, alpha, gamma, rho);

        num_cpn = nPay - first_cpn_index[jEx] - 1;

        Rf  = SParms.long_swps_strikes[jEx];
        err = srt_f_autocal_swpts_greeks(
            jEx, num_cpn, Rf, long_swps_payts, Zeta1[jEx - 1], &diag_swaption_pv);
        if (err)
            return err;

        (*levenberg_value) += diag_swaption_pv;

        if (zeta1_param[1] == 0)
            param_shift = zero_shift;
        else
            param_shift = prop_shift * fabs(zeta1_param[1]);

        zeta1_ex += param_shift;
        Zeta1[jEx - 1] = zeta1_ex;

        Zeta2[jEx - 1] = New_Zeta2_Func(
            ex_dates[jEx], ex_dates, nEx, date_today, Zeta1, Zeta2, alpha, gamma, rho);
        Zeta12[jEx - 1] = New_Zeta12_Func(
            ex_dates[jEx], ex_dates, nEx, date_today, Zeta1, Zeta12, alpha, gamma, rho);

        err = srt_f_autocal_swpts_greeks(
            jEx, num_cpn, Rf, long_swps_payts, Zeta1[jEx - 1], &diag_swaption_pv);
        if (err)
            return err;

        levenberg_value_der[1] = 0.0;
        levenberg_value_der[1] += diag_swaption_pv;
        levenberg_value_der[1] -= (*levenberg_value);
        levenberg_value_der[1] /= param_shift;

        zeta1_ex -= param_shift;
        zeta1_param[1] = Zeta1[jEx - 1] = zeta1_ex;

        Zeta2[jEx - 1] = New_Zeta2_Func(
            ex_dates[jEx], ex_dates, nEx, date_today, Zeta1, Zeta2, alpha, gamma, rho);
        Zeta12[jEx - 1] = New_Zeta12_Func(
            ex_dates[jEx], ex_dates, nEx, date_today, Zeta1, Zeta12, alpha, gamma, rho);
    }
    else
        (*levenberg_value) += 10.;

    return err;
}
static Err srt_f_autocal_calboot(double* optim_params, LGMCalMeth cal_method)
{
    Err   err         = NULL;
    Date *start_dates = SParms.start_dates, *ex_dates = SParms.ex_dates,
         *caplet_end_dates = SParms.caplet_end_dates, last_swp_pay_date = SParms.last_swp_pay_date,
         date_today = SParms.date_today;
    long    i, nEx = SParms.nEx, nPay = SParms.nPay;
    double *H1 = SParms.H1, *H2 = SParms.H2, alpha = SParms.alpha, gamma = SParms.gamma,
           rho = SParms.rho, *Zeta1 = SParms.Zeta1, *Zeta2 = SParms.Zeta2, *Zeta12 = SParms.Zeta12,
           *df_start = SParms.df_start, *df_cap = SParms.df_cap, *cvg_cap = SParms.cvg_cap;
    double      H1tStart, H2tStart, H1tEnd, H2tEnd, zeta_den, chisq;
    double      nstop, zeta1_min, zeta1_max;
    long        niter, j;
    SRT_Boolean apply_criteria = SParms.apply_criteria;
    double      min_zeta_diff, zeta_bound_post_fact, *lgm_cap_vol = SParms.lgm_cap_vol;

    /* BUILD H1 */
    switch (cal_method)
    {
    case TenorAndDiag:
        for (i = 1; i <= nEx; i++)
        {
            if (i == nEx)
                H1[i - 1] = optim_params[i];
            else
                H1[i - 1] = H1[i] + optim_params[i];
        }
        break;
    default:
        for (i = 1; i <= nEx; i++)
            H1[i - 1] = (exp(-s_yr_to_start_ptr[i] * optim_params[1]) -
                         exp(-s_yr_to_last_pay_date * optim_params[1])) /
                        optim_params[1];

        for (i = 1; i <= nEx; i++)
            H2[i - 1] = H2_Func(start_dates[i], start_dates, date_today, nEx, H1, gamma);

        break;
    }
    /* FIND THE FIRST GUESS */

    if (SParms.use_start_point_for_zeta1s == SRT_NO)
    {
        H1tStart = H1[0];
        H1tEnd   = H1_Func(caplet_end_dates[1], start_dates, date_today, nEx, H1);
        H2tStart = H2[0];
        H2tEnd   = H2_Func(caplet_end_dates[1], start_dates, date_today, nEx, H1, gamma);

        zeta_den = (H1tEnd - H1tStart) * (H1tEnd - H1tStart) +
                   (alpha * alpha) / (2 * gamma * s_yr_to_ex_dates_ptr[1]) *
                       (exp(2 * gamma * s_yr_to_ex_dates_ptr[1]) - 1) * (H2tEnd - H2tStart) *
                       (H2tEnd - H2tStart) +
                   2 * (alpha * rho) / (gamma * s_yr_to_ex_dates_ptr[1]) *
                       (exp(gamma * s_yr_to_ex_dates_ptr[1]) - 1) * (H1tEnd - H1tStart) *
                       (H2tEnd - H2tStart);
    }

    niter = 2;
    /* BOOTSTRAP ON DIAG SWAPTIONS */
    for (i = 1; i <= nEx; i++)
    {
        SParms.ex_index = i;
        /* GET THE GUESS FOR LEVENBERG */
        if (i == 1)
        {
            zeta_bound_post_fact = ONE_HALF *
                                   (exp(2 * s_yr_to_ex_dates_ptr[1] * optim_params[1]) - 1) /
                                   optim_params[1];

            SParms.zeta1_min = zeta1_min = FIRST_EX_SIGMA_MIN_SQ * zeta_bound_post_fact;
            SParms.zeta1_max = zeta1_max = SIGMA_MAX_SQ * zeta_bound_post_fact;

            if (SParms.use_start_point_for_zeta1s == SRT_NO)
                s_zeta1_ptr[1] =
                    min(max(s_yr_to_ex_dates_ptr[1] * lgm_cap_vol[1] * lgm_cap_vol[1] / zeta_den,
                            zeta1_min),
                        zeta1_max);

            else
                s_zeta1_ptr[1] = min(max(SParms.start_zeta1s[0], zeta1_min), zeta1_max);
        }
        else
        {
            if (SParms.use_start_point_for_zeta1s == SRT_NO)
            {
                H1tStart = H1[i - 1];
                H1tEnd   = H1_Func(caplet_end_dates[i], start_dates, date_today, nEx, H1);
                H2tStart = H2[i - 1];
                H2tEnd   = H2_Func(caplet_end_dates[i], start_dates, date_today, nEx, H1, gamma);

                zeta_den = (H1tEnd - H1tStart) * (H1tEnd - H1tStart) +
                           (alpha * alpha) / (2 * gamma * s_yr_to_ex_dates_ptr[i]) *
                               (exp(2 * gamma * s_yr_to_ex_dates_ptr[i]) - 1) *
                               (H2tEnd - H2tStart) * (H2tEnd - H2tStart) +
                           2 * (alpha * rho) / (gamma * s_yr_to_ex_dates_ptr[i]) *
                               (exp(gamma * s_yr_to_ex_dates_ptr[i]) - 1) * (H1tEnd - H1tStart) *
                               (H2tEnd - H2tStart);
            }

            min_zeta_diff = ONE_HALF * SIGMA_MIN_SQ *
                            (exp(2 * s_yr_to_ex_dates_ptr[i] * optim_params[1]) -
                             exp(2 * s_yr_to_ex_dates_ptr[i - 1] * optim_params[1])) /
                            optim_params[1];
            SParms.zeta1_min = zeta1_min = Zeta1[i - 2] + min_zeta_diff;
            SParms.zeta1_max             = zeta1_max =
                ONE_HALF * SIGMA_MAX_SQ * (exp(2 * s_yr_to_ex_dates_ptr[i] * optim_params[1]) - 1) /
                optim_params[1];

            if (SParms.use_start_point_for_zeta1s == SRT_NO)
                s_zeta1_ptr[1] =
                    min(max(s_yr_to_ex_dates_ptr[i] * lgm_cap_vol[i] * lgm_cap_vol[i] / zeta_den,
                            zeta1_min),
                        zeta1_max);

            else
                s_zeta1_ptr[1] = min(max(SParms.start_zeta1s[i - 1], zeta1_min), zeta1_max);

        } /* END OF GETTING THE GUESS FOR LEVENBERG */

        s_levenberg_marquard_target_for_swpcalboot[1] = SParms.long_swps_pv[i];
        s_levenberg_marquard_weight_for_swpcalboot[1] = SParms.long_swps_vega[i];

        nstop = 0.;
        j     = 0;
        while ((nstop < 1.0) && (j < MAX_ITERATION))
        {
            j++;
            err = levenberg_marquardt(
                s_levenberg_marquard_data_for_swpcalboot,
                s_levenberg_marquard_target_for_swpcalboot,
                s_levenberg_marquard_weight_for_swpcalboot,
                1,
                s_zeta1_ptr,
                1,
                niter,
                &levenberg_func_for_tau,
                &chisq);
            if (err)
            {
                err = delete_levenberg_struct(nEx);
                if (err)
                    return err;
                return err;
            }

            nstop = 1.0;

        } /* END OF WHILE */

        if ((chisq > SWAPTION_SQ_VOL_TOL) && (apply_criteria == SRT_TRUE) && (0))
        {
            err = make_static_desallocations(nEx, nPay);
            if (err)
                return err;

            err = delete_levenberg_struct(nEx);
            if (err)
                return err;

            return serror("Can not match swaption price");
        }

        Zeta1[i - 1] = s_zeta1_ptr[1];
    }

    return err;
}

Err XLLGMZeta12Func(double t, char* UndName, double* Zeta)
{
    Err         err = NULL;
    TermStruct* ts;
    SrtTFTSMat  G;
    SrtUndPtr   sUndPtr;

    sUndPtr = lookup_und(UndName);
    err     = get_underlying_ts(sUndPtr, &ts);
    if (err)
        return err;

    err = get_2f_G_funcs(t, ts, &G);
    if (err)
        return err;

    *Zeta = G[0][1];

    return err;
}

Err XLLGMH1Func(double t, char* UndName, double* H)
{
    Err         err = NULL;
    TermStruct* ts;
    SrtUndPtr   Und;
    SrtTFTSVec  Psi;
    long        mdl_dim;

    Und = lookup_und(UndName);
    if (!Und)
        return serror("Underlying %s not defined", UndName);

    err = get_underlying_ts(Und, &ts);
    if (err)
        return err;

    err = get_underlying_mdldim(Und, &mdl_dim);
    if (err)
        return err;

    if (mdl_dim == TWO_FAC)
    {
        err = get_2f_Psi_funcs(t, ts, &Psi);
        if (err)
            return err;

        (*H) = Psi[0];
    }
    else
    {
        (*H) = Psi_func(t, ts);
    }

    return err;
}

Err XLLGMH2Func(double t, char* UndName, double* H)
{
    Err         err = NULL;
    TermStruct* ts;
    SrtUndPtr   Und;
    SrtTFTSVec  Psi;

    Und = lookup_und(UndName);
    if (!Und)
        return serror("Underlying %s not defined", UndName);

    err = get_underlying_ts(Und, &ts);
    if (err)
        return err;

    err = get_2f_Psi_funcs(t, ts, &Psi);
    if (err)
        return err;

    (*H) = Psi[1];

    return err;
}

static Err make_swaptions_payts(
    long     nEx,
    long     nPay,
    long*    first_cpn_index,
    double*  cvg_first,
    double*  cvg_pay,
    double*  df_pay,
    double*  K,
    double** SwapPayts)
{
    Err  err = NULL;
    long j, k, num_cpn;

    for (j = 1; j <= nEx; j++) /* Loop on all the swaption, i.e loop on all the exercise dates */
    {
        num_cpn = nPay - first_cpn_index[j] - 1;

        for (k = 0; k <= num_cpn; k++) /* Loop on Cpns */
        {
            if ((k == 0) && (num_cpn > 0))
                SwapPayts[j][k] = df_pay[first_cpn_index[j]] * cvg_first[j] * K[j];
            else if ((k > 0) && (k < num_cpn))
                SwapPayts[j][k] =
                    df_pay[first_cpn_index[j] + k] * cvg_pay[first_cpn_index[j] + k - 1] * K[j];
            else if (k == num_cpn)
            {
                if (num_cpn > 0)
                    SwapPayts[j][k] = df_pay[first_cpn_index[j] + k] *
                                      (1 + K[j] * cvg_pay[first_cpn_index[j] + k - 1]);
                else if (num_cpn == 0)
                    SwapPayts[j][k] = df_pay[first_cpn_index[j] + k] * (1 + K[j] * cvg_first[j]);
            }
        } /* End of loop on k */

    } /* End of loop on j */

    return err;
}

static Err make_short_swaptions_payts(
    long     nEx,
    long*    num_cpn,
    long*    first_cpn_index,
    double*  cvg_first,
    double*  cvg_pay,
    double*  df_pay,
    double*  K,
    double** SwapPayts)
{
    Err  err = NULL;
    long j, k;

    for (j = 1; j <= nEx; j++) /* Loop on all the swaption, i.e loop on all the exercise dates */
    {
        for (k = 0; k <= num_cpn[j]; k++) /* Loop on Cpns */
        {
            if ((k == 0) && (num_cpn[j] > 0))
                SwapPayts[j][k] = df_pay[first_cpn_index[j]] * cvg_first[j] * K[j];
            else if ((k > 0) && (k < num_cpn[j]))
                SwapPayts[j][k] =
                    df_pay[first_cpn_index[j] + k] * cvg_pay[first_cpn_index[j] + k - 1] * K[j];
            else if (k == num_cpn[j])
            {
                if (num_cpn[j] > 0)
                    SwapPayts[j][k] = df_pay[first_cpn_index[j] + k] *
                                      (1 + K[j] * cvg_pay[first_cpn_index[j] + k - 1]);
                else if (num_cpn[j] == 0)
                    SwapPayts[j][k] = df_pay[first_cpn_index[j] + k] * (1 + K[j] * cvg_first[j]);
            }
        } /* End of loop on k */

    } /* End of loop on j */

    return err;
}

static Err make_static_struct_allocations(
    long      nEx,
    long      nPay,
    Date**    caplet_end_dates,
    long**    first_cpn_index,
    double**  df_cap,
    double**  df_start,
    double**  df_pay,
    double**  cvg_pay,
    double**  cvg_cap,
    double**  cvg_first,
    double**  long_swps_strikes,
    double*** long_swps_payts,
    double**  long_swps_fwd,
    double**  long_swps_cev,
    double**  long_swps_beta,
    double**  long_swps_pv,
    double**  long_swps_vega,
    double**  long_swps_stdev,
    double**  short_swps_strikes,
    double*** short_swps_payts,
    double**  short_swps_fwd,
    double**  short_swps_cev,
    double**  short_swps_beta,
    double**  short_swps_pv,
    double**  short_swps_vega,
    long**    num_short_swp_cpn,
    double**  caplets_pv,
    double**  caplets_vega,
    double**  caplets_strikes,
    double**  lgm_cap_vol)
{
    Err err = NULL;

    (*caplet_end_dates)  = lngvector(1, nEx);
    (*first_cpn_index)   = lngvector(1, nEx);
    (*df_cap)            = dvector(1, nEx);
    (*df_start)          = dvector(1, nEx);
    (*df_pay)            = dvector(0, nPay - 1);
    (*cvg_pay)           = dvector(0, nPay - 1);
    (*cvg_cap)           = dvector(1, nEx);
    (*cvg_first)         = dvector(1, nEx);
    (*long_swps_strikes) = dvector(1, nEx);
    (*long_swps_payts)   = dmatrix(1, nEx, 0, nPay - 1);
    (*long_swps_fwd)     = dvector(1, nEx);
    (*long_swps_cev)     = dvector(1, nEx);
    (*long_swps_beta)    = dvector(1, nEx);
    (*long_swps_pv)      = dvector(1, nEx);
    (*long_swps_vega)    = dvector(1, nEx);
    (*long_swps_stdev)   = dvector(1, nEx);

    (*short_swps_strikes) = dvector(1, nEx);
    (*short_swps_payts)   = dmatrix(1, nEx, 0, nPay - 1);
    (*short_swps_fwd)     = dvector(1, nEx);
    (*short_swps_cev)     = dvector(1, nEx);
    (*short_swps_beta)    = dvector(1, nEx);
    (*short_swps_pv)      = dvector(1, nEx);
    (*short_swps_vega)    = dvector(1, nEx);
    (*num_short_swp_cpn)  = lngvector(1, nEx);
    (*caplets_pv)         = dvector(1, nEx);
    (*caplets_vega)       = dvector(1, nEx);
    (*caplets_strikes)    = dvector(1, nEx);

    (*lgm_cap_vol) = dvector(1, nEx);

    return NULL;
}

static Err make_static_struct_desallocations(long nEx, long nPay)
{
    Err err = NULL;

    if (SParms.caplet_end_dates)
        free_lngvector(SParms.caplet_end_dates, 1, nEx);
    SParms.caplet_end_dates = NULL;
    if (SParms.first_cpn_index)
        free_lngvector(SParms.first_cpn_index, 1, nEx);
    SParms.first_cpn_index = NULL;
    if (SParms.df_cap)
        free_dvector(SParms.df_cap, 1, nEx);
    SParms.df_cap = NULL;
    if (SParms.df_start)
        free_dvector(SParms.df_start, 1, nEx);
    SParms.df_start = NULL;
    if (SParms.df_pay)
        free_dvector(SParms.df_pay, 0, nPay - 1);
    SParms.df_pay = NULL;
    if (SParms.cvg_pay)
        free_dvector(SParms.cvg_pay, 0, nPay - 1);
    SParms.cvg_pay = NULL;
    if (SParms.cvg_cap)
        free_dvector(SParms.cvg_cap, 1, nEx);
    SParms.cvg_cap = NULL;
    if (SParms.cvg_first)
        free_dvector(SParms.cvg_first, 1, nEx);
    SParms.cvg_first = NULL;

    if (SParms.long_swps_strikes)
        free_dvector(SParms.long_swps_strikes, 1, nEx);
    SParms.long_swps_strikes = NULL;
    if (SParms.long_swps_payts)
        free_dmatrix(SParms.long_swps_payts, 1, nEx, 0, nPay - 1);
    SParms.long_swps_payts = NULL;
    if (SParms.long_swps_fwd)
        free_dvector(SParms.long_swps_fwd, 1, nEx);
    SParms.long_swps_fwd = NULL;
    if (SParms.long_swps_cev)
        free_dvector(SParms.long_swps_cev, 1, nEx);
    SParms.long_swps_cev = NULL;
    if (SParms.long_swps_beta)
        free_dvector(SParms.long_swps_beta, 1, nEx);
    SParms.long_swps_beta = NULL;
    if (SParms.long_swps_pv)
        free_dvector(SParms.long_swps_pv, 1, nEx);
    SParms.long_swps_pv = NULL;
    if (SParms.long_swps_vega)
        free_dvector(SParms.long_swps_vega, 1, nEx);
    SParms.long_swps_vega = NULL;
    if (SParms.long_swps_stdev)
        free_dvector(SParms.long_swps_stdev, 1, nEx);
    SParms.long_swps_stdev = NULL;

    if (SParms.short_swps_strikes)
        free_dvector(SParms.short_swps_strikes, 1, nEx);
    SParms.short_swps_strikes = NULL;
    if (SParms.short_swps_payts)
        free_dmatrix(SParms.short_swps_payts, 1, nEx, 0, nPay - 1);
    SParms.short_swps_payts = NULL;
    if (SParms.short_swps_fwd)
        free_dvector(SParms.short_swps_fwd, 1, nEx);
    SParms.short_swps_fwd = NULL;

    if (SParms.short_swps_cev)
        free_dvector(SParms.short_swps_cev, 1, nEx);
    SParms.short_swps_cev = NULL;
    if (SParms.short_swps_beta)
        free_dvector(SParms.short_swps_beta, 1, nEx);
    SParms.short_swps_beta = NULL;
    if (SParms.short_swps_pv)
        free_dvector(SParms.short_swps_pv, 1, nEx);
    SParms.short_swps_pv = NULL;
    if (SParms.short_swps_vega)
        free_dvector(SParms.short_swps_vega, 1, nEx);
    SParms.short_swps_vega = NULL;
    if (SParms.caplets_pv)
        free_dvector(SParms.caplets_pv, 1, nEx);
    SParms.caplets_pv = NULL;
    if (SParms.caplets_vega)
        free_dvector(SParms.caplets_vega, 1, nEx);
    SParms.caplets_vega = NULL;
    if (SParms.caplets_strikes)
        free_dvector(SParms.caplets_strikes, 1, nEx);
    SParms.caplets_strikes = NULL;

    if (SParms.num_short_swp_cpn)
        free_lngvector(SParms.num_short_swp_cpn, 1, nEx);
    SParms.num_short_swp_cpn = NULL;
    if (SParms.lgm_cap_vol)
        free_dvector(SParms.lgm_cap_vol, 1, nEx);
    SParms.lgm_cap_vol = NULL;

    return NULL;
}

Err srt_f_get_autocal_ref_instrs(
    long  nEx,
    long  nPay,
    Date* start_dates,
    Date* swp_pay_dates,
    char* ycname,
    Err (*srt_f_get_vol)(Date, Date, double, SRT_Boolean, double*),
    Err (*srt_f_get_beta)(Date, Date, double*),
    double*             LongKs,
    LGMRMeth            strike_method,
    LGMCalMeth          calib_method,
    long*               likely_most_expensive_index,
    SrtLgmRefSwptnData* lgm_ref_swp_data)

{
    SrtCurvePtr      yld_crv;
    String           ccy_str;
    SrtCcyParam*     ccy_param  = NULL;
    SrtCompounding   float_freq = SRT_SEMIANNUAL;
    SrtBasisCode     fixed_basis, float_basis;
    SrtBusDayConv    float_bus_day_conv, fixed_bus_day_conv;
    SrtDiffusionType vol_type;
    SwapDP*          float_leg_DP;

    Date date_today, *ex_dates, *caplet_end_dates, cap_theo_end_date, cap_real_end_date,
        *float_pay_dates, *float_fixing_dates, *float_start_dates, *float_end_dates;
    Err  err = NULL;
    long i, j, k, num_cpn, num_cap_fp, last_cap_ex_index, spot_lag, most_exp_index,
        *first_cpn_index, *num_short_swp_cpn;
    int    num_float_dates, num_float_pay_dates;
    double cev_vol, swp_level, swp_rate, fra_rate, yr_to_exp, atm_vol, norm_vol, beta, caplet_pv,
        caplet_vega, cap_pv, cap_vega, hybrid_swpts_pv, hybrid_swpts_vega, yr_from_prev_ex_date,
        nbr_of_stdev_from_prev_strike, sgn, loc_max_stdev, *float_coverages,
        *float_start_dates_discount_factors, *float_end_dates_discount_factors;

    double *df_cap, *df_start, *df_pay, *cvg_pay, *cvg_cap, *cvg_first, *long_swps_strikes,
        **long_swps_payts, *long_swps_fwd, *long_swps_cev, *long_swps_beta, *long_swps_pv,
        *long_swps_vega, *long_swps_stdev, *short_swps_strikes, **short_swps_payts, *lgm_cap_vol,
        *short_swps_fwd, *short_swps_cev, *short_swps_beta, *short_swps_pv, *short_swps_vega,
        *caplets_pv, *caplets_vega, *caplets_strikes;

    /* Make the allocations  (the static first) */
    err = make_static_struct_allocations(
        nEx,
        nPay,
        &caplet_end_dates,
        &first_cpn_index,
        &df_cap,
        &df_start,
        &df_pay,
        &cvg_pay,
        &cvg_cap,
        &cvg_first,
        &long_swps_strikes,
        &long_swps_payts,
        &long_swps_fwd,
        &long_swps_cev,
        &long_swps_beta,
        &long_swps_pv,
        &long_swps_vega,
        &long_swps_stdev,
        &short_swps_strikes,
        &short_swps_payts,
        &short_swps_fwd,
        &short_swps_cev,
        &short_swps_beta,
        &short_swps_pv,
        &short_swps_vega,
        &num_short_swp_cpn,
        &caplets_pv,
        &caplets_vega,
        &caplets_strikes,
        &lgm_cap_vol);
    if (err)
        return err;

    ex_dates = lngvector(1, nEx);
    if (lgm_ref_swp_data != NULL)
    {
        lgm_ref_swp_data->NrefLongSwptn = nEx;
        lgm_ref_swp_data->refLongSwptnArr =
            (SrtLgmRefSwptn*)srt_calloc(lgm_ref_swp_data->NrefLongSwptn, sizeof(SrtLgmRefSwptn));

        if (lgm_ref_swp_data->refLongSwptnArr == NULL)
        {
            lgm_ref_swp_data->NrefLongSwptn = 0;
        }
    }

    /* Get some market informations from the yc name */
    yld_crv    = lookup_curve(ycname);
    date_today = get_clcndate_from_yldcrv(yld_crv);
    ccy_param  = get_ccyparam_from_yldcrv(yld_crv);

    if (!ccy_param)
    {
        ccy_str = get_curve_ccy(yld_crv);
        err     = swp_f_get_CcyParam_from_CcyStr(ccy_str, &ccy_param);
        if (err)
            return (err);
    }

    fixed_basis        = ccy_param->swap_basis_code;
    fixed_bus_day_conv = ccy_param->swap_bus_day_conv;

    float_basis        = ccy_param->cash_basis_code;
    float_bus_day_conv = ccy_param->cash_bus_day_conv;

    spot_lag = get_spotlag_from_curve(yld_crv);

    /* Compute the exercise dates from the start dates - ex dates is spot_lag bd before the start
     * dates  */
    for (i = 1; i <= nEx; i++)
        ex_dates[i] = add_unit(start_dates[i], -(spot_lag), SRT_BDAY, SUCCEEDING);

    /* Compute the discount factor at the midat payts dates */
    for (i = 0; i < nPay; i++)
    {
        df_pay[i] = swp_f_df(date_today, swp_pay_dates[i], ycname);
        if (df_pay[i] == SRT_DF_ERROR)
            return ("No Discount Factor");
    }
    /* Compute the coverages between two payts dates */
    for (i = 0; i < (nPay - 1); i++)
        cvg_pay[i] = coverage(swp_pay_dates[i], swp_pay_dates[i + 1], fixed_basis);

    /* Compute the caplets end dates (caplet_end_dates) ,the caplets coverages (cvg_cap), the long
      swaptions first coverages (cvg_first),
      the discount factors at the start dates (df_start) and the discount factors at the caplets end
      dates (df_cap) */
    for (j = 1; j <= nEx; j++)
    {
        /* Find the first pay date strictly after start_dates[j] */
        i = BasicDichotomie(swp_pay_dates, 0, (nPay - 1), start_dates[j]);

        if (swp_pay_dates[i] == start_dates[j])
            first_cpn_index[j] = i + 1;
        else
            first_cpn_index[j] = i;

        caplet_end_dates[j] =
            add_unit(start_dates[j], 12 / (int)float_freq, SRT_MONTH, float_bus_day_conv);

        cvg_cap[j] = coverage(start_dates[j], caplet_end_dates[j], float_basis);

        cvg_first[j] = coverage(start_dates[j], swp_pay_dates[first_cpn_index[j]], fixed_basis);

        df_start[j] = swp_f_df(date_today, start_dates[j], ycname);
        if (df_start[j] == SRT_DF_ERROR)
            return ("No Discount Factor ");

        df_cap[j] = swp_f_df(date_today, caplet_end_dates[j], ycname);
        if (df_cap[j] == SRT_DF_ERROR)
            return ("No Discount Factor");
    }

    /* For midat with monthly exercise dates, reduce the max stdev to 1.5 */
    loc_max_stdev = MAX_STDEV;
    if (nEx > 1)
    {
        if ((ex_dates[nEx] - ex_dates[nEx - 1]) < DAYS_IN_QUARTER)
            loc_max_stdev = 1.5;
        else
            loc_max_stdev = MAX_STDEV;
    }

    /* Get the market instrument prices */
    for (i = 1; i <= nEx; i++)
    {
        /* Get the long swaptions strikes and cut them to loc_max_stdev */

        long_swps_strikes[i] = LongKs[i]; /* initialisation */
        /* Compute the cash at money swap rate (long_swps_fwd)*/
        swp_level = 0.;
        num_cpn =
            nPay - first_cpn_index[i] - 1; /* num_cpn is the number of coupon after the first one */

        swp_level = 0.0;
        for (k = 0; k <= num_cpn; k++) /* k = 0 corresponds to the first coupon */
        {
            if ((k == 0) && (num_cpn > 0))
                swp_level += cvg_first[i] * df_pay[first_cpn_index[i]];
            else if ((k > 0) && (k < num_cpn))
                swp_level += cvg_pay[first_cpn_index[i] + k - 1] * df_pay[first_cpn_index[i] + k];
            else if (k == num_cpn)
            {
                if (num_cpn > 0)
                    swp_level +=
                        cvg_pay[first_cpn_index[i] + k - 1] * df_pay[first_cpn_index[i] + k];
                else if (num_cpn == 0)
                    swp_level += cvg_first[i] * df_pay[first_cpn_index[i] + k];
            }
        }

        long_swps_fwd[i] = swp_rate =
            (df_start[i] - df_pay[first_cpn_index[i] + num_cpn]) / swp_level;
        yr_to_exp = (double)(ex_dates[i] - date_today) * YEARS_IN_DAY;

        /* COMPUTE THE STANDARD DEVIATION AND CUT IT IF NECESSARY */
        err = srt_f_get_vol(start_dates[i], swp_pay_dates[nPay - 1], swp_rate, SRT_FALSE, &atm_vol);
        if (err)
            return err;

        err = srt_f_get_beta(start_dates[i], swp_pay_dates[nPay - 1], &beta);
        if (err)
            return err;
        long_swps_beta[i] = beta;

        if (beta == 1)
            vol_type = SRT_LOGNORMAL;
        else
            vol_type = SRT_NORMAL;

        err = srt_f_optsarbvol(
            swp_rate, swp_rate, yr_to_exp, atm_vol, 0.0, 0.0, 0.0, vol_type, SRT_NORMAL, &norm_vol);
        if (err)
            return err;

        long_swps_stdev[i] = (long_swps_strikes[i] - swp_rate) / (norm_vol * sqrt(yr_to_exp));

        if (fabs(long_swps_stdev[i]) > loc_max_stdev)
        {
            if (long_swps_stdev[i] > 0)
                long_swps_stdev[i] = min(loc_max_stdev, long_swps_stdev[i]);
            else
                long_swps_stdev[i] = max(-loc_max_stdev, long_swps_stdev[i]);
        }

        long_swps_strikes[i] = long_swps_stdev[i] * norm_vol * sqrt(yr_to_exp) +
                               swp_rate; /* END OF COMPUTATION OF THE STDEV */

        /* VERIFY IF THE STRIKE IS AT LESS THAN MAX_NR_STDEV_FROM_PREV_STRIKE FROM THE PREVIOUS ONE
        - FOR STEP UP AND ZERO COUPON MIDAT */
        if (i > 1)
        {
            yr_from_prev_ex_date          = YEARS_IN_DAY * (ex_dates[i] - ex_dates[i - 1]);
            nbr_of_stdev_from_prev_strike = (long_swps_strikes[i] - long_swps_strikes[i - 1]) /
                                            (norm_vol * sqrt(yr_from_prev_ex_date));

            if (fabs(nbr_of_stdev_from_prev_strike) > MAX_NR_STDEV_FROM_PREV_STRIKE)
            {
                if (long_swps_strikes[i] > long_swps_strikes[i - 1])
                    sgn = 1.0;
                else
                    sgn = -1.0;

                long_swps_strikes[i] =
                    sgn * MAX_NR_STDEV_FROM_PREV_STRIKE * norm_vol * sqrt(yr_from_prev_ex_date) +
                    long_swps_strikes[i - 1];
            }
        }

        /* END OF COMPUTATION OF THE STRIKE */

        /* GET THE PV AND THE VEGA OF THE LONG SWAPTIONS */

        /* GET VOL, RECORD VOL & GET CEV EXPONENT */
        err = srt_f_get_vol(
            start_dates[i], swp_pay_dates[nPay - 1], long_swps_strikes[i], SRT_FALSE, &cev_vol);
        if (err)
            return err;
        long_swps_cev[i] = cev_vol;

        /*COMPUTATION OF THE REFS SWAPTION VEGA */

        long_swps_pv[i] = 0.;
        if (beta == 1.0) /* THEN THE DIFFUSION IS LOGNORMAL */
        {
            vol_type = SRT_LOGNORMAL;
            /* COMPUTE THE PREMIUM */
            long_swps_pv[i] = srt_f_optblksch(
                swp_rate, long_swps_strikes[i], cev_vol, yr_to_exp, 1.0, SRT_PUT, SRT_PREMIUM);
            long_swps_pv[i] *= swp_level;

            if (long_swps_pv[i] < 0.)
                return ("Fatal Error while Computing Swaption Premium ");

            /* COMPUTE THE VEGA */
            long_swps_vega[i] = 0.;
            long_swps_vega[i] = srt_f_optblksch(
                swp_rate,
                long_swps_strikes[i],
                cev_vol + CEV_VOL_SHIFT,
                yr_to_exp,
                1.0,
                SRT_PUT,
                SRT_PREMIUM);
            long_swps_vega[i] *= swp_level;

            long_swps_vega[i] -= long_swps_pv[i];
            long_swps_vega[i] /= CEV_VOL_SHIFT;
        }
        else if (beta == 0.0)
        {
            vol_type = SRT_NORMAL;

            long_swps_pv[i] = srt_f_optblknrm(
                swp_rate, long_swps_strikes[i], cev_vol, yr_to_exp, 1.0, SRT_PUT, SRT_PREMIUM);
            long_swps_pv[i] *= swp_level;

            if (long_swps_pv[i] < 0.)
                return ("Fatal Error while Computing Swaption Premium ");

            /* COMPUTE THE VEGA */
            long_swps_vega[i] = 0.;
            long_swps_vega[i] = srt_f_optblknrm(
                swp_rate,
                long_swps_strikes[i],
                cev_vol + CEV_VOL_SHIFT,
                yr_to_exp,
                1.0,
                SRT_PUT,
                SRT_PREMIUM);
            long_swps_vega[i] *= swp_level;

            long_swps_vega[i] -= long_swps_pv[i];
            long_swps_vega[i] /= CEV_VOL_SHIFT;
        }
        else
            return serror("Expect a beta of 0 (NORMAL) or 1.0 (LOGNORMAL)");

        /* store the ref swpts details */
        if (lgm_ref_swp_data != NULL)
        {
            lgm_ref_swp_data->refLongSwptnArr[i - 1].bgnDt  = start_dates[i];
            lgm_ref_swp_data->refLongSwptnArr[i - 1].endDt  = swp_pay_dates[nPay - 1];
            lgm_ref_swp_data->refLongSwptnArr[i - 1].strike = long_swps_strikes[i];
            lgm_ref_swp_data->refLongSwptnArr[i - 1].bsPv   = long_swps_pv[i];
            lgm_ref_swp_data->refLongSwptnArr[i - 1].bsVol  = cev_vol;
        }

        /* COMPUTE THE SRIKE OF THE CAPLETS */
        fra_rate = (df_start[i] / df_cap[i] - 1) / cvg_cap[i];

        err = srt_f_get_vol(start_dates[i], caplet_end_dates[i], fra_rate, SRT_TRUE, &atm_vol);
        if (err)
            return err;

        err = srt_f_optsarbvol(
            fra_rate, fra_rate, yr_to_exp, atm_vol, 0.0, 0.0, 0.0, vol_type, SRT_NORMAL, &norm_vol);
        if (err)
            return err;

        if (strike_method == EMK1)
            short_swps_strikes[i] = fra_rate;
        else
            short_swps_strikes[i] = fra_rate + long_swps_stdev[i] * norm_vol * sqrt(yr_to_exp);

        /* VERIFY IF THE STRIKE IS AT LESS THAN MAX_NR_STDEV_FROM_PREV_STRIKE FROM THE PREVIOUS ONE
         */
        if (i > 1)
        {
            yr_from_prev_ex_date          = YEARS_IN_DAY * (ex_dates[i] - ex_dates[i - 1]);
            nbr_of_stdev_from_prev_strike = (short_swps_strikes[i] - short_swps_strikes[i - 1]) /
                                            (norm_vol * sqrt(yr_from_prev_ex_date));

            if (fabs(nbr_of_stdev_from_prev_strike) > MAX_NR_STDEV_FROM_PREV_STRIKE)
            {
                if (short_swps_strikes[i] > short_swps_strikes[i - 1])
                    sgn = 1.0;
                else
                    sgn = -1.0;

                short_swps_strikes[i] =
                    sgn * MAX_NR_STDEV_FROM_PREV_STRIKE * norm_vol * sqrt(yr_from_prev_ex_date) +
                    short_swps_strikes[i - 1];
            }
        }
        /* GET THE VOLATILITY AT THIS STRIKE AND THE CAPLET PRICE */
        err = srt_f_get_vol(
            start_dates[i], caplet_end_dates[i], short_swps_strikes[i], SRT_TRUE, &cev_vol);
        if (err)
            return err;

        if (beta == 1.0)
        {
            short_swps_pv[i] = 0.;
            short_swps_pv[i] = cvg_cap[i] * srt_f_optblksch(
                                                fra_rate,
                                                short_swps_strikes[i],
                                                cev_vol,
                                                yr_to_exp,
                                                df_cap[i],
                                                SRT_PUT,
                                                SRT_PREMIUM);

            short_swps_vega[i] = 0.;
            short_swps_vega[i] = cvg_cap[i] * srt_f_optblksch(
                                                  fra_rate,
                                                  short_swps_strikes[i],
                                                  cev_vol + CEV_VOL_SHIFT,
                                                  yr_to_exp,
                                                  df_cap[i],
                                                  SRT_PUT,
                                                  SRT_PREMIUM);
            short_swps_vega[i] -= short_swps_pv[i];
            short_swps_vega[i] /= CEV_VOL_SHIFT;
        }
        else if (beta == 0.0)
        {
            short_swps_pv[i] = 0.;
            short_swps_pv[i] =
                cvg_cap[i] * df_cap[i] *
                srt_f_optblknrm(
                    fra_rate, short_swps_strikes[i], cev_vol, yr_to_exp, 1.0, SRT_PUT, SRT_PREMIUM);

            short_swps_vega[i] = 0.;
            short_swps_vega[i] = cvg_cap[i] * df_cap[i] *
                                 srt_f_optblknrm(
                                     fra_rate,
                                     short_swps_strikes[i],
                                     cev_vol + CEV_VOL_SHIFT,
                                     yr_to_exp,
                                     df_cap[i],
                                     SRT_PUT,
                                     SRT_PREMIUM);
            short_swps_vega[i] -= short_swps_pv[i];
            short_swps_vega[i] /= CEV_VOL_SHIFT;
        }
        else
            return serror("Expect a beta of 0 (NORMAL) or 1.0 (LOGNORMAL)");
    } /* END OF LOOP ON i FOR CAP */

    /* GET THE LONG SWAPTIONS PAYTS */
    err = make_swaptions_payts(
        nEx, nPay, first_cpn_index, cvg_first, cvg_pay, df_pay, long_swps_strikes, long_swps_payts);
    if (err)
        return err;

    if (calib_method == FullCapAndDiag)
    {
        /* GET THE FULL CAP DETAILS - GET THE SCHEDULE FIRST */
        float_leg_DP = (SwapDP*)malloc(sizeof(SwapDP));

        last_cap_ex_index = nEx;
        cap_real_end_date = add_unit(
            start_dates[last_cap_ex_index], 12 / (int)float_freq, SRT_MONTH, MODIFIED_SUCCEEDING);
        if (nEx > 1)
        {
            while (cap_real_end_date > swp_pay_dates[nPay - 1])
            {
                last_cap_ex_index -= 1;
                cap_real_end_date = add_unit(
                    start_dates[last_cap_ex_index],
                    12 / (int)float_freq,
                    SRT_MONTH,
                    MODIFIED_SUCCEEDING);
            }
        }

        cap_theo_end_date = add_unit(
            start_dates[last_cap_ex_index], 12 / (int)float_freq, SRT_MONTH, NO_BUSDAY_CONVENTION);

        err = swp_f_setSwapDP(
            start_dates[1], cap_theo_end_date, float_freq, float_basis, float_leg_DP);
        if (err)
            return err;

        err = swp_f_make_FloatLegDatesAndCoverages(
            float_leg_DP,
            date_today,
            &float_pay_dates,
            &num_float_pay_dates,
            &float_fixing_dates,
            &float_start_dates,
            &float_end_dates,
            &float_coverages,
            &num_float_dates);
        if (err)
            return err;

        if (float_start_dates[0] == float_end_dates[0])
        {
            num_cap_fp = num_float_dates;
            while (add_unit(
                       start_dates[1],
                       12 * num_cap_fp / (int)float_freq,
                       SRT_MONTH,
                       NO_BUSDAY_CONVENTION) > start_dates[last_cap_ex_index])
            {
                num_cap_fp -= 1;
            }

            num_cap_fp += 1;

            cap_theo_end_date = add_unit(
                start_dates[1], 12 * num_cap_fp / (int)float_freq, SRT_MONTH, NO_BUSDAY_CONVENTION);

            err = swp_f_setSwapDP(
                start_dates[1], cap_theo_end_date, float_freq, float_basis, float_leg_DP);
            if (err)
                return err;

            err = swp_f_make_FloatLegDatesAndCoverages(
                float_leg_DP,
                date_today,
                &float_pay_dates,
                &num_float_pay_dates,
                &float_fixing_dates,
                &float_start_dates,
                &float_end_dates,
                &float_coverages,
                &num_float_dates);
            if (err)
                return err;
        }

        float_start_dates_discount_factors = dvector(0, num_float_pay_dates - 1);
        float_end_dates_discount_factors   = dvector(0, num_float_pay_dates - 1);

        for (i = 1; i <= num_float_dates; i++)
        {
            float_start_dates_discount_factors[i - 1] =
                swp_f_df(date_today, float_start_dates[i - 1], ycname);
            if (float_start_dates_discount_factors[i - 1] == SRT_DF_ERROR)
                return ("No Discount Factor");

            float_end_dates_discount_factors[i - 1] =
                swp_f_df(date_today, float_end_dates[i - 1], ycname);
            if (float_end_dates_discount_factors[i - 1] == SRT_DF_ERROR)
                return ("No Discount Factor");

            float_fixing_dates[i - 1] =
                add_unit(float_start_dates[i - 1], -spot_lag, SRT_BDAY, SUCCEEDING);
        }

        if (lgm_ref_swp_data != NULL)
        {
            lgm_ref_swp_data->NrefCapSwptn = num_float_dates;
            lgm_ref_swp_data->refCapSwptnArr =
                (SrtLgmRefSwptn*)srt_calloc(lgm_ref_swp_data->NrefCapSwptn, sizeof(SrtLgmRefSwptn));

            if (lgm_ref_swp_data->refCapSwptnArr == NULL)
            {
                lgm_ref_swp_data->NrefCapSwptn = 0;
            }
        }

        /* GET THE FULL CAP PRICE */
        cap_pv   = 0.;
        cap_vega = 0.;

        for (i = 1; i <= num_float_dates; i++)
        {
            j = BasicDichotomie(start_dates, 1, nEx, float_start_dates[i - 1]);

            err = srt_f_get_vol(
                float_start_dates[i - 1],
                float_end_dates[i - 1],
                short_swps_strikes[j],
                SRT_TRUE,
                &cev_vol);
            if (err)
                return err;

            yr_to_exp = (double)(float_fixing_dates[i - 1] - date_today) * YEARS_IN_DAY;
            fra_rate  = (float_start_dates_discount_factors[i - 1] /
                            float_end_dates_discount_factors[i - 1] -
                        1) /
                       float_coverages[i - 1];

            caplet_pv = 0.;
            caplet_pv = LGMCEVCapletPrice(
                date_today,
                float_fixing_dates[i - 1],
                SRT_RECEIVER,
                short_swps_strikes[j],
                float_coverages[i - 1],
                float_start_dates_discount_factors[i - 1],
                float_end_dates_discount_factors[i - 1],
                cev_vol,
                beta);

            cap_pv += caplet_pv;

            caplet_vega = 0.;
            caplet_vega = LGMCEVCapletPrice(
                date_today,
                float_fixing_dates[i - 1],
                SRT_RECEIVER,
                short_swps_strikes[j],
                float_coverages[i - 1],
                float_start_dates_discount_factors[i - 1],
                float_end_dates_discount_factors[i - 1],
                cev_vol + CEV_VOL_SHIFT,
                beta);

            caplet_vega -= caplet_pv;
            caplet_vega /= CEV_VOL_SHIFT;

            cap_vega += caplet_vega;

            /* store the ref caplet details */
            if (lgm_ref_swp_data != NULL)
            {
                lgm_ref_swp_data->refCapSwptnArr[i - 1].bgnDt  = float_start_dates[i - 1];
                lgm_ref_swp_data->refCapSwptnArr[i - 1].endDt  = float_end_dates[i - 1];
                lgm_ref_swp_data->refCapSwptnArr[i - 1].strike = short_swps_strikes[j];
                lgm_ref_swp_data->refCapSwptnArr[i - 1].bsPv   = caplet_pv;
                lgm_ref_swp_data->refCapSwptnArr[i - 1].bsVol  = cev_vol;
            }
        }
    }
    else if (calib_method == HybridShortAndDiag)
    {
        /* allocate the structure for the storage of ref swpts */
        if (lgm_ref_swp_data != NULL)
        {
            lgm_ref_swp_data->NrefShortSwptn   = nEx;
            lgm_ref_swp_data->refShortSwptnArr = (SrtLgmRefSwptn*)srt_calloc(
                lgm_ref_swp_data->NrefShortSwptn, sizeof(SrtLgmRefSwptn));

            if (lgm_ref_swp_data->refShortSwptnArr == NULL)
            {
                lgm_ref_swp_data->NrefShortSwptn = 0;
            }
        }

        /* FIRST COPY CAPLETS PRICES AND VEGAS */

        for (i = 1; i <= nEx; i++)
        {
            caplets_pv[i]      = short_swps_pv[i];
            caplets_vega[i]    = short_swps_vega[i];
            caplets_strikes[i] = short_swps_strikes[i];
        }

        for (j = 1; j < nEx; j++)
        {
            most_exp_index = j + likely_most_expensive_index[j];
            if (most_exp_index <= nEx)
            {
                i = BasicDichotomie(
                    swp_pay_dates,
                    0,
                    (nPay - 1),
                    start_dates[most_exp_index]); /* start_dates[most_exp_index] IS THE REAL END
                                                     DATE OF THE SHORT SWAPTIONS */
                num_short_swp_cpn[j] = i - first_cpn_index[j]; /* num_short_swp_cpn IS THE NUMBER OF
                                                                  CPN AFTER THE FIRST ONE */
            }

            swp_level = 0.0;
            for (k = 0; k <= num_short_swp_cpn[j]; k++) /* k = 0 REPRESENTS THE FIRST COUPON */
            {
                if ((k == 0) && (num_short_swp_cpn[j] > 0))
                    swp_level += cvg_first[j] * df_pay[first_cpn_index[j]];
                else if ((k > 0) && (k < num_short_swp_cpn[j]))
                    swp_level +=
                        cvg_pay[first_cpn_index[j] + k - 1] * df_pay[first_cpn_index[j] + k];
                else if (k == num_short_swp_cpn[j])
                {
                    if (num_short_swp_cpn[j] > 0)
                        swp_level +=
                            cvg_pay[first_cpn_index[j] + k - 1] * df_pay[first_cpn_index[j] + k];
                    else if (num_short_swp_cpn[j] == 0)
                        swp_level += cvg_first[j] * df_pay[first_cpn_index[j] + k];
                }
            } /* END OF LOOP ON k */

            /* GET THE FORWARD */
            short_swps_fwd[j] =
                (df_start[j] - df_pay[first_cpn_index[j] + num_short_swp_cpn[j]]) / swp_level;
            yr_to_exp = (double)(ex_dates[j] - date_today) * YEARS_IN_DAY;

            /* GET THE SHORT SWPTIONS STRIKES */
            err = srt_f_get_vol(
                start_dates[j],
                start_dates[most_exp_index],
                short_swps_fwd[j],
                SRT_FALSE,
                &atm_vol);
            if (err)
                return err;

            err = srt_f_optsarbvol(
                short_swps_fwd[j],
                short_swps_fwd[j],
                yr_to_exp,
                atm_vol,
                0.0,
                0.0,
                0.0,
                vol_type,
                SRT_NORMAL,
                &norm_vol);
            if (err)
                return err;

            if (strike_method == EMK1)
                short_swps_strikes[j] = short_swps_fwd[j];
            else
                short_swps_strikes[j] =
                    short_swps_fwd[j] + long_swps_stdev[j] * norm_vol * sqrt(yr_to_exp);

            /* GET THE SHORT SWAPTIONS PRICES AND VEGAS */
            err = srt_f_get_vol(
                start_dates[j],
                start_dates[most_exp_index],
                short_swps_strikes[j],
                SRT_FALSE,
                &cev_vol);
            if (err)
                return err;

            if (vol_type == SRT_LOGNORMAL)
            {
                short_swps_pv[j] = 0.0;
                short_swps_pv[j] = srt_f_optblksch(
                    short_swps_fwd[j],
                    short_swps_strikes[j],
                    cev_vol,
                    yr_to_exp,
                    1.0,
                    SRT_PUT,
                    SRT_PREMIUM);
                short_swps_pv[j] *= swp_level;

                short_swps_vega[j] = 0;
                short_swps_vega[j] = srt_f_optblksch(
                    short_swps_fwd[j],
                    short_swps_strikes[j],
                    cev_vol + CEV_VOL_SHIFT,
                    yr_to_exp,
                    1.0,
                    SRT_PUT,
                    SRT_PREMIUM);
                short_swps_vega[j] *= swp_level;

                short_swps_vega[j] -= short_swps_pv[j];
                short_swps_vega[j] /= CEV_VOL_SHIFT;
            }

            else if (vol_type == SRT_NORMAL)
            {
                short_swps_pv[j] = 0.0;
                short_swps_pv[j] = srt_f_optblknrm(
                    short_swps_fwd[j],
                    short_swps_strikes[j],
                    cev_vol,
                    yr_to_exp,
                    1.0,
                    SRT_PUT,
                    SRT_PREMIUM);
                short_swps_pv[j] *= swp_level;

                short_swps_vega[j] = 0;
                short_swps_vega[j] = srt_f_optblknrm(
                    short_swps_fwd[j],
                    short_swps_strikes[j],
                    cev_vol + CEV_VOL_SHIFT,
                    yr_to_exp,
                    1.0,
                    SRT_PUT,
                    SRT_PREMIUM);
                short_swps_vega[j] *= swp_level;

                short_swps_vega[j] -= short_swps_pv[j];
                short_swps_vega[j] /= CEV_VOL_SHIFT;
            }

            else
                return serror("Expect a beta of 0 (NORMAL) or 1.0 (LOGNORMAL)");

            err = make_short_swaptions_payts(
                nEx,
                num_short_swp_cpn,
                first_cpn_index,
                cvg_first,
                cvg_pay,
                df_pay,
                short_swps_strikes,
                short_swps_payts);
            if (err)
                return err;

            if (lgm_ref_swp_data != NULL)
            {
                lgm_ref_swp_data->refShortSwptnArr[j - 1].bgnDt  = start_dates[j];
                lgm_ref_swp_data->refShortSwptnArr[j - 1].endDt  = start_dates[most_exp_index];
                lgm_ref_swp_data->refShortSwptnArr[j - 1].strike = short_swps_strikes[j];
                lgm_ref_swp_data->refShortSwptnArr[j - 1].bsPv   = short_swps_pv[j];
                lgm_ref_swp_data->refShortSwptnArr[j - 1].bsVol  = cev_vol;
            }

        } /* END OF LOOP ON j */

        hybrid_swpts_pv   = 0.0;
        hybrid_swpts_vega = 0.0;
        for (j = 1; j < nEx; j++)
        {
            hybrid_swpts_pv += short_swps_pv[j];
            hybrid_swpts_vega += short_swps_vega[j];
        }
    }

    for (i = 1; i <= nEx; i++)
    {
        if ((SParms.cal_method == FullCapAndDiag) || (SParms.cal_method == FixKappa))
        {
            err = srt_f_optimpvol(
                (short_swps_pv[i] / df_cap[i]),
                df_start[i] / df_cap[i],
                (1 + cvg_cap[i] * short_swps_strikes[i]),
                (double)(ex_dates[i] - date_today) * YEARS_IN_DAY,
                1.0,
                SRT_PUT,
                SRT_LOGNORMAL,
                &lgm_cap_vol[i]);

            if (err)
                return err;
        }
        else if (SParms.cal_method == HybridShortAndDiag)
        {
            err = srt_f_optimpvol(
                (caplets_pv[i] / df_cap[i]),
                df_start[i] / df_cap[i],
                (1 + cvg_cap[i] * caplets_strikes[i]),
                (double)(ex_dates[i] - date_today) * YEARS_IN_DAY,
                1.0,
                SRT_PUT,
                SRT_LOGNORMAL,
                &lgm_cap_vol[i]);

            if (err)
                return err;
        }
    }

    /* STORE THE INFORMATION IN THE STATIC PARAMS */
    SParms.date_today  = date_today;
    SParms.ex_dates    = ex_dates;
    SParms.start_dates = start_dates;

    SParms.swp_pay_dates      = swp_pay_dates;
    SParms.last_swp_pay_date  = swp_pay_dates[nPay - 1];
    SParms.cap_theo_end_date  = cap_theo_end_date;
    SParms.float_start_dates  = float_start_dates;
    SParms.float_end_dates    = float_end_dates;
    SParms.float_fixing_dates = float_fixing_dates;

    SParms.caplet_end_dates = caplet_end_dates;
    SParms.first_cpn_index  = first_cpn_index;
    SParms.df_start         = df_start;
    SParms.df_cap           = df_cap;
    SParms.df_pay           = df_pay;
    SParms.cvg_pay          = cvg_pay;
    SParms.cvg_cap          = cvg_cap;
    SParms.cvg_first        = cvg_first;
    SParms.cvg_first        = cvg_first;

    SParms.long_swps_strikes = long_swps_strikes;
    SParms.long_swps_payts   = long_swps_payts;
    SParms.long_swps_fwd     = long_swps_fwd;
    SParms.long_swps_cev     = long_swps_cev;
    SParms.long_swps_beta    = long_swps_beta;
    SParms.long_swps_pv      = long_swps_pv;
    SParms.long_swps_vega    = long_swps_vega;
    SParms.long_swps_stdev   = long_swps_stdev;

    SParms.short_swps_strikes = short_swps_strikes;
    SParms.short_swps_payts   = short_swps_payts;
    SParms.short_swps_fwd     = short_swps_fwd;
    SParms.short_swps_cev     = short_swps_cev;
    SParms.short_swps_beta    = short_swps_beta;
    SParms.short_swps_pv      = short_swps_pv;
    SParms.short_swps_vega    = short_swps_vega;
    SParms.num_short_swp_cpn  = num_short_swp_cpn;
    SParms.caplets_pv         = caplets_pv;
    SParms.caplets_vega       = caplets_vega;
    SParms.caplets_strikes    = caplets_strikes;

    SParms.hybrid_swpts_pv   = hybrid_swpts_pv;
    SParms.hybrid_swpts_vega = hybrid_swpts_vega;

    SParms.float_coverages                    = float_coverages;
    SParms.float_start_dates_discount_factors = float_start_dates_discount_factors;
    SParms.float_end_dates_discount_factors   = float_end_dates_discount_factors;
    SParms.pv_full_cap                        = cap_pv;
    SParms.vega_full_cap                      = cap_vega;

    SParms.nPay            = nPay;
    SParms.nEx             = nEx;
    SParms.num_float_dates = num_float_dates;

    SParms.float_basis = float_basis;
    SParms.float_freq  = float_freq;

    SParms.lgm_cap_vol = lgm_cap_vol;

    return err;
}

static Err newtonfuncforcalbootonfullcap(double lambda, double* f)
{
    Err    err = NULL;
    double theo_hybrid_swpts_pv, theo_cap_pv, mkt_cap_pv = SParms.pv_full_cap,
                                              mkt_cap_vega          = SParms.vega_full_cap,
                                              mkt_hybrid_swpts_pv   = SParms.hybrid_swpts_pv,
                                              mkt_hybrid_swpts_vega = SParms.hybrid_swpts_vega;
    LGMCalMeth cal_method                                           = SParms.cal_method;

    s_optim_lambda[1] = lambda;

    err = srt_f_autocal_calboot(s_optim_lambda, cal_method);
    if (err)
        return err;

    (*f) = 0.;
    if (cal_method == FullCapAndDiag)
    {
        err = srt_f_autocal_cap_greeks(&theo_cap_pv);
        if (err)
            return err;

        (*f) += (theo_cap_pv - mkt_cap_pv);
    }
    else if (cal_method == HybridShortAndDiag)
    {
        err = srt_f_autocal_hybrid_swpts_greeks(&theo_hybrid_swpts_pv);
        if (err)
            return err;

        (*f) += (theo_hybrid_swpts_pv - mkt_hybrid_swpts_pv);
    }

    return err;
}

Err lgm_2f_autocal_calibrate(
    void*               dealPtr,
    LGM_TS**            lgm_ts_ptr_ptr,
    SrtLgmRefSwptnData* lgm_ref_swp_data,
    LGMCalParm*         cal_req,
    Err (*srt_f_get_vol)(Date, Date, double, SRT_Boolean, double*),
    Err (*srt_f_get_beta)(Date, Date, double*),
    String yc_name)
{
    Err          err   = NULL;
    SrtSimMidAt* MidAt = (SrtSimMidAt*)dealPtr;
    long         i, j, k, nEx = MidAt->nEx, nPay = MidAt->nPay;
    LGM_TS*      tsPtr = LGMCreateLGM2F_TS(nEx + 1, nEx + 1);
    Date         date_today, *ex_dates, *start_dates, *swp_pay_dates = MidAt->tPay;

    double *H1 = (tsPtr->H1), *H2 = (tsPtr->H2), *Zeta1 = (tsPtr->Zeta1), *Zeta2 = (tsPtr->Zeta2),
           *Zeta12 = (tsPtr->Zeta12), alpha = cal_req->alpha, gamma = cal_req->gamma,
           rho = cal_req->rho, *lambdaPtr;

    double *    gi, *Wi, *LongKs;
    double      nstop, a[3], b[3], kappa;
    long        lStartTime = 0, lEndTime = 0;
    SrtCurvePtr yld_crv;

    long* likely_most_expensive_index;

    double *start_zeta1s = cal_req->StartZeta1s, *start_taus = cal_req->StartTaus,
           **HybridShortInstrsIndex = cal_req->HybridShortInstrsIndex;

    SRT_Boolean swith_to_fix_kappa_calibration = SRT_NO;
    double *    kappa_iteration_values         = NULL, two_last_iter_diff_tau, last_iter_diff_tau;

    yld_crv    = lookup_curve(yc_name);
    date_today = get_clcndate_from_yldcrv(yld_crv);

    start_dates = lngvector(1, nEx);
    for (i = 1; i <= nEx; i++)
        start_dates[i] = (MidAt->tStart)[i - 1];

    gi                          = dvector(1, NBR_HERMITE_POINTS);
    Wi                          = dvector(1, NBR_HERMITE_POINTS);
    lambdaPtr                   = dvector(1, 1);
    LongKs                      = dvector(1, nEx);
    likely_most_expensive_index = lngvector(1, nEx);

    err = make_static_allocations(nEx, nPay);
    if (err)
        return err;

    err = update_autocal_midat_struct(MidAt, yc_name, srt_f_get_vol, srt_f_get_beta);
    if (err)
        return err;

    err = HermiteStandard(gi, Wi, NBR_HERMITE_POINTS);
    if (err)
        return err;

    for (i = 1; i <= nEx; i++)
        LongKs[i] = (MidAt->MidAtStrike[i - 1]);

    if (cal_req->calmeth == HybridShortAndDiag)
    {
        for (i = 1; i <= nEx; i++)
            likely_most_expensive_index[i] = (long)HybridShortInstrsIndex[i - 1][1];
    }

    SParms.cal_method   = cal_req->calmeth;
    SParms.start_zeta1s = start_zeta1s;
    SParms.start_taus   = start_taus;
    if (SParms.start_zeta1s == NULL)
        SParms.use_start_point_for_zeta1s = SRT_NO;
    else
        SParms.use_start_point_for_zeta1s = SRT_YES;

    if (SParms.start_taus == NULL)
        SParms.use_start_point_for_taus = SRT_NO;
    else
        SParms.use_start_point_for_taus = SRT_YES;

    err = srt_f_get_autocal_ref_instrs(
        nEx,
        nPay,
        start_dates,
        swp_pay_dates,
        yc_name,
        srt_f_get_vol,
        srt_f_get_beta,
        LongKs,
        (cal_req->Rmeth), /* TO DETERMINE THE SHORT SWAPTIONS INSTRUMENT  */
        (cal_req->calmeth),
        likely_most_expensive_index,
        lgm_ref_swp_data);

    if (err)
        return err;

    err = make_levenberg_struct(nEx);
    if (err)
        return err;

    /* STORE SOME INFO IN THE STATIC STRUCTURE  */
    SParms.gi     = gi;
    SParms.Wi     = Wi;
    SParms.Zeta1  = Zeta1;
    SParms.Zeta2  = Zeta2;
    SParms.Zeta12 = Zeta12;
    SParms.H1     = H1;
    SParms.H2     = H2;
    SParms.alpha  = alpha;
    SParms.gamma  = gamma;
    SParms.rho    = rho;

    /* GET INFO FROM STATIC STRUCTURE */
    ex_dates   = SParms.ex_dates;
    date_today = SParms.date_today;

    for (i = 1; i <= nEx; i++)
        s_yr_to_ex_dates_ptr[i] = (double)(ex_dates[i] - date_today) * YEARS_IN_DAY;
    for (i = 1; i <= nEx; i++)
        s_yr_to_start_ptr[i] = (double)(start_dates[i] - date_today) * YEARS_IN_DAY;
    s_yr_to_last_pay_date = (double)(swp_pay_dates[nPay - 1] - date_today) * YEARS_IN_DAY;

    /* START OF THE CALIBRATION CORE */
    if (nEx == 1)
    {
        SParms.apply_criteria = SRT_TRUE;
        if (SParms.use_start_point_for_taus == SRT_YES)
            lambdaPtr[1] = 1.0 / SParms.start_taus[0];
        else
            lambdaPtr[1] = 1.0 / FIXED_TAU;

        lStartTime = time(NULL);
        err        = srt_f_autocal_calboot(lambdaPtr, FixKappa);
        if (err)
            return err;

        kappa = lambdaPtr[1];

        lEndTime = time(NULL);
        smessage("Calibration achieved in %.2f s", difftime(lEndTime, lStartTime));
    }
    else
    {
        switch (cal_req->calmeth)
        {
        case FixKappa:

            lambdaPtr[1]          = cal_req->kap;
            SParms.apply_criteria = SRT_TRUE;

            lStartTime = time(NULL);

            err = srt_f_autocal_calboot(lambdaPtr, FixKappa);
            if (err)
                return err;

            kappa    = lambdaPtr[1];
            lEndTime = time(NULL);
            smessage("Calibration achieved in %.2f s", difftime(lEndTime, lStartTime));

            break;
        case FullCapAndDiag:
        case HybridShortAndDiag:

            lStartTime = time(NULL);

            /* NEWTON ROUTINE TO CALIBRATE THE MEANREVERSION */
            nstop = 0.0;
            if (SParms.use_start_point_for_taus == SRT_NO)
                a[0] = 1.0 / FIXED_TAU;
            else
                a[0] = 1.0 / SParms.start_taus[0];

            err = newtonfuncforcalbootonfullcap(a[0], &b[0]);
            if (err)
                return err;

            if (SParms.use_start_point_for_taus == SRT_NO)
                a[1] = 1.0 / (FIXED_TAU + 1.0);
            else
                a[1] = 1.0 / (SParms.start_taus[0] + 1.0);

            err = newtonfuncforcalbootonfullcap(a[1], &b[1]);
            if (err)
                return err;

            if (SParms.use_start_point_for_taus == SRT_NO)
                a[2] = 1.0 / (FIXED_TAU + 2.0);
            else
                a[2] = 1.0 / (SParms.start_taus[0] + 2);
            SParms.apply_criteria = SRT_FALSE;

            kappa_iteration_values = dvector(1, MAX_ITERATION);

            k = 0;
            while ((nstop < 1.0) && (k < MAX_ITERATION) &&
                   (swith_to_fix_kappa_calibration == SRT_NO))
            {
                err = newtonfuncforcalbootonfullcap(a[2], &b[2]);
                if (err)
                    return err;

                newton(0.0, 2.0, a, b, &nstop);

                k++;

                if ((a[2] < MIN_MEAN_REVERSION) || (a[2] > MAX_MEAN_REVERSION))
                    swith_to_fix_kappa_calibration = SRT_YES;

                if (fabs(b[2]) < WEAK_CAP_VOL_TOL)
                    SParms.apply_criteria = SRT_TRUE;
                if ((fabs(b[2]) < CAP_VOL_TOL) && (SParms.apply_criteria == SRT_TRUE))
                    nstop = 1.0;

                kappa_iteration_values[k] = a[2]; /* store the kappa values */

                if ((k > 2) && (nstop < 1.0) && (k < MAX_ITERATION) &&
                    (swith_to_fix_kappa_calibration == SRT_NO))
                {
                    last_iter_diff_tau =
                        1.0 / kappa_iteration_values[k] - 1.0 / kappa_iteration_values[k - 1];
                    two_last_iter_diff_tau =
                        1.0 / kappa_iteration_values[k] - 1.0 / kappa_iteration_values[k - 2];

                    if ((fabs(last_iter_diff_tau) <= MIN_ITER_TAU_DIFF) &&
                        (fabs(two_last_iter_diff_tau) <= MIN_ITER_TAU_DIFF))
                        swith_to_fix_kappa_calibration = SRT_YES;
                }
            }
            /* END OF CALIBRATION OF THE TAU - THE OPTIMAL TAU IS 1/a[2] */

            lEndTime = time(NULL);
            smessage("Calibration achieved in %.2f s", difftime(lEndTime, lStartTime));

            kappa = a[2];

            if (swith_to_fix_kappa_calibration == SRT_YES)
            {
                if (SParms.use_start_point_for_taus == SRT_YES)
                    kappa = 1.0 / SParms.start_taus[0];
                else
                    kappa = 1.0 / FIXED_TAU;

                lambdaPtr[1] = kappa;
                err          = srt_f_autocal_calboot(lambdaPtr, FixKappa);
                if (err)
                    return err;
            }

            for (j = nEx; j > 0; j--)
                H1[j - 1] =
                    (exp(-s_yr_to_start_ptr[j] * kappa) - exp(-s_yr_to_last_pay_date * kappa)) /
                    kappa;
            for (j = 1; j <= nEx; j++)
                H2[j - 1] = H2_Func(start_dates[j], start_dates, date_today, nEx, H1, gamma);

            /* COMPUTE ZETA2, ZETA12 AND H2 */
            for (j = 1; j <= nEx; j++)
            {
                Zeta2[j - 1] = New_Zeta2_Func(
                    ex_dates[j], ex_dates, nEx, date_today, Zeta1, Zeta2, alpha, gamma, rho);
                Zeta12[j - 1] = New_Zeta12_Func(
                    ex_dates[j], ex_dates, nEx, date_today, Zeta1, Zeta12, alpha, gamma, rho);
            }

            if (kappa_iteration_values)
                free_dvector(kappa_iteration_values, 1, MAX_ITERATION);
            kappa_iteration_values = NULL;

            break;
        }
    }

    (tsPtr->alpha)     = alpha;
    (tsPtr->gamma)     = gamma;
    (tsPtr->one_kappa) = kappa;
    (tsPtr->rho)       = rho;
    (*lgm_ts_ptr_ptr)  = tsPtr;

    err = make_static_desallocations(nEx, nPay);
    if (err)
        return err;

    err = make_static_struct_desallocations(nEx, nPay);
    if (err)
        return err;

    err = delete_levenberg_struct(nEx);
    if (err)
        return err;

    if (lambdaPtr)
        free_dvector(lambdaPtr, 1, 1);
    lambdaPtr = NULL;
    if (gi)
        free_dvector(gi, 1, NBR_HERMITE_POINTS);
    gi = NULL;
    if (Wi)
        free_dvector(Wi, 1, NBR_HERMITE_POINTS);
    Wi = NULL;
    if (likely_most_expensive_index)
        free_lngvector(likely_most_expensive_index, 1, nEx);
    likely_most_expensive_index = NULL;
    if (start_dates)
        free_lngvector(start_dates, 1, nEx);
    start_dates = NULL;

    return (err);
}

#undef MAX_ITERATION
#undef FIXED_TAU
#undef MIN_MEAN_REVERSION
#undef MAX_MEAN_REVERSION
#undef MIN_ITER_TAU_DIFF
#undef NBR_HERMITE_POINTS
#undef MAX_STDEV
#undef SIGMA_MIN_SQ
#undef FIRST_EX_SIGMA_MIN_SQ
#undef SIGMA_MAX_SQ
#undef ONE_HALF
#undef WEAK_CAP_VOL_TOL
#undef DAYS_IN_QUARTER
