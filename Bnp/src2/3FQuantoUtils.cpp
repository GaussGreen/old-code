
/* ==========================================================================
   FILE_NAME:	3FQuantoUtils.c

   PURPOSE:		Useful functions for the 3FQuanto model

   DATE:		26 Sep 2003

   AUTHOR:		J.M.L.
   ========================================================================== */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_all3FQuanto.h"
#include "srt_h_allFx3F.h"

// Selects the dates in r1_dates, r3_dates, fx_dates which are <= last_time
// and put them in a single array. (r1 then r3 then fx)
Err compute_vol_times_3FQuanto(
    double*  r1_dates,
    double*  r1_sigma,
    int      n_r1_dates,
    double*  r3_dates,
    double*  r3_sigma,
    int      n_r3_dates,
    double*  fx_dates,
    double*  fx_sigma,
    int      n_fx_dates,
    int*     num_vol_times,
    double** vol_times,
    double   last_time)
{
    Err err = NULL;
    int i;

    *vol_times = (double*)calloc(n_r1_dates + n_r3_dates + n_fx_dates, sizeof(double));
    if (!vol_times)
    {
        err = "Memory allocation error in compute_vol_times_3FQuanto";
        goto FREE_RETURN;
    }

    *num_vol_times = 0;

    i = 1;
    while ((i < n_r1_dates) && (r1_dates[i - 1] < last_time))
    {
        if (fabs(r1_sigma[i] - r1_sigma[i - 1]) > EPS)
        {
            (*vol_times)[*num_vol_times] = r1_dates[i - 1];
            (*num_vol_times)++;
        }
        i++;
    }

    i = 1;
    while ((i < n_r3_dates) && (r3_dates[i - 1] < last_time))
    {
        if (fabs(r3_sigma[i] - r3_sigma[i - 1]) > EPS)
        {
            (*vol_times)[*num_vol_times] = r3_dates[i - 1];
            (*num_vol_times)++;
        }
        i++;
    }

    i = 1;
    while ((i < n_fx_dates) && (fx_dates[i - 1] < last_time))
    {
        if (fabs(fx_sigma[i] - fx_sigma[i - 1]) > EPS)
        {
            (*vol_times)[*num_vol_times] = fx_dates[i - 1];
            (*num_vol_times)++;
        }
        i++;
    }

    if (!(*num_vol_times))
    {
        free(*vol_times);
        *vol_times = NULL;
    }

FREE_RETURN:

    return err;
}

// Merge rates, fx and corr term structures and returns all variables in lists
Err merge_rates_fx_corr_ts_3FQuanto(
    double*   sig_r1_mat,
    double*   sig_r1,
    int       sig_n_r1,
    double*   sig_r3_mat,
    double*   sig_r3,
    int       sig_n_r3,
    double*   sig_fx_mat,
    double*   sig_fx,
    int       sig_n_fx,
    double*   corr_4x4_mat,
    double*** corr_4x4,
    int       corr_4x4_n_mat,

    // Outputs
    double** sig_merge_mat,
    double** sig_merge_r1,
    double** sig_merge_r3,
    double** sig_merge_fx,
    double** corr_merge_r1_r2,
    double** corr_merge_r1_r3,
    double** corr_merge_r1_fx,
    double** corr_merge_r2_r3,
    double** corr_merge_r2_fx,
    double** corr_merge_r3_fx,
    int*     sig_merge_n)
{
    int i, index;
    Err err = NULL;

    *sig_merge_mat    = NULL;
    *sig_merge_r1     = NULL;
    *sig_merge_r3     = NULL;
    *sig_merge_fx     = NULL;
    *corr_merge_r1_r2 = NULL;
    *corr_merge_r1_r3 = NULL;
    *corr_merge_r1_fx = NULL;
    *corr_merge_r2_r3 = NULL;
    *corr_merge_r2_fx = NULL;
    *corr_merge_r3_fx = NULL;

    *sig_merge_mat = (double*)calloc(sig_n_r1, sizeof(double));
    if (!(*sig_merge_mat))
    {
        err = "Memory allocation error in merge_rates_fx_corr_ts_3FQuanto";
        goto FREE_RETURN;
    }

    memcpy(*sig_merge_mat, sig_r1_mat, sig_n_r1 * sizeof(double));
    *sig_merge_n = sig_n_r1;
    num_f_concat_vector(sig_merge_n, sig_merge_mat, sig_n_r3, sig_r3_mat);
    num_f_concat_vector(sig_merge_n, sig_merge_mat, sig_n_fx, sig_fx_mat);
    num_f_concat_vector(sig_merge_n, sig_merge_mat, corr_4x4_n_mat, corr_4x4_mat);
    num_f_sort_vector(*sig_merge_n, *sig_merge_mat);
    num_f_unique_vector(sig_merge_n, *sig_merge_mat);

    *sig_merge_r1     = (double*)calloc(*sig_merge_n, sizeof(double));
    *sig_merge_r3     = (double*)calloc(*sig_merge_n, sizeof(double));
    *sig_merge_fx     = (double*)calloc(*sig_merge_n, sizeof(double));
    *corr_merge_r1_r2 = (double*)calloc(*sig_merge_n, sizeof(double));
    *corr_merge_r1_r3 = (double*)calloc(*sig_merge_n, sizeof(double));
    *corr_merge_r1_fx = (double*)calloc(*sig_merge_n, sizeof(double));
    *corr_merge_r2_r3 = (double*)calloc(*sig_merge_n, sizeof(double));
    *corr_merge_r2_fx = (double*)calloc(*sig_merge_n, sizeof(double));
    *corr_merge_r3_fx = (double*)calloc(*sig_merge_n, sizeof(double));

    if (!(*sig_merge_r1) || !(*sig_merge_r3) || !(*sig_merge_fx) || !(*corr_merge_r1_r2) ||
        !(*corr_merge_r1_r3) || !(*corr_merge_r1_fx) || !(*corr_merge_r2_r3) ||
        !(*corr_merge_r2_fx) || !(*corr_merge_r3_fx))
    {
        err = "Memory allocation error (2) in merge_rates_fx_corr_ts_3FQuanto";
        goto FREE_RETURN;
    }

    for (i = *sig_merge_n - 1; i >= 0; i--)
    {
        (*sig_merge_r1)[i] = sig_r1[Get_Index((*sig_merge_mat)[i], sig_r1_mat, sig_n_r1)];
        (*sig_merge_r3)[i] = sig_r3[Get_Index((*sig_merge_mat)[i], sig_r3_mat, sig_n_r3)];
        (*sig_merge_fx)[i] = sig_fx[Get_Index((*sig_merge_mat)[i], sig_fx_mat, sig_n_fx)];

        index                  = Get_Index((*sig_merge_mat)[i], corr_4x4_mat, corr_4x4_n_mat);
        (*corr_merge_r1_r2)[i] = corr_4x4[index][LGM2F_W1][LGM2F_W2];
        (*corr_merge_r1_r3)[i] = corr_4x4[index][LGM2F_W1][LGM1F_W3];
        (*corr_merge_r1_fx)[i] = corr_4x4[index][LGM2F_W1][FX_3F_QUANTO];
        (*corr_merge_r2_r3)[i] = corr_4x4[index][LGM2F_W2][LGM1F_W3];
        (*corr_merge_r2_fx)[i] = corr_4x4[index][LGM2F_W2][FX_3F_QUANTO];
        (*corr_merge_r3_fx)[i] = corr_4x4[index][LGM1F_W3][FX_3F_QUANTO];
    }

FREE_RETURN:

    if (err)
    {
        if (*sig_merge_mat)
            free(*sig_merge_mat);
        if (*sig_merge_r1)
            free(*sig_merge_r1);
        if (*sig_merge_r3)
            free(*sig_merge_r3);
        if (*sig_merge_fx)
            free(*sig_merge_fx);
        if (*corr_merge_r1_r2)
            free(*corr_merge_r1_r2);
        if (*corr_merge_r1_r3)
            free(*corr_merge_r1_r3);
        if (*corr_merge_r1_fx)
            free(*corr_merge_r1_fx);
        if (*corr_merge_r2_r3)
            free(*corr_merge_r2_r3);
        if (*corr_merge_r2_fx)
            free(*corr_merge_r2_fx);
        if (*corr_merge_r3_fx)
            free(*corr_merge_r3_fx);
    }
    return err;
}

// Fills the expectations (under Q_Beta_Domestic) and variances for a 3FQuanto model
void fill_fwd_var_corr_3FQuanto(
    long    nstp,
    double* time,
    double* date,
    double* sigma_R1,
    double* sigma_R2,
    double* sigma_R3,
    double* sigma_FX,
    double  lambda_R1,
    double  lambda_R2,
    double  lambda_R3,
    double* correl_R1_R2,
    double* correl_R1_R3,
    double* correl_R1_FX,
    double* correl_R2_R3,
    double* correl_R2_FX,
    double* correl_R3_FX,
    char*   lgm2F_yc,
    char*   lgm1F_yc,
    int     lgm2F_is_dom_for,
    /*	To be allocated by caller, filled by the function */
    double* fwd_R1,
    double* var_R1,
    double* fwd_R2,
    double* var_R2,
    double* covar_R1_R2,
    double* fwd_R3,
    double* var_R3,
    double* lgm2F_ifr,
    double* lgm1F_ifr)
{
    int    i, j;
    double vol_r1, vol_r2, vol_r3, vol_fx;
    double var_r1, var_r2, covar_r1_r2, var_r3;
    double fwd_r1, fwd_r2, fwd_r3;
    double correl_r1_r2, correl_r1_r3, correl_r1_fx, correl_r2_r3, correl_r2_fx, correl_r3_fx;
    double t, t1, t2, lambda_R1R2;
    double aux1_1, aux1_2, aux2_1, aux2_2, aux3_1, aux3_2, L1, L2, L3, P;

    lgm2F_ifr[0] = swp_f_zr(date[0], date[1], lgm2F_yc);
    lgm1F_ifr[0] = swp_f_zr(date[0], date[1], lgm1F_yc);

    lambda_R1R2 = lambda_R1 + lambda_R2;

    for (i = 1; i < nstp; i++)
    {
        if (i < nstp - 1)
        {
            lgm2F_ifr[i] = swp_f_zr(date[i], date[i + 1], lgm2F_yc);
            lgm1F_ifr[i] = swp_f_zr(date[i], date[i + 1], lgm1F_yc);
        }
        else
        {
            lgm2F_ifr[i] = swp_f_zr(date[i], date[i] + 1, lgm2F_yc);
            lgm1F_ifr[i] = swp_f_zr(date[i], date[i] + 1, lgm1F_yc);
        }

        t           = time[i];
        fwd_r1      = 0.0;
        var_r1      = 0.0;
        fwd_r2      = 0.0;
        var_r2      = 0.0;
        covar_r1_r2 = 0.0;
        fwd_r3      = 0.0;
        var_r3      = 0.0;
        t1          = 0.0;

        /* integrates accross all Tj in [0, i] */
        for (j = 0; j < i; j++)
        {
            t2 = time[j + 1];

            /* on ]j, j+1], we use the volatility at time j+1 */
            vol_r1       = sigma_R1[j + 1];
            vol_r2       = sigma_R2[j + 1];
            vol_r3       = sigma_R3[j + 1];
            vol_fx       = sigma_FX[j + 1];
            correl_r1_r2 = correl_R1_R2[j + 1];
            correl_r1_r3 = correl_R1_R3[j + 1];
            correl_r1_fx = correl_R1_FX[j + 1];
            correl_r2_r3 = correl_R2_R3[j + 1];
            correl_r2_fx = correl_R2_FX[j + 1];
            correl_r3_fx = correl_R3_FX[j + 1];
            aux1_1       = exp(-lambda_R1 * (t - t1));
            aux1_2       = exp(-lambda_R1 * (t - t2));
            aux2_1       = exp(-lambda_R2 * (t - t1));
            aux2_2       = exp(-lambda_R2 * (t - t2));
            aux3_1       = exp(-lambda_R3 * (t - t1));
            aux3_2       = exp(-lambda_R3 * (t - t2));

            /* R1 non quanto fwd */
            L1 = 1.0 / (lambda_R1 * lambda_R1) *
                 (aux1_2 * (1.0 - 0.5 * aux1_2) - aux1_1 * (1.0 - 0.5 * aux1_1));
            L2 = 1.0 / lambda_R2 *
                 (aux1_2 * (1.0 / lambda_R1 - aux2_2 / lambda_R1R2) -
                  aux1_1 * (1.0 / lambda_R1 - aux2_1 / lambda_R1R2));
            fwd_r1 += vol_r1 * vol_r1 * L1 + correl_r1_r2 * vol_r1 * vol_r2 * L2;

            /* R2 non quanto fwd */
            L1 = 1.0 / (lambda_R2 * lambda_R2) *
                 (aux2_2 * (1.0 - 0.5 * aux2_2) - aux2_1 * (1.0 - 0.5 * aux2_1));
            L2 = 1.0 / lambda_R1 *
                 (aux2_2 * (1.0 / lambda_R2 - aux1_2 / lambda_R1R2) -
                  aux2_1 * (1.0 / lambda_R2 - aux1_1 / lambda_R1R2));
            fwd_r2 += vol_r2 * vol_r2 * L1 + correl_r1_r2 * vol_r1 * vol_r2 * L2;

            /* R3 non quanto fwd */
            L1 = 1.0 / (lambda_R3 * lambda_R3) *
                 (aux3_2 * (1.0 - 0.5 * aux3_2) - aux3_1 * (1.0 - 0.5 * aux3_1));
            fwd_r3 += vol_r3 * vol_r3 * L1;

            /* R1 cumulated var */
            P = 1.0 / (2.0 * lambda_R1) * (aux1_2 * aux1_2 - aux1_1 * aux1_1);
            var_r1 += vol_r1 * vol_r1 * P;

            /* R2 cumulated var */
            P = 1.0 / (2.0 * lambda_R2) * (aux2_2 * aux2_2 - aux2_1 * aux2_1);
            var_r2 += vol_r2 * vol_r2 * P;

            /* R1_R2 cumulated covar */
            P = 1.0 / (lambda_R1R2) * (aux2_2 * aux1_2 - aux2_1 * aux1_1);
            covar_r1_r2 += vol_r1 * vol_r2 * P;

            /* R3 cumulated var */
            P = 1.0 / (2.0 * lambda_R3) * (aux3_2 * aux3_2 - aux3_1 * aux3_1);
            var_r3 += vol_r3 * vol_r3 * P;

            /* Quanto correction */
            if (lgm2F_is_dom_for == 1)
            {
                L3 = 1.0 / lambda_R1 * (aux1_2 - aux1_1);
                fwd_r1 += -correl_r1_fx * vol_r1 * vol_fx * L3;
                L3 = 1.0 / lambda_R2 * (aux2_2 - aux2_1);
                fwd_r2 += -correl_r2_fx * vol_r2 * vol_fx * L3;
            }
            else
            {
                L3 = 1.0 / lambda_R3 * (aux3_2 - aux3_1);
                fwd_r3 += -correl_r3_fx * vol_r3 * vol_fx * L3;
            }

            t1 = t2;
        }

        /* Writes the results into the term structures */
        fwd_R1[i] = fwd_r1;
        var_R1[i] = var_r1;

        fwd_R2[i] = fwd_r2;
        var_R2[i] = var_r2;

        covar_R1_R2[i] = covar_r1_r2;

        fwd_R3[i] = fwd_r3;
        var_R3[i] = var_r3;
    }
}